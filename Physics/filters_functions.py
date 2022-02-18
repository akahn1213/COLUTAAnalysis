#!/usr/bin/env python2
# -*- coding: utf-8 -*-

"""Module to wrap code for the CR-RC² transfer function."""

from __future__ import absolute_import
from __future__ import print_function
from __future__ import division

import numpy
from scipy import linalg

#k's imports
import pylab as pl 
import sys 
import matplotlib
import matplotlib.pyplot as plt
import math
from decimal import *
from scipy.stats import * #norm

class FirFilter(object):

    """Simple implementation of a finite-impulse-response filter.

    From a possibly infinite input sequence `x_n` and its `N`
    coefficients `c_n`, an FIR filter calculates its output sequence
    `y_n` as follows:

        y_n = sum(i=0, N-1) c_i x_{n-i},

    i.e. the output of an FIR filter is the convolution of its input
    with its coefficients.

    This class allows applying the filter on a sample-by-sample basis
    (by keeping a history of the last `N` passed samples) and
    batch-applying it.
    """

    def __init__(self, coeffs):
        """Creates a new FIR filter and sets its coefficients."""
        self.coeffs = numpy.ravel(coeffs)
        self.history = numpy.zeros_like(coeffs)
        print(coeffs)

    def __len__(self):
        """Returns the number of coefficients of this filter."""
        return len(self.coeffs)

    def clear(self):
        """Set all elements of the input history to zero."""
        self.history[:] = 0.0

    def update(self, sample):
        """Update the input history and calculate the next output."""
        self.history = numpy.roll(self.history, 1)
        self.history[0] = sample
        return numpy.dot(self.coeffs, self.history)

    __call__ = update

    def process(self, samples):
        """Like `update()`, but applied to an entire array of samples."""
        # Take all but the oldest history sample and reverse them so
        # that the newest sample comes last.
        history = self.history[-2::-1]
        # Prepend them to the given samples.
        samples = numpy.hstack((history, numpy.ravel(samples)))
        # Update the history.
        self.history = samples[-len(self):]
        # Now mode='valid' gives us exactly one output for each input.
        return numpy.convolve(self.coeffs, samples, mode='valid')


def calc_of_coeffs(ac, g, dg=None):
    # pylint: disable=invalid-name,too-many-locals
    """Calculate Optimal Filter coefficients.

    The [Optimal Filter][1] (OF) is an FIR filter that minimizes the
    mean square error between a sample sequence `S_n` and the
    amplitudes `A_n` and arrival times `tau_n` that caused this
    sequence. This error is weighted by `R_mn`, the autocorrelation of
    a noise sequence smearing out the samples.

    Convolving the filter coefficients `a_n` with an input sequence
    `x_n` gives, for each sample, an estimate of the signal amplitude
    `A_n` contributing this sample. Convolving the coefficients `b_n`
    with `x_n` gives `A_n*tau_n`, where `tau_n` is an estimate for the
    arrival time of the hit `A_n`:

        A_n = sum(i=0, D-1) a_i x_{n-i},
        A_n*tau_n = sum(i=0, D-1) b_i x_{n-i}.

    Note that this function returns the filter coefficients in reverse
    order to the original paper. The paper defines them in such an
    order that each output sample y_n is calculated from the input
    samples x_n and the filter coefficients a_n as follows:

        y_n = sum(i=0, D-1) a_i x_{n-D+i+1},

    i.e. the last coefficient is applied to the "current" sample. In
    contrast, the former definition applies the *first* coefficient to
    the current sample. This has the advantage that the definition is
    the same as for the convolution operation.

    Args:
        ac: 1D array of the noise's autocorrelation coefficients. They
            are normalized in such a way that `ac[0]==1.0`.
        g: 1D array of the ideal pulse shape samples. Must have the
            same length as `ac`.
        dg: If passed, a 1D array of the ideal pulse shape's
            derivative. If passed, it must have the same length as
            `ac`. If not passed, the derivative is calculated using
            `numpy.gradient`.

    Returns:
        A tuple `(a_coeffs, b_coeffs)` of 1D arrays with the same
        length as `pulse`.

    [1]: Cleland, Stern: Signal processing considerations for liquid
         ionization calorimeters in high rate environment

    """
    ac = numpy.ravel(ac) / ac[0]
    g = numpy.ravel(g)
    dg = numpy.ravel(dg) if dg is not None else numpy.gradient(g)
    scale = max(g)  
    g = [i/float(scale) for i in g] #g /= scale does not work?!?!
    dg /= scale
    # Calculate V = R^{-1}.
    inv_ac = linalg.inv(linalg.toeplitz(ac))
    # Calculate V*g and V*dg only once.
    vg = numpy.dot(inv_ac, g)
    vdg = numpy.dot(inv_ac, dg)
    # Calculate helper variables.
    q1 = numpy.dot(g, vg)
    q2 = numpy.dot(dg, vdg)
    q3 = numpy.dot(dg, vg)
    delta = q1*q2 - q3*q3
    # Calculate Lagrange multipliers
    lm_lambda = q2/delta
    lm_kappa = -q3/delta
    lm_mu = q3/delta
    lm_rho = -q1/delta
    # Calculate filter coefficients.
    a_coeffs = lm_lambda*vg + lm_kappa*vdg
    b_coeffs = lm_mu*vg + lm_rho*vdg
    # DON'T Reverse order to get canonical coefficient order.
    return a_coeffs[::1], b_coeffs[::1] #a_coeffs, b_coeffs # 

def calc_wiener_coeffs(pulse, length, offset=1, prepeak=False, postpeak=False):
    """Calculate Wiener filter coefficients.

    The [Wiener filter][1] as used in this module is an FIR filter that
    minimizes the mean square error between a measured pulse and a
    given ideal pulse. It does so by effectively calculating
    coefficients for a deconvolution matrix and applying this matrix to
    its input. By assuming a delta peak as ideal pulse, the Wiener
    filter can deconvolute a sample sequence from the known pulse
    shape.

    The ideal pulse need not be an exact delta peak. It may also be
    slightly widened by adding a pre- or a post-peak sample (or both)
    of half height. Widening the ideal pulse may improve performance in
    the face of noise.

    Args:
        pulse: 1D array of the ideal pulse shape samples. Must have a
            length greater than `length`.
        length: The length of the Wiener filter. A longer filter has
            more coefficients and can determine the deconvolution
            matrix more accurately.
        offset: An offset between `pulse`'s peak and the delta peak.
            A positive value gives the filter more information to base
            its estimate on, a negative value makes the filter
            predictive.
        prepeak: If passed and `True`, the ideal pulse is given a
            pre-peak sample of height 0.5.
        prepeak: If passed and `True`, the ideal pulse is given a
            post-peak sample of height 0.5.

    Returns:
        A 1D array of with length `length` containing the filter
        coefficients.

    [1]: https://en.wikipedia.org/wiki/Wiener_filter
    """
    pulse = numpy.ravel(pulse)
    i_max = numpy.argmax(pulse)
    assert 0 < i_max < len(pulse)-1, 'maximum must not lie on edge'

    target = numpy.zeros_like(pulse)
    target[i_max+offset] = 1.0
    target[i_max+offset-1] = 0.5 if prepeak else -0.0
    target[i_max+offset+1] = 0.5 if postpeak else -0.0
    target *= pulse[i_max]

    autocorr = numpy.correlate(pulse, pulse, mode='full')[-len(pulse):]
    print(autocorr)
    print(sum(pulse * pulse))
    crosscorr = numpy.correlate(target, pulse, mode='full')[-len(pulse):]

    autocorr = autocorr[:length]
    crosscorr = crosscorr[:length]

    coeffs = linalg.inv(linalg.toeplitz(autocorr)).dot(crosscorr)
    coeffs *= pulse[i_max] / max(run_fir_filter(coeffs, pulse))
    return coeffs


def run_fir_filter(coeffs, samples):
    """Process an array of samples with an FIR filter."""
    samples = numpy.ravel(samples)
    return numpy.convolve(coeffs, samples, mode='full')[:len(samples)]


def __make_pulse(**kwargs):
    # pylint: disable=too-many-locals
    """Creates a bipolar pulse for testing purposes.

    This is based on the analytic model presented [here][1]. It
    essentially approximates the cell output pulse as triangular, the
    preamplifier with a first-order low-pass filter, and the analog
    shaper as an ideal CR-(RC)^2 shaping element. An additional
    low-pass filter is added to simulate a small parasitic capacitance
    at the level of the preamplifier.

    All given times are in nanoseconds, all capacitances in picofarads,
    all impedances in ohms, all gains without unit.

    Args:
        sampling_period: Time between two samples.
        t_shift: Time before the first non-zero sample.
        t_end: The duration of the entire pulse.
        t_drift: The calo cell drift time, i.e. the length of the
            triangular pulse.
        tau: The shaping time of the analog shaper, i.e. 1/(R*C).
        c_det: The capacitance of the calo cell.
        c_paras: The parasitic capacitance at the preamp.
        z_in: The preamp's input impedance.
        z_trans: The preamp's transimpedance (V_{out}/I_{in}).
        g_sh: The analog shaper's gain.
        g_lsb: The Layer Sum Board's gain.

    Returns:
        A tuple `(times, samples)` of same-length 1D arrays.
            times: The time coordinate of each sample in seconds.
            samples: The amplitude at each given time.

    [1]: http://cds.cern.ch/record/685385
    """
    from areustools.common import transfer
    # Read all parameters.
    try:
        sampling_period = kwargs.pop('sampling_period')*1e-9
        t_shift = kwargs.pop('t_shift')*1e-9
        t_end = kwargs.pop('t_end')*1e-9
        t_drift = kwargs.pop('t_drift')*1e-9
        tau = kwargs.pop('tau')*1e-9
        c_det = kwargs.pop('c_det')*1e-12
        c_paras = kwargs.pop('c_paras')*1e-12
        z_in = kwargs.pop('z_in')
        z_trans = kwargs.pop('z_trans')
        g_sh = kwargs.pop('g_sh')
        g_lsb = kwargs.pop('g_lsb')
        # `1+` to account for `endpoint=True` further down.
        num_samples = 1+int(round(t_end/sampling_period))
    except KeyError as exc:
        (key,) = exc.args
        raise TypeError('Missing argument: '+repr(key))
    if kwargs:
        (key, _) = kwargs.popitem()
        raise TypeError('Unrecognized argument: '+repr(key))
    # Basic consistency checks.
    assert 1.5 * t_drift + t_shift < t_end
    assert num_samples > 10
    # Create a triangular pulse.
    times = numpy.linspace(0.0, t_end, num_samples, endpoint=True)
    pulse = 1.0 - (times - t_shift) / t_drift
    pulse[numpy.where(numpy.logical_or(pulse > 1.0, pulse < 0.0))] = 0.0
    # Run it through a rough estimation of the analog chain.
    preamp = transfer.AnalogShaper(num_low_passes=1)
    pulse = preamp.shape_time(pulse, sampling_period, c_det * z_in)
    pulse = preamp.shape_time(pulse, sampling_period, c_paras * z_trans)
    pulse = transfer.time_cr_rc2(pulse, sampling_period, tau)
    pulse *= g_sh * g_lsb
    return times, pulse

####################################################################################################

####################################################################################################


def simple_read_in(path_, col_, length_):
    
    if _print_ == 1: print("Reading in "+ path_)
    datafile = open(path_,'r')

    data_ = numpy.empty(length_)

    i = 0
    for line in datafile:
        if i == length_: break
        data_[i] = line.split(",")[col_]
        i += 1

    datafile.close()

    return data_

def simple_plot(array_, color_, save_fig_):
    
    return

def scale_function(func_):

    floor = numpy.mean(func_[:80])
    func_2 = [(val - floor) for val in func_]
    func_2 /= max(func_2)
    return func_2


def getmax_function(func_):

    floor = numpy.mean(func_[:80])
    func_2 = [(val - floor) for val in func_]
    return max(func_2)

def plot_fit_hist(array_, show_, saveto_ = None):

    #Sandard Deviation Histogram with Fit
    mu, std = norm.fit(array_)

    mean_flr = math.floor(numpy.mean(array_))

    xmin = math.floor(mean_flr) - 20
    xmax = math.floor(mean_flr) + 20

    #y_, x_, _ = plt.hist(array_)
    
    y_, x_, _ = plt.hist(array_, bins=19, normed=False, color='b') #range=[xmin,xmax], alpha=0.6, color='b')
    
    #xmin, xmax, ymin, ymax = plt.axis()
    
    #xmin, xmax = plt.xlim()
    x = numpy.linspace(xmin, xmax, 100)
    p = norm.pdf(x, mu, std)*y_.sum()*1.2
    plt.plot(x, p, 'k', linewidth=2)
    
    #title = "Fit results: sample %d,  sigma = %.2f, mu = %.2f" % (index, mu, std)
    title = "Fit results: sigma = %.2f, mu = %.2f" % (std, mu)
    plt.title(title)
    if saveto_ is not None: 
        plt.savefig(saveto_)
    if show_ == 1 or show_ == True: plt.show()
    plt.clf()
    
    return

def get_acf_coeffs(ped):
    
    if _print_ == 1: print("Getting Autocorrelation Coeffs... ")
    
    n = len(ped)    
    mean = numpy.mean(ped)
    c0 = numpy.sum((ped - mean) ** 2) 
    lag = numpy.arange(n) #including lag 0 calculation

    return [((ped[:n-k] - mean) * (ped[k:] - mean)).sum() / float(c0) for k in lag]

def get_overlayed_pulse(path_, colpulse_, isInverted):
    #This overlays pulses/aligns pulses WITHIN each file NOT between files

    data_corrected = numpy.empty((n_blocks, n_points_per_block))
    
    if _print_ == 1: print("Reading in "+ path_)
    datafile = open(path_,'r')

    if _print_ == 1: print("Getting Overlayed Pulse... ")
    
    i = 0
    j = 0
    for line in datafile:

        #Exclude initial 240 points and first "different" pulse
        j += 1
        if j <= exclude_n_points: continue
        
        if i == n_points: break

        #pulse and sample space
        i_p = math.floor(i/n_points_per_pulse)
        i_s1 = i%n_points_per_pulse

        #block and overlayed pulse space
        i_b = math.floor(i_p/n_pulses_per_block)
        corr = i_b*n_pulses_per_block
        i_s2 = i_s1*n_pulses_per_block + i_p%n_pulses_per_block + corr
        i_s2 = i_s2 % n_points_per_block

        #if inverted is not None:
        if isInverted: 
            data_corrected[i_b][i_s2] = -float(line.split(",")[colpulse_])
        else: 
            data_corrected[i_b][i_s2] = float(line.split(",")[colpulse_])
        
        i += 1

    datafile.close()
    
    pulses_blocks_ = numpy.transpose(data_corrected)
    if isInverted: 
        pulse_overlayed_ = [-numpy.mean(pulses_blocks_[i]) for i in xrange(n_points_per_block)]
    else:
        pulse_overlayed_ = [numpy.mean(pulses_blocks_[i]) for i in xrange(n_points_per_block)] 

    return pulse_overlayed_, pulses_blocks_

def align_pulses(baseline_pulse_, shifted_pulse_, d_ij_range, pulses_blocks_ = None):
    #This aligns pulses BETWEEN files but NOT within files
    
    ref_idx1 = numpy.argmin(-baseline_pulse_)
    ref_idx2 = numpy.argmin(baseline_pulse_)


    #ref_idx1 = numpy.argmin(-numpy.gradient(baseline_pulse_))
    #ref_idx2 = numpy.argmin(numpy.gradient(baseline_pulse_))
    

    if _print_ == 1: print('slope max and min indices = ', ref_idx1, ref_idx2)
    
    sum_sqr = 1e12
    #d_ij_range = 75
    sum_sqr_temp = 0
    offset_ = 0
    extra_points = 10
    
    for d_ij in xrange(-d_ij_range, d_ij_range+1):
            
        step_i = (ref_idx1 - extra_points) % len(shifted_pulse_)
        step_f = (ref_idx1 + extra_points) % len(shifted_pulse_)
        sum_sqr_temp = ((baseline_pulse_[step_i : step_f] - shifted_pulse_[step_i+d_ij : step_f+d_ij])**2).sum()

        step_i = (ref_idx2 - extra_points) % len(shifted_pulse_)
        step_f = (ref_idx2 + extra_points) % len(shifted_pulse_)
        sum_sqr_temp += ((baseline_pulse_[step_i : step_f] - shifted_pulse_[step_i+d_ij : step_f+d_ij])**2).sum()
            
        if sum_sqr_temp < sum_sqr:
            sum_sqr = sum_sqr_temp
            offset_ = d_ij
            
        #print(d_ij, sum_sqr_temp)

    #Apply corrections
    pulse_overlayed_corrected_ = numpy.roll(shifted_pulse_, -offset_)
    g_pulse_overlayed_corrected_ = numpy.gradient(pulse_overlayed_corrected_) #numpy.roll(shifted_pulse_, -offset_)
    
    if pulses_blocks_ == None:
        return pulse_overlayed_corrected_, offset_

    else:

        pulses_blocks_corrected_ = numpy.empty((len(pulses_blocks_), len(pulses_blocks_[0])))
        for i_s1 in xrange(n_points_per_block):
            i_s2 = i_s1 - offset_
            i_s2 = i_s2 % n_points_per_block
                     
            pulses_blocks_corrected_[i_s2] = pulses_blocks_[i_s1]

        #These two other methods don't work (?)
        #pulses_blocks_corrected_ = [ pulses_blocks_[(i_s1 - offset_) % n_points_per_block ] for i_s1 in xrange(n_points_per_block)]
        #pulses_blocks_corrected_ = numpy.roll(pulses_blocks_, -offset_)

        return pulse_overlayed_corrected_, pulses_blocks_corrected_, offset_


#### GLOBAL VARIABLES #############################################
_print_ = 0
_show_all_figs_ = 0

exclude_n_points = 800 #initial points to exclude

n_blocks = 10 #18
n_pulses_per_block = 15 #keep
n_points_per_pulse = 400 #240
n_points_per_block = n_points_per_pulse*n_pulses_per_block #6000

n_pulses = n_blocks*n_pulses_per_block
n_points = n_pulses*n_points_per_pulse #points read in

colpulse = 2
colnoise = 1

path_figures = 'figures_temp/'
path_output = 'output/'

isInverted = False
###################################################################

def _testmod():
    """Test this module."""
    import doctest
    doctest.testmod()

    n_inputs = len(sys.argv) - 1

    if _print_ == 1: print("Number of input files = "+ n_inputs)    

    #data = simple_read_in(sys.argv[1], colpulse)
    #x = numpy.arange(len(data))
    #plt.plot(x, data, 'b.')
    #plt.show()
        
    pulse_overlayed = numpy.empty((n_inputs, n_points_per_block))
    g_pulse_overlayed = numpy.empty((n_inputs, n_points_per_block))
    pulses_blocks = numpy.empty((n_points_per_block, n_blocks*n_inputs))

    offset_hist = numpy.empty(n_inputs)
    
    for i_f in xrange(n_inputs):

        #Overlay for 1 file        
        pulse_overlayed[i_f], g_pulse_overlayed[i_f], pulses_blocks_temp = get_overlayed_pulse(sys.argv[i_f+1], colpulse)

        #Adjust timing of pulse to line up with first pulse
        pulse_overlayed[i_f], g_pulse_overlayed[i_f], pulses_blocks_temp, offset = align_pulses(pulse_overlayed[0], pulse_overlayed[i_f], pulses_blocks_temp)

        #Append blocks of pulses in each file to one including all blocks
        for i_s in xrange(n_points_per_block):
            for i_b1 in xrange(n_blocks):
                i_b2 = n_blocks*i_f + i_b1
                pulses_blocks[i_s][i_b2] = pulses_blocks_temp[i_s][i_b1]

        print('File ', i_f, ': Offset --> ', offset)
        
        offset_hist[i_f] = offset

        #Plot all waveforms on eachother
        x = numpy.arange(len(pulse_overlayed[i_f]))
        plt.plot(x[250:1250], pulse_overlayed[i_f][250:1250], 'b.')
        plt.plot(x[250:1250], g_pulse_overlayed[i_f][250:1250], 'g.')
    
    plt.savefig('figures/plot_overlay_100.png')
    if _show_all_figs_ == 1: plt.show()
    plt.clf()
    
    #plt.plot(x[250:1250], pulse_overlayed[0][250:1250], 'b.-')
    #plt.plot(x[250:1250], pulse_overlayed[1][250:1250], 'g.-')
    #plt.plot(x[250:1250], pulse_overlayed[2][250:1250], 'r.-')
    #plt.plot(x[250:1250], pulse_overlayed[3][250:1250], 'k.-')
    #plt.show()
    #plt.clf()
    
    plt.hist(offset_hist, bins=15, normed=True, alpha=0.6, color='b')
    plt.savefig('figures/timing_diffs_100.png')
    if _show_all_figs_ == 1: plt.show()
    plt.clf()
    
    #Get all pulses overlayed ("Master")
    
    pulse_overlayed_master = numpy.empty(len(pulses_blocks))
    pulse_stdev_master = numpy.empty(len(pulses_blocks))
    for i_s2 in xrange(n_points_per_block):
        pulse_overlayed_master[i_s2] = numpy.mean(pulses_blocks[i_s2])                      
        pulse_stdev_master[i_s2] = numpy.std(pulses_blocks[i_s2])

    #gradient
    g_pulse_overlayed_master = numpy.gradient(pulse_overlayed_master, 1)


    print(numpy.argmin(-pulse_overlayed_master))

    #plot master
    x = numpy.arange(len(pulse_overlayed_master))
    plt.plot(x, pulse_overlayed_master, 'b.-'), plt.plot(x, g_pulse_overlayed_master, 'g.-')
    plt.plot(x, pulse_stdev_master, 'r.-')
    plt.savefig('figures/overlayed_master.png')
    if _show_all_figs_ == 1: plt.show()
    plt.show()
    plt.clf()

    plot_fit_hist(pulses_blocks[40])
    
    if _print_ == 1: print("File not written... ")

    #Do OFC Stuff Here
    ped = simple_read_in(sys.argv[1], colnoise, n_points_per_pulse)
    testsample = simple_read_in(sys.argv[12], colpulse, n_points_per_pulse)
    acf_coeffs = get_acf_coeffs(ped)

    E = numpy.empty(n_pulses_per_block)
    ExT = numpy.empty(n_pulses_per_block)

    a = numpy.empty((n_pulses_per_block, n_points_per_pulse))
    b = numpy.empty((n_pulses_per_block, n_points_per_pulse))
    
    acf_coeffs = numpy.zeros(n_points_per_pulse)
    acf_coeffs[0] = 1.
    for i in xrange(n_pulses_per_block):
        a[i], b[i] = calc_of_coeffs(acf_coeffs, pulse_overlayed_master[i::n_pulses_per_block], g_pulse_overlayed_master[i::n_pulses_per_block])
        E[i] = numpy.dot(a[i], pulse_overlayed_master[i::n_pulses_per_block])
        ExT[i] = numpy.dot(b[i], pulse_overlayed_master[i::n_pulses_per_block])
        print(E[i], ExT[i])

    print(max(pulse_overlayed_master))

    #if _show_all_figs_ == 1:
    #plt.show()
    plt.clf()    
    
if __name__ == "__main__":
    _testmod()
