from __future__ import absolute_import
from __future__ import print_function
from __future__ import division

import numpy as np
import csv
import matplotlib.pyplot as plt
from matplotlib.font_manager import FontProperties
from scipy import linalg
import h5py
import sys, getopt

from helperfunctions import *


##########GLOBAL##########
run_to_use = "Feb6"
amp_to_use = "0p8"
gain_to_use = "1x"
make_plots = False
##########################



def calc_of_coeffs(ac, g, amp, dg=None):
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
            `np.gradient`.

    Returns:
        A tuple `(a_coeffs, b_coeffs)` of 1D arrays with the same
        length as `pulse`.

    [1]: Cleland, Stern: Signal processing considerations for liquid
         ionization calorimeters in high rate environment

    """
    ac = np.ravel(ac) / ac[0]
    g = np.ravel(g)
    dg = np.ravel(dg) if dg is not None else np.gradient(g)
    #scale = max(g)
    scale = amp
    g /= scale
    dg /= scale
    # Calculate V = R^{-1}.
    inv_ac = linalg.inv(linalg.toeplitz(ac))
    # Calculate V*g and V*dg only once.
    vg = np.dot(inv_ac, g)
    vdg = np.dot(inv_ac, dg)
    # Calculate helper variables.
    q1 = np.dot(g, vg)
    q2 = np.dot(dg, vdg)
    q3 = np.dot(dg, vg)
    delta = q1*q2 - q3*q3
    # Calculate Lagrange multipliers
    lm_lambda = q2/delta
    lm_kappa = -q3/delta
    lm_mu = q3/delta
    lm_rho = -q1/delta
    # Calculate filter coefficients.
    a_coeffs = lm_lambda*vg + lm_kappa*vdg
    b_coeffs = lm_mu*vg + lm_rho*vdg
    # Reverse order to get canonical coefficient order.
    #return a_coeffs[::-1], b_coeffs[::-1]
    return a_coeffs, b_coeffs



def autocorr(x):
  result = np.correlate(x, x, mode='full')
  return result[result.size // 2:]


if __name__ == "__main__":


    try:
      opts, args = getopt.getopt(sys.argv[1:], "hr:a:g:d:c:p", ["run=", "amp=", "amplitude=", "gain=", "adc=", "channel=", "plots"])
    except getopt.GetoptError:
      print ("Usage: calcOFCs.py -r <run> -a <amplitude [e.g. 0p8]> -g <gain [e.g. 1x]> -d <ADC/COLUTA Number [e.g. 1,2]> -c <Channel Number [e.g 1, 2]>")
      sys.exit(2)
    for opt, arg in opts:
      if (opt == "-h"):
        print ("calcOFCs.py")
        print ("")
        print ("Usage: calcOFCs.py -v <run> -a <amplitude [e.g. 0p8]> -p(optional)")
        print ("")
        print ("Arguments:")
        print ("-r, --run \t\t\t Specify run")
        print ("-a, --amp, --amplitude \t\t Specify amplitude")
        print ("-p, --plots \t\t\t Save Plots")
        sys.exit()
      elif opt in ("-r", "--run"):
        run_to_use = arg
      elif opt in ("-a", "--amp", "--amplitude"):
        amp_to_use = arg
      elif opt in ("-g", "--gain"):
        gain_to_use = arg
      elif opt in ("-p", "--plots"):
        make_plots = True
    print ("Run to use: "+run_to_use)
    print ("Amplitude to use: "+amp_to_use)
    print ("Make Plots: "+str(make_plots))
    print ("These settings can be changed from the command line. Run with the option -h for help")


    plotsDir = checkDir(os.path.abspath(os.path.join(sys.path[0], os.pardir))+"/Plotting/PlotsSaved/Pulses/OFCs_Calculation/"+run_to_use+"/"+amp_to_use+"/")

    f = h5py.File(getRootDir()+"/Physics/Blocks_Final/"+run_to_use+"/Blocks_Amp"+amp_to_use+"_"+gain_to_use+".hdf5", "r")
    block_length = f.attrs["block_length"]
    awg_freq = f.attrs["awg_freq"]
    adc_freq = f.attrs["adc_freq"]
    n_adcs = f.attrs["n_adcs"]
    n_channels = f.attrs["n_channels"]
    adc_freq = f.attrs["adc_freq"]
    awg_freq = f.attrs["awg_freq"]
    n_offset = 1
    outFile = makeH5File(checkDir(sys.path[0]+"/OFCs/"+run_to_use+"/"+amp_to_use+"/"+gain_to_use+"/")+"COLUTA_OFCs.hdf5")

    min_samples = 2
    max_samples = 9
    outFile.attrs.create("min_samples", min_samples)
    outFile.attrs.create("max_samples", max_samples)
    outFile.attrs.create("n_adcs", n_adcs)
    outFile.attrs.create("n_channels", n_channels)
    outFile.attrs.create("adc_freq", adc_freq)
    outFile.attrs.create("awg_freq", awg_freq)

    #Autocorrelation File
    acFile = h5py.File(getRootDir()+"/Physics/OFCs/Autocorrelation/Run_1204/Autocorrelation.hdf5", "r")

    first = True
    for adc in range(n_adcs):
      for channel in range(n_channels):
  
        blocks = f["coluta"+str(adc+1)+"/channel"+str(channel+1)+"/samples"][()] #ADC 1, Channel 1

        g = np.empty(block_length)
        for i in range(block_length):
          g[i] = np.mean(np.transpose(blocks)[i])
          
        dg = np.gradient(g)

        amp = np.max(g)
        fontP = FontProperties()
        fontP.set_size('x-small')
        fig=plt.figure()

        timeStep = (1/awg_freq)*1000 #in ns
        n_phases = int((awg_freq/n_offset)/adc_freq)


        ac = acFile["coluta"+str(adc+1)+"/channel"+str(channel+1)][()]
  


        start_sample = 120
        #start_sample = 142
        print ("Saving OFCs for "+run_to_use.replace("_", " ")+" COLUTA "+str(adc+1)+" Channel "+str(channel+1))
        for nSamples in range(min_samples, max_samples + 1):
            g_25 = [0]*n_phases   #Finely sampled pulse(g) and its derivative(dg) sampled at 40MHz
            dg_25 = [0]*n_phases  #n_phases Copies of each array, one for each of the n_phases phases
            acCoeffs = ac[0:nSamples] #Autocorrelation Function
            #acCoeffs = [0]*nSamples
            #acCoeffs[0] = 1 #Will make the autocorrelation matrix the identity
            for i in range(0, n_phases):
              g_25[i] = g[start_sample + i:start_sample + i + n_phases*nSamples:n_phases]   #Every n_phases samples -> 25ns
              dg_25[i] = dg[start_sample + i:start_sample + i + n_phases*nSamples:n_phases]   #Every n_phases samples -> 25ns


            a = [0]*n_phases  #a[phase][sample]
            b = [0]*n_phases  #b[phase][sample]
            A = [0]*n_phases  #A[phase]
            AT = [0]*n_phases #AT[phase]

            a_master = [0]*8 #a_master[nSamples][phase][sample]

            for i in range(0, n_phases):
              try:
                a[i], b[i] = calc_of_coeffs(acCoeffs, g_25[i], amp, dg_25[i])
              except ValueError:
                a[i] = np.zeros(nSamples)
                b[i] = np.zeros(nSamples)
              for j in range(0, nSamples):
                A[i] += a[i][j]*g_25[i][j]
                AT[i] += b[i][j]*g_25[i][j]
              

            outFile.create_dataset("coluta"+str(adc+1)+"/channel"+str(channel+1)+"/"+str(nSamples)+"S/a", data=a)
            outFile.create_dataset("coluta"+str(adc+1)+"/channel"+str(channel+1)+"/"+str(nSamples)+"S/b", data=b)
            #print ("OFCs saved for "+str(nSamples)+ " samples")
          
            if(make_plots and first):
              #View results of sampling the fine pulse at 40MHz with different combinations of nSamples and phase
              #First set of plots on https://www.nevis.columbia.edu/~akahn/OFCs/Jun6/main.html
              dg_view = [x * 20 for x in dg_25]
              a_view = [x for x in a]
              b_view = [x for x in b]
              for nPhase in range(0, n_phases):
                  ax=fig.add_subplot(111)
                  ax.plot(np.linspace(0, 449*timeStep, 450), g[0:450], 'b.-', label='Finely Sampled Pulse')
                  ax.plot(np.linspace((start_sample +  nPhase)*timeStep, (start_sample + nPhase + n_phases*(nSamples-1))*timeStep, nSamples), g_25[nPhase][0:nSamples], 'r.-', label='40MHz Pulse')
                  ax.plot(np.linspace((start_sample +  nPhase)*timeStep, (start_sample + nPhase + n_phases*(nSamples-1))*timeStep, nSamples), dg_25[nPhase][0:nSamples], 'g.-', label='40MHz Pulse Derivative')
                  ax.plot(np.linspace((start_sample +  nPhase)*timeStep, (start_sample + nPhase + n_phases*(nSamples-1))*timeStep, nSamples), a_view[nPhase][0:nSamples], 'm.-', label='OFCs: a')
                  ax.plot(np.linspace((start_sample +  nPhase)*timeStep, (start_sample + nPhase + n_phases*(nSamples-1))*timeStep, nSamples), b_view[nPhase][0:nSamples], 'y.-', label='OFCs: b')

                  ax.grid()
                  plt.title('Samples: '+str(nSamples)+' Phase: '+str(nPhase))
                  plt.ylabel('Counts')
                  plt.xlabel('Time [ns]')
                  ax.legend(loc='upper right', prop = fontP)
                  ax.text(0.8, 0.78, 'Amp: '+str(np.round(amp, 2)), verticalalignment='top', horizontalalignment='left', transform=ax.transAxes)
                  ax.text(0.8, 0.74, 'Amp OFC: '+str(np.round(A[nPhase], 2)), verticalalignment='top', horizontalalignment='left', transform=ax.transAxes)
                  ax.text(0.8, 0.70, 'AT: '+'{0:1.3g}'.format(AT[nPhase]), verticalalignment='top', horizontalalignment='left', transform=ax.transAxes)
                  plt.savefig(plotsDir+"Pulse_S"+str(nSamples)+"_P"+str(nPhase)+".png")
                  print("Created file: "+plotsDir+"Pulse_S"+str(nSamples)+"_P"+str(nPhase)+".png")
                  plt.clf()
        first=False

    print ("File Created: "+checkDir(sys.path[0]+"/OFCs/"+run_to_use+"/"+amp_to_use+"/"+gain_to_use+"/")+"COLUTA_OFCs.hdf5")





