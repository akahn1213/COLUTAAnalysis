from scipy import stats 
import numpy as np
from scipy.stats import norm
import csv
import os
import matplotlib.pyplot as plt
from matplotlib.font_manager import FontProperties
from scipy import linalg
import sys
import glob
from optparse import OptionParser
import argparse
import h5py
from builtins import input

from helperfunctions import *

#########Global##########
nSampleRange = 7
lowest_nSamples = 4
nMaxSamples = nSampleRange+lowest_nSamples

#ofcs_to_use = "Feb6/0p8/"
ofcs_to_use = "Run_1232/0p75/1x/"

dataDir = checkDir(getRootDir()+"/Physics/Blocks_Final/")
ofcsDir = checkDir(getRootDir()+"/Physics/OFCs/"+ofcs_to_use)
outDir = checkDir(getRootDir()+"/Physics/Filter_Output/")

filterDict = {"OF": "Optimal Filter", "WF": "Wiener Filter"}

#phase = 8
#phase = 16

#gains = ["1x", "AG"]
# gains = ["1x"]
gains = ['1x','4x','AG']

##########################



def apply_filter(blocks, adc, channel, amp, gain):
  global filter_data
  global filter_out 
  global filter_out_peak 
  global filterType
  global filterDict
  global min_samples
  global max_samples
  global awg_freq

  #start_sample = 142
  #n_phases = 15
  start_sample = 120
  n_phases = 30

  timeUnit = (1./awg_freq)*1000. #Sampling rate of finely sampled pulse in ns

  if filterType == "OF": #Optimal Filter
    #Apply OFCs, store energy and timing


    for phase in range(n_phases):    
      #Iterate through #Samples in ofcs file
      for n_samples in range(min_samples, max_samples+1): 
        #Collect OFCs 
        ofcs_a = filter_data["/coluta"+str(adc+1)+"/channel"+str(channel+1)+"/"+str(n_samples)+"S/a"][(phase)] 
        ofcs_b = filter_data["/coluta"+str(adc+1)+"/channel"+str(channel+1)+"/"+str(n_samples)+"S/b"][(phase)] 
        #Apply OFCs to the blocks
        energy = np.empty(len(blocks))
        timing = np.empty(len(blocks))
        energy_peak = np.empty(len(blocks))
        timing_peak = np.empty(len(blocks))
        for b in range(len(blocks)):
          tmp_energy = 0
          tmp_time = 0
          for s in range(n_samples):
            tmp_energy += blocks[b][start_sample+phase+n_phases*s]*ofcs_a[s]
            tmp_time += blocks[b][start_sample+phase+n_phases*s]*ofcs_b[s]

          tmp_energy_peak = blocks[b][165]

          tmp_time = tmp_time/tmp_energy
          tmp_time_peak = tmp_time/tmp_energy_peak
          energy[b] = tmp_energy
          timing[b] = tmp_time*timeUnit
          energy_peak[b] = tmp_energy_peak
          timing_peak[b] = tmp_time_peak*timeUnit
        filter_out.create_dataset("/coluta"+str(adc+1)+"/channel"+str(channel+1)+"/"+gain+"/Amp"+amp+"/"+str(n_samples)+"S/P"+str(phase)+"/energy", data = energy)
        filter_out.create_dataset("/coluta"+str(adc+1)+"/channel"+str(channel+1)+"/"+gain+"/Amp"+amp+"/"+str(n_samples)+"S/P"+str(phase)+"/timing", data = timing)
        filter_out_peak.create_dataset("/coluta"+str(adc+1)+"/channel"+str(channel+1)+"/"+gain+"/Amp"+amp+"/"+str(n_samples)+"S/P"+str(phase)+"/energy", data = energy_peak)
        filter_out_peak.create_dataset("/coluta"+str(adc+1)+"/channel"+str(channel+1)+"/"+gain+"/Amp"+amp+"/"+str(n_samples)+"S/P"+str(phase)+"/timing", data = timing_peak)


if __name__ == "__main__":

  #Get directories to run in
  parser = argparse.ArgumentParser()
  parser.add_argument("-r", "--runs", default = [], type=str, nargs='+',
                     help="list of runs (directories in Physics/Blocks_Final/) to include")
  parser.add_argument("-o", "--ofcrun", default = "", type=str, nargs=1,
                     help="Run to use OFCs from")
  parser.add_argument("-a", "--amp", default = "0p75", type=str, nargs=1,
                     help="Amplitude to use")
  parser.add_argument("-g", "--gain", default = "1x", type=str, nargs=1,
                     help="Gain to use")
  #Allow user input for choice of runs
  args = parser.parse_args()
  #runs = args.runs
  runs = [dataDir + run + "/" for run in args.runs]
  directories = []

  #Search for runs if none are provided
  if runs == []:
    directories = glob.glob(dataDir+"*/")
    if directories == []:
      print("No runs found in "+dataDir+" exiting...")
      sys.exit(0)
    for directory in directories:
      if glob.glob(directory+"*.hdf5"):
        runs.append(directory)
    print("Found blocks in the following directories:\n")
    for run in runs:
      print("-%s" % short(run))
    continue_script = input("\nContinue running? [y/n] (default=y): ")
    if(continue_script == "n"):
      sys.exit(0)

  use_analysis_run_for_ofcs = False
  ofc_run = args.ofcrun
  if(ofc_run == ""):
    use_analysis_run_for_ofcs = True


  filterType = "OF"

  
  
  for run in runs: #Create performance plots
    #Load OFCs
    if(use_analysis_run_for_ofcs): ofc_run = run
    ofcsDir = checkDir(getRootDir()+"/Physics/OFCs/"+short(ofc_run)+"/"+str(args.amp[0])+"/"+str(args.gain[0])+"/")
    filter_data = h5py.File(ofcsDir+'COLUTA_OFCs.hdf5', 'r') 

    print("Applying: " + filterDict[filterType] + " using blocks for run " + short(run))

    amps = getAmps(run)    
    n_adcs = filter_data.attrs["n_adcs"]
    n_channels = filter_data.attrs["n_channels"]
    min_samples = filter_data.attrs["min_samples"]
    max_samples = filter_data.attrs["max_samples"]
    adc_freq = filter_data.attrs["adc_freq"]
    awg_freq = filter_data.attrs["awg_freq"]
    filter_out = makeH5File(checkDir(outDir+short(run)+filterType+"/")+filterType+'_Output.hdf5', 'core')
    filter_out.attrs.create("min_samples", min_samples)
    filter_out.attrs.create("max_samples", max_samples)
    filter_out.attrs.create("n_adcs", n_adcs)
    filter_out.attrs.create("n_channels", n_channels)
    filter_out.attrs.create("adc_freq", adc_freq)
    filter_out.attrs.create("awg_freq", awg_freq)
    filter_out.attrs.create("amps", amps, dtype = h5py.special_dtype(vlen=str))
    filter_out.attrs.create("gains", gains, dtype = h5py.special_dtype(vlen=str))


    filter_out_peak = makeH5File(checkDir(outDir+short(run)+filterType+"/")+filterType+'_Output_Peak.hdf5', 'core')
    filter_out_peak.attrs.create("min_samples", min_samples)
    filter_out_peak.attrs.create("max_samples", max_samples)
    filter_out_peak.attrs.create("n_adcs", n_adcs)
    filter_out_peak.attrs.create("n_channels", n_channels)
    filter_out_peak.attrs.create("adc_freq", adc_freq)
    filter_out_peak.attrs.create("awg_freq", awg_freq)
    filter_out_peak.attrs.create("amps", amps, dtype = h5py.special_dtype(vlen=str))
    filter_out_peak.attrs.create("gains", gains, dtype = h5py.special_dtype(vlen=str))



    print("Sample Range: " + str(min_samples) + " to " + str(max_samples))
    print("Amplitudes Found: "+str(amps).replace("p","."))
    for amp in amps:
      for gain in gains:
        #Initialize output file
        blocks_file = h5py.File(run+'Blocks_Amp'+amp+'_'+gain+'.hdf5', 'r') #Load a file given the amplitude, gain, and run
        #Iterate through the number of adcs
        for adc in range(n_adcs):
          for channel in range(n_channels):
            #Apply the filter to the current blocks
            blocks = blocks_file["/coluta"+str(adc+1)+"/channel"+str(channel+1)+"/samples"][()]
            apply_filter(blocks, adc, channel, amp, gain) #Apply the filter to the blocks and save the energy and timing information in a separate file
      print(filterDict[filterType] + " applied for Amplitude " + amp.replace("p", "."))
    filter_out.close()
    filter_out_peak.close()
    print("Filter Application Completed")
