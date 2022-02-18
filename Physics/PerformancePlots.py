from scipy import stats 
import numpy as np
from scipy.stats import norm
import csv
import os
from scipy import linalg
import sys
from optparse import OptionParser
import argparse
import h5py
import atlas_plotter as plotter
from helperfunctions import *
import copy

#########Global##########
nSampleRange = 7
lowest_nSamples = 4
nMaxSamples = nSampleRange+lowest_nSamples


filterDir = checkDir(getRootDir()+"/Physics/Filter_Output/")


if __name__ == "__main__":

  #Get directories to run in
  parser = argparse.ArgumentParser()
  parser.add_argument("-r", "--runs", default = [], type=str, nargs='+',
                     help="list of runs (directories in Physics/Filter_Output/) to include")
  #Allow user input for choice of runs
  args = parser.parse_args()
  #runs = args.runs
  runs = [filterDir + run + "/" for run in args.runs]
  directories = []

  #Search for runs if none are provided
  if runs == []:
    directories = glob.glob(filterDir+"*/")
    if directories == []:
      print("No runs found in "+filterDir+" exiting...")

      sys.exit(0)
    for directory in directories:
      if glob.glob(directory+"*/*.hdf5"):
        runs.append(directory)
    print("Found filter output in the following directories:\n")
    for run in runs:
      print("-%s" % short(run))
    #continue_script = input("\nContinue running? [y/n] (default=y): ")
    continue_script = "y"
    if(continue_script == "n"):
      sys.exit(0)

  filterType = "OF"

  for run in runs: #Create performance plots
    #Load the filter output
    filter_out = h5py.File(filterDir+short(run)+filterType+'/'+filterType+'_Output.hdf5', 'r') 
    #filter_out_p16 = h5py.File(filterDir+short(run)+filterType+'/'+filterType+'_Output_P16.hdf5', 'r') 
    # filter_out_p16 = h5py.File(filterDir+short(run)+filterType+'/'+filterType+'_Output_P16.hdf5', 'r')
    #filter_out_lo = h5py.File(filterDir+short(run)+filterType+'/'+filterType+'_Output_LO.hdf5', 'r') 
    #filter_out_peak = h5py.File(filterDir+short(run)+filterType+'/'+'Peak_Output.hdf5', 'r') 
    filter_out_peak = h5py.File(filterDir+short(run)+filterType+'/'+filterType+'_Output_Peak.hdf5', 'r') 
    amps = filter_out.attrs["amps"]    
    gains = filter_out.attrs["gains"]    
    amps = amps.astype("U13")
    gains = gains.astype("U13")
    #Get the number of ADCs and channels from the OFCs file. MAY BREAK in the future
    n_adcs = filter_out.attrs["n_adcs"]
    n_channels = filter_out.attrs["n_channels"]
    max_samples = filter_out.attrs["max_samples"]
    min_samples = filter_out.attrs["min_samples"]
    adc_freq = filter_out.attrs["adc_freq"]
    awg_freq = filter_out.attrs["awg_freq"]
    awg_amps = getAmpNums(amps)
    energy_mean = np.empty(len(amps))  
    energy_mean_peak = np.empty(len(amps))  
    energy_std_peak = np.empty(len(amps))  
    timing_mean = np.empty(len(amps))
    energy_std = np.empty(len(amps))  
    timing_std = np.empty(len(amps))
    setup_params = {"amps":amps,
                    "awg_amps":awg_amps,
                    "run":short(run),  
                    "adc_freq":adc_freq,
                    "awg_freq":awg_freq,
                    "filter_name": "Optimal Filter",
                    "filter_info": str(min_samples)+" Samples"}
    data = [{"energy_mean":np.empty(len(amps)), 
            "energy_std":np.empty(len(amps)),
            "timing_mean":np.empty(len(amps)), 
            "timing_std":np.empty(len(amps)),
            "gain":""}]

    phase = 15
    for adc in range(n_adcs):
      for channel in range(n_channels):
        if adc==1 and channel==1: continue
        #for n_samples in range(min_samples, max_samples):
        for n_samples in range(5, 6):
          setup_params.update({"filter_info": str(n_samples) + " Samples"})
          multidata=[]
          for gain in gains:
            #makePlots(adc, channel, n_samples, gain)
            #Data made for each amplitude available
            for amp in range(len(amps)):
              energy_mean[amp] = np.mean(filter_out["/coluta"+str(adc+1)+"/channel"+str(channel+1)+"/"+gain+"/Amp"+amps[amp]+"/"+str(n_samples)+"S/P"+str(phase)+"/energy"][()])
              energy_mean_peak[amp] = np.mean(filter_out_peak["/coluta"+str(adc+1)+"/channel"+str(channel+1)+"/"+gain+"/Amp"+amps[amp]+"/"+str(n_samples)+"S/P"+str(phase)+"/energy"][()])
              energy_std_peak[amp] = np.std(filter_out_peak["/coluta"+str(adc+1)+"/channel"+str(channel+1)+"/"+gain+"/Amp"+amps[amp]+"/"+str(n_samples)+"S/P"+str(phase)+"/energy"][()])
              energy_std[amp] = np.std(filter_out["/coluta"+str(adc+1)+"/channel"+str(channel+1)+"/"+gain+"/Amp"+amps[amp]+"/"+str(n_samples)+"S/P"+str(phase)+"/energy"][()])
              timing_mean[amp] = np.mean(filter_out["/coluta"+str(adc+1)+"/channel"+str(channel+1)+"/"+gain+"/Amp"+amps[amp]+"/"+str(n_samples)+"S/P"+str(phase)+"/timing"][()])
              timing_std[amp] = np.std(filter_out["/coluta"+str(adc+1)+"/channel"+str(channel+1)+"/"+gain+"/Amp"+amps[amp]+"/"+str(n_samples)+"S/P"+str(phase)+"/timing"][()])
            #Data using only the highest available amplitude
            timing = filter_out["/coluta"+str(adc+1)+"/channel"+str(channel+1)+"/"+gain+"/Amp"+amps[len(amps)-1]+"/"+str(n_samples)+"S/P"+str(phase)+"/timing"][()]
            timing_p16 = filter_out["/coluta"+str(adc+1)+"/channel"+str(channel+1)+"/"+gain+"/Amp"+amps[len(amps)-1]+"/"+str(n_samples)+"S/P"+str(phase+1)+"/timing"][()]
            #Store data into a dictionary to be plotted
            data[0].update({"energy_mean":energy_mean})
            data[0].update({"energy_mean_peak":energy_mean_peak})
            data[0].update({"energy_std_peak":energy_std_peak})
            data[0].update({"energy_std":energy_std})
            data[0].update({"timing_mean":timing_mean})
            data[0].update({"timing_std":timing_std})
            data[0].update({"timing":timing})
            data[0].update({"timing_p16":timing_p16})
            data[0].update({"gain":gain})
            multidata.append(copy.deepcopy(data[0]))
            plotter.makePlots(data, adc, channel, setup_params, short(run))
          plotter.makePlots(multidata, adc, channel, setup_params, short(run))

