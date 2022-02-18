import numpy as np
import csv
import matplotlib.pyplot as plt
from matplotlib.font_manager import FontProperties
from scipy import linalg
import os
import glob
import sys
import h5py
from scipy import stats 
from scipy.stats import norm
from optparse import OptionParser
import argparse
import atlas_plotter as plotter
from helperfunctions import *
import time

dataDir = checkDir(getRootDir()+"/Data/Processed/Pedestal/")


def save_blocks(f, s, b):
  for adc in range(len(s)):
    f.create_group("adc"+str(adc))
    f["/adc"+str(adc)].create_dataset("samples", data = s[adc])
    f["/adc"+str(adc)].create_dataset("bits", data = b[adc])
  f.close()

def get_data_from_hdf5(f, n_measurements):
  global n_adcs
  global n_channels
  global n_samples_per_run
  data = np.zeros(shape=(n_adcs, n_channels, n_measurements, n_samples_per_run))
  for adc in range(n_adcs):
    for channel in range(n_channels):
      if(adc == 1 and channel == 1): continue
      #if(f["Measurement_0/coluta"+str(adc+1)].attrs["channels"][channel] == 'frame'): continue
      for measurement in range(n_measurements):
        data[adc][channel][measurement] = f["Measurement_"+str(measurement)+"/coluta"+str(adc+1)+"/channel"+str(channel+1)+"/samples"][()]
  return data

def makeAttributes(f):
  global block_length
  global adc_freq
  global awg_freq
  global n_adcs
  global n_channels
  f.attrs.create("block_length", block_length)
  f.attrs.create("adc_freq", adc_freq)
  f.attrs.create("awg_freq", awg_freq)
  f.attrs.create("n_adcs", n_adcs)
  f.attrs.create("n_channels", n_channels)









if __name__ == "__main__":

  t_start = time.time()

  print("\n\nCOLUTA Pulse Analysis\n\n")

  #Get directories to run in
  parser = argparse.ArgumentParser()
  parser.add_argument("-r", "--runs", default = [], type=str, nargs='+',
                     help="list of runs (directories in Data/Processed/) to include")
  #Allow user input for choice of runs
  args = parser.parse_args()
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
    print("Found data in the following directories:\n")
    for run in runs:
      print("-%s" % short(run))
    continue_script = input("\nContinue running? [y/n] (default=y): ")
    if(continue_script == "n"):
      sys.exit(0)
     

  #Run over data found in each run's directory
  for run in runs:
    #n_files = getNFiles(run+"Pulse_Amp"+amps[0]+"_1x") #TODO Less hard coding
    tmpFile = h5py.File(glob.glob(run+"Pedestal*.hdf5")[0], 'r') #TODO Less hard coding
    acFile = makeH5File(checkDir(sys.path[0]+"/OFCs/Autocorrelation/"+short(run))+"Autocorrelation.hdf5")

    #Collect data from attributes
    #n_adcs = tmpFile.attrs["n_adcs"]
    n_adcs = 2
    #n_channels = tmpFile.attrs["channels"]
    n_channels = 2
    n_samples_per_run = tmpFile.attrs["n_samples"]
    n_measurements = tmpFile.attrs["n_measurements"]
    
    print("=================================================================")
    print("= Processing Pedestal Run")
    print("= Setup: "+str(n_adcs)+" ADCs per board, "+str(n_channels)+" channels per ADC")
    for adc in range(n_adcs):
      for channel in range(n_channels):
        print("ADC"+str(adc)+"CHANNEL"+str(channel))
        data = get_data_from_hdf5(tmpFile, n_measurements)
        ac = []
        for meas in range(n_measurements):
          tmp_data = data[adc][channel][meas]
          tmp_mean = np.mean(tmp_data)
          tmp_ac = np.zeros(int(len(tmp_data)/2))
          for i in range(len(ac)):
            n = 0
            d = 0
            for j in range(len(tmp_data)):
              x = tmp_data[j] - tmp_mean
              n += x * (tmp_data[(i+j)%len(tmp_data)] - tmp_mean)
              d += x*x
            tmp_ac[i] = n/d
          ac.append(tmp_ac)
        ac = np.mean(ac, axis=0)
        acFile.create_dataset("coluta"+str(adc+1)+"/channel"+str(channel+1), data=ac) 
        #plt.plot(ac[0:150], '.-')
        #plt.title("Autocorrelation function: COLUTA"+str(adc+1) + " Channel "+str(channel+1))
        #plt.show()


        print("Running analysis for ADC "+str(adc)+", Channel "+str(channel)+" ...")
        #The input file has a 1D array for each channel, the output file with have a 2D array of blocks for each channel, but the file structure will remain the same
        plotter.make_pedestal_plots(data[adc][channel], adc, channel, run) #Run the analysis, and store the pulse blocks
        

  print("This process took: " + str(np.round((time.time() - t_start), 3))+" seconds to complete")


