import h5py
import numpy as np
import helperfunctions as hf
import os
import glob
import sys
import argparse
from builtins import input
from collections import defaultdict



def get_gain(bits):
  gain = ""
  if(bits[0] == 1): gain = "1x"
  elif(bits[0] == 0): gain = "4x"
  else: print("Bits found to have values != 0 or 1. Gain type cannot be determined")
  for i in range(2048):
    if(bits[i] == 0 and gain == "1x"): gain = "AG"
    if(bits[i] == 1 and gain == "4x"): gain = "AG"
  return gain

dataDir = hf.checkDir(hf.getRootDir()+"/Data/Raw/Pulses/")

#Get directories to run in
parser = argparse.ArgumentParser()
parser.add_argument("-r", "--runs", default = [], type=str, nargs='+',
                   help="list of runs (directories in Data/Processed/) to include")
#Allow user input for choice of runs
args = parser.parse_args()
debug_plots=False
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
    print("-%s" % hf.short(run))
  continue_script = input("\nContinue running? [y/n] (default=y): ")
  if(continue_script == "n"):
    sys.exit(0)

#Raw Data path
raw_dir = hf.getRootDir()+"/Data/Raw/Pulses/"

for run in runs:
  run_name = hf.short(run)[:-1]
  #Choose Input File
  #----------------------------------
  infile_name = run_name + "_Output.hdf5"
  #----------------------------------

  #Try to find the file in either the raw data path, or in this directory
  try:
    input_file = h5py.File(raw_dir+infile_name, "r")
  except:
    try:
      input_file = h5py.File(infile_name, "r")
    except:
      print("Input file not found in either ../Data/Raw/Pulses/ or ../Scripts/")

  #Make the output directory
  output_dir = hf.getRootDir()+"/Data/Processed/Pulses/"+run_name
  if not os.path.exists(output_dir):
    os.mkdir(output_dir)
  print("Writing processed output to: " + output_dir)

  #Define measurement ranges per amplitude
  meas_dict = defaultdict(list)
  for measurement in range(input_file.attrs['n_measurements']):
    pulse_amp = input_file['Measurement_'+str(measurement)].attrs['pulse_amplitude'].replace('.','p')
    meas_dict[pulse_amp].append(measurement)
  #meas_dict.update({"1p0": np.arange(1, 101)}) #Measurements 1 -> 100 : Amplitude 1.0
  # meas_dict.update({"0p75": np.arange(0, 100)})
  #meas_dict.update({"0p5": np.arange(151, 201)})
  #meas_dict.update({"0p1": np.arange(201, 231)})


  for amp, meas_range in meas_dict.items():
    print(amp)
    print(meas_range)

    #Choose the gain by investigating the bits data
    gain = get_gain(input_file["Measurement_"+str(meas_range[0])+"/coluta1/channel1/bits"][()])


    #Output name format: Pulse_Amp#p###_(gain).hdf5
    out_file = h5py.File(output_dir+"/Pulse_Amp"+amp+"_"+gain+".hdf5", "w")
    
    #Attributes:
    out_file.attrs.create("n_adcs", input_file.attrs["n_adcs"])
    #out_file.attrs.create("n_channels", 2)
    out_file.attrs.create("n_measurements", len(meas_range))
    out_file.attrs.create("n_pulses", 30)
    out_file.attrs.create("n_samples_per_pulse", 64)
    out_file.attrs.create("adc_freq", 40)
    out_file.attrs.create("awg_freq", 1200)
    
    
    for meas_i in range(len(meas_range)):
      for adc in range(out_file.attrs["n_adcs"]):
        #for channel in range(out_file.attrs["n_channels"]):
        out_file.create_group("Measurement_"+str(meas_i)+"/coluta"+str(adc+1))
        out_file["Measurement_"+str(meas_i)+"/coluta"+str(adc+1)].attrs.create("channels", input_file["Measurement_"+str(meas_range[meas_i])+"/coluta"+str(adc+1)].attrs["channels"])
        #for channel in input_file["Measurement_"+str(meas_i)+"/coluta"+str(adc+1)].attrs["channels"]:
        for channel in out_file["Measurement_"+str(meas_i)+"/coluta"+str(adc+1)].attrs["channels"]:
          if channel == 'frame': continue
          address = "Measurement_"+str(meas_range[meas_i])+"/coluta"+str(adc+1)+"/"+channel   #"Measurement_(9->13)/coluta1/channel(1,2)/samples"
          # print(address)
          out_file.create_dataset("Measurement_"+str(meas_i)+"/coluta"+str(adc+1)+"/"+channel+"/samples", data=input_file[address+"/samples"][()])
          out_file.create_dataset("Measurement_"+str(meas_i)+"/coluta"+str(adc+1)+"/"+channel+"/bits", data=input_file[address+"/bits"][()])
          #The above lines maps measurements in the current meas_range in the input file to measurements 0->len(meas_range) in the output file
    
    out_file.close()
