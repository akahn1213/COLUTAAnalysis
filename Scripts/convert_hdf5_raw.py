import h5py
import numpy as np
import helperfunctions as hf
import os
import glob
import sys
import argparse
import time
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

def isSequence(arg):
  """Determines if arg is a sequence. See https://stackoverflow.com/questions/1835018/"""
  return (not hasattr(arg, "strip") and
          (hasattr(arg, "__getitem__") or
           hasattr(arg, "__iter__")))

def setHDF5Attributes(hdf5File,**kwargs):
  """Create and attach attributes specified in kwargs to the HDF5 group or dataset hdf5Object"""
  # with h5py.File(hdf5Object,'a') as hdf5File:
  for key, value in kwargs.items():
    if isSequence(value) and not isinstance(value, np.generic):
      if type(value[0]) is str:
        hdf5File.attrs.create(key,value,dtype=h5py.special_dtype(vlen=str))
      else:
        hdf5File.attrs.create(key,value)
    else:
      if type(value) is str:
        hdf5File.attrs.create(key,value,dtype=h5py.special_dtype(vlen=str))
      else:
        hdf5File.attrs.create(key,value)

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
      sys.exit(0)

  #Make the output directory
  output_dir = hf.getRootDir()+"/Data/Processed/Pulses/"+run_name
  if not os.path.exists(output_dir):
    os.mkdir(output_dir)
  print("Writing processed output to: " + output_dir)

  time_start = time.time()

  #Define measurement ranges per amplitude
  gain_array = ['AG', '4x', '1x']
  meas_dict = defaultdict(lambda : defaultdict(list))
  for measurement in range(input_file.attrs['n_measurements']):
    gain = gain_array[input_file['Measurement_'+str(measurement)+'/coluta1/channel1'].attrs['gain']]
    pulse_amp = input_file['Measurement_'+str(measurement)].attrs['pulse_amplitude'].replace('.','p')
    meas_dict[pulse_amp][gain].append(measurement)

  default_attrs = ['laurocDynamicRange','n_pulse_timestep','n_pulses','n_samples','n_samples_per_pulse','pulse_amplitude','run_type','timestamp','dc_offset']

  for amp, gain_dict in meas_dict.items():
    for gain, meas_range in gain_dict.items():

      #Output name format: Pulse_Amp#p###_(gain).hdf5
      out_file = h5py.File(output_dir+"/Pulse_Amp"+amp+"_"+gain+".hdf5", "w")

      #Attributes:
      setHDF5Attributes(out_file, n_adcs = input_file.attrs["n_adcs"],
                                  n_measurements = len(meas_range),
                                  n_pulses = 30,
                                  n_samples_per_pulse = 64,
                                  adc_freq = 40,
                                  awg_freq = 1200
                        )
      # out_file.attrs.create("n_adcs", input_file.attrs["n_adcs"])
      # #out_file.attrs.create("n_channels", 2)
      # out_file.attrs.create("n_measurements", len(meas_range))
      # out_file.attrs.create("n_pulses", 30)
      # out_file.attrs.create("n_samples_per_pulse", 64)
      # out_file.attrs.create("adc_freq", 40)
      # out_file.attrs.create("awg_freq", 1200)

      for meas_i in range(len(meas_range)):
        cut_attrs = [attr for attr in input_file["Measurement_"+str(meas_range[meas_i])].attrs if attr not in default_attrs]
        cut_attrs_dict = {attr:input_file["Measurement_"+str(meas_range[meas_i])].attrs[attr] for attr in cut_attrs}
        out_file.create_group("Measurement_"+str(meas_i))
        setHDF5Attributes(out_file["Measurement_"+str(meas_i)], **cut_attrs_dict)
                          # shaper_constants = input_file["Measurement_"+str(meas_range[meas_i])].attrs["shaper_constants"],
                          # hg_lg_c2 = input_file["Measurement_"+str(meas_range[meas_i])].attrs["hg_lg_c2"])
        for adc in range(out_file.attrs["n_adcs"]):
          #for channel in range(out_file.attrs["n_channels"]):
          out_file.create_group("Measurement_"+str(meas_i)+"/coluta"+str(adc+1))
          # out_file["Measurement_"+str(meas_i)+"/coluta"+str(adc+1)].attrs.create("channels", input_file["Measurement_"+str(meas_range[meas_i])+"/coluta"+str(adc+1)].attrs["channels"])
          # out_file["Measurement_"+str(meas_i)].attrs.create("shaper_constants", input_file["Measurement_"+str(meas_range[meas_i])].attrs["shaper_constants"])
          # out_file["Measurement_"+str(meas_i)].attrs.create("hg_lg_c2", input_file["Measurement_"+str(meas_range[meas_i])].attrs["hg_lg_c2"])
          setHDF5Attributes(out_file["Measurement_"+str(meas_i)+"/coluta"+str(adc+1)],
                            channels = input_file["Measurement_"+str(meas_range[meas_i])+"/coluta"+str(adc+1)].attrs["channels"])
          for channel in out_file["Measurement_"+str(meas_i)+"/coluta"+str(adc+1)].attrs["channels"]:
            if channel == 'frame': continue
            address = "Measurement_"+str(meas_range[meas_i])+"/coluta"+str(adc+1)+"/"+channel   #"Measurement_(9->13)/coluta1/channel(1,2)/samples"
            SAR_weights = input_file[address].attrs["SAR_weights"]
            raw_data = input_file[address+"/raw_data"]
            bits = raw_data[:,1]
            samples = np.dot(raw_data,SAR_weights)
            out_file.create_dataset("Measurement_"+str(meas_i)+"/coluta"+str(adc+1)+"/"+channel+"/samples", data=samples)
            out_file.create_dataset("Measurement_"+str(meas_i)+"/coluta"+str(adc+1)+"/"+channel+"/bits", data=bits)
            #The above lines maps measurements in the current meas_range in the input file to measurements 0->len(meas_range) in the output file
      out_file.close()
  time_total = time.time() - time_start
  print('Took ' + str(np.round(time_total,3)) + ' seconds')
