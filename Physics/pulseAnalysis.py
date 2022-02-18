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
from helperfunctions import *
import time
from builtins import input

dataDir = checkDir(getRootDir()+"/Data/Processed/Pulses/")

def interlace(samples, bits):
  global block_length
  global pulse_length
  new_samples = np.empty(block_length)
  new_bits = np.empty(block_length)
  for s in range(block_length-1):
    new_samples[s] = samples[(s*pulse_length)%(block_length-1)]
    new_bits[s] = bits[(s*pulse_length)%(block_length-1)]
  new_samples[block_length-1] = samples[block_length-1]
  new_bits[block_length-1] = bits[block_length-1]
  return new_samples, new_bits

def get_blocks(samples, bits):
  #Create a block for each full cycle of pulses
  global block_length
  global n_files
  global n_blocks
  global n_pulses
  global pulse_length
  n_blocks_total = n_blocks * len(samples)
  pulse_blocks = np.empty(shape=(n_blocks_total, block_length))
  bit_blocks = np.empty(shape=(n_blocks_total, block_length))
  for file_num in range(len(samples)):
    #Level the data (only modifying 4x samples)
    debug_means = []
    mean = np.mean(samples[file_num][0:5])
    debug_means.append(mean)
    for i in range(1, n_pulses):
      current_pulse_mean = np.mean(samples[file_num][pulse_length*i:pulse_length*i + 5])
      debug_means.append(current_pulse_mean)
      for s in range(pulse_length):
        if int(bits[file_num][pulse_length*i + s]) == 0:
          samples[file_num][pulse_length*i + s] += mean - current_pulse_mean
    if(debug_plots):
      plt.plot(debug_means, '.-')
      plt.title("get_blocks: Pulse Means, Std: "+str(np.std(debug_means)))
      plt.show()
    for i in range(n_blocks):
      #Interlace the pulses
      pulse_blocks[file_num*n_blocks + i], bit_blocks[file_num*n_blocks + i] = interlace(samples[file_num][block_length*i:block_length*(i+1)], bits[file_num][block_length*i:block_length*(i+1)])
    if(debug_plots):
      #plt.plot(np.roll(samples[file_num], 1000), '.-')
      plt.plot(samples[file_num], '.-')
      plt.title("get_blocks "+str(file_num)+": Original Samples, std[0:5] = "+str(np.std(samples[file_num][0:5])))
      plt.show()
      plt.plot(pulse_blocks[file_num], '.-')
      plt.title("get_blocks "+str(file_num)+": Interlaced Pulse, std[0:5] = "+str(np.std(pulse_blocks[file_num][0:5])))
      plt.show()

  return np.asarray(pulse_blocks), np.round(np.asarray(bit_blocks)).astype(int)

def get_gain_constant(pb1, pb4, s1, s4, pbA = None, sA = None):
  #Returns the ratio between the negative lobe samples of 1x -> 4x
  avg_1x = getMatrixAvg(np.add(pb1, -s1), 700, 800, 0, len(pb1))
  avg_4x = getMatrixAvg(np.add(pb4, -s4), 700, 800, 0, len(pb4))
  gc = (avg_4x/avg_1x)

  if pbA is not None and sA is not None: #Applying minor 4x->AG gain correction
    avg_AG = getMatrixAvg(np.add(pbA, -sA), 700, 800, 0, len(pbA))
    gc = gc*(avg_AG/avg_4x)
    print("Gain Constant: " + str(np.round(gc, 2)) + ", uncorrected: "+ str(np.round(avg_4x/avg_1x, 2)))
  else:
    print("Gain Constant: " + str(np.round(gc, 2)))
  if(np.abs(gc - 4.) > 0.5):
    print("!! Gain constant far from 4.0! Forcing it to 4.0 instead")
    gc = 4.0
  return gc

#Shift and Scale Blocks. Applies 1x, 4x, and AG shifts (s1, s4, sA) to their respective samples and scales 1x samples by the gain constant(gc)
def shift_and_scale_blocks(pb, bb, sh, gc, neg=False):
  for gain in pb.keys():
    n_blocks = len(pb[gain])
    n_samples = len(pb[gain][0])
    for b in range(n_blocks): # b = block index
      for s in range(n_samples): # s = sample index
        #1x: Shift and Scale
        #if(gain == "1x"): pb[gain][b][s] = (pb[gain][b][s] - sh[gain])*gc
        if(gain == "1x"): pb[gain][b][s] = (pb[gain][b][s] - sh[gain])
        #if(gain == "1x"): pb[gain][b][s] = (pb[gain][b][s] + int(neg)*sh[gain])
        #4x: Shift
        #elif(gain == "4x"): pb[gain][b][s] = (pb[gain][b][s] - sh[gain])
        elif(gain == "4x"): pb[gain][b][s] = (pb[gain][b][s] - sh[gain])
        #AG: Apply AG shift, scale when decision bit = 1
        elif(gain == "AG"):
          if(bb[gain][b][s] == 1):
            pb[gain][b][s] = (pb[gain][b][s] - sh["1x"])*gc
          else:
            pb[gain][b][s] = (pb[gain][b][s] - sh[gain]) 
  return pb

def scale_pulse_blocks(pb, gc):
  for b in range(len(pb)): #n_blocks
    for s in range(len(pb[0])): #n_samples
      pb[b][s] = pb[b][s]/gc
  return pb 


#This aligns each finely sample pulse by minimizing the sum of difference squared of samples near the highest derivative point, then sets the average max index to 165
def align_blocks(samples, bits):

  #Shift blocks near 165  
  for b in range(len(samples)):
    max_ind = np.argmax(samples[b])
    samples[b] = np.roll(samples[b], 165-max_ind)
    bits[b] = np.roll(bits[b], 165-max_ind)

  n_SSQ = 5 #Samples to use in sum of difference squared calculation
  n_shift = 5 #Search range: How far apart the samples will be shifted from each other until the minimization ends
  mean_index = np.argmax(samples[0])
  #d_idx = np.argmax(np.gradient(samples[0])) #Index of the max derivative of the first block
  d_idx = np.argmax(samples[0]) - 30 #Index of the max derivative of the first block
  for b in range(1, len(samples)): #Iterate through blocks, but can skip the first
    index = 0
    sds = calcSumSquares(samples[0][d_idx - n_SSQ:d_idx + n_SSQ], samples[b][d_idx - n_SSQ - n_shift:d_idx + n_SSQ - n_shift])
    for shift in range(2*n_shift):
      if(calcSumSquares(samples[0][d_idx - n_SSQ:d_idx + n_SSQ], samples[b][d_idx - n_SSQ - n_shift + shift:d_idx + n_SSQ - n_shift + shift]) < sds):
        index = shift
        sds = calcSumSquares(samples[0][d_idx - n_SSQ:d_idx + n_SSQ], samples[b][d_idx - n_SSQ - n_shift + shift:d_idx + n_SSQ - n_shift + shift])
    samples[b] = np.roll(samples[b], n_shift - index)
    bits[b] = np.roll(bits[b], n_shift - index)
    mean_index += np.argmax(samples[b])

  #Set the peak at index 165, if it isn't already
  #mean_index = int(np.round(mean_index/len(samples)))
  mean_index = np.argmax(np.mean(samples, axis=0))
  if(mean_index != 165):
    for b in range(len(samples)):
      samples[b] = np.roll(samples[b], 165 - mean_index)
      bits[b] = np.roll(bits[b], 165 - mean_index)
     

  return samples, bits

def runAnalysis(data, adc, channel):
  
  #Steps:
  #1: Generate interlaced blocks ONCE
  #2: Find the 1x, 4x, AG shifts
  #3: Find the gain constant
  #     -Shift and align the 1x and 4x samples
  #     -Calculate the ratio between them
  #4: Apply the gain constant
  #     -Multiply the 1x samples by the gain constant
  #5: Align 1x, 4x, and AG samples
  #6: Save the blocks 


  global samples_to_save
  global bits_to_save 
  global gain_exist

  pulse_blocks = {}
  bit_blocks = {}
  shifts = {}

  available_gains = data.keys()
  #for gain in ("1x", "4x", "AG"):
  for gain in available_gains:

    #1: Interlace the blocks
    tmp_samples, tmp_bits = get_blocks(data[gain][0], data[gain][1])
    pulse_blocks.update({gain:tmp_samples})
    bit_blocks.update({gain:tmp_bits})

    if(debug_plots):
      plt.plot(pulse_blocks[gain][0])
      plt.title("Sample Interlaced Pulse")
      plt.show()
    
    #2: Find the 1x,4x,AG shifts by taking the mean of the first 5 samples in all files
    shifts.update({gain:getMatrixAvg(pulse_blocks[gain], 0, 5, 0, len(pulse_blocks[gain]))})
    print("Shift "+gain+": " + str(np.round(shifts[gain], 1)))

  #3: If the 4x samples exist and are not saturated, align the 1x and 4x pulses and calculate the gain ratio. Otherwise, use 4.0
  if not "4x" in available_gains or not "1x" in available_gains or np.min(pulse_blocks["4x"]) < 20: #Need both 1x and 4x. Last condition is to check for saturation
    gain_const = 4.0
  else:
    gain_const = get_gain_constant(pulse_blocks["1x"], pulse_blocks["4x"], shifts["1x"], shifts["4x"], pulse_blocks["AG"], shifts["AG"]) #Supply the AG arguments to calculate the corrected (small 4x->AG adjustment) gain constant
    
  #4: Apply the gain constant and shift blocks
  pulse_blocks = shift_and_scale_blocks(pulse_blocks, bit_blocks, shifts, gain_const)

  #5: Align the blocks 
  for gain in available_gains:
    pulse_blocks[gain], bit_blocks[gain] = align_blocks(pulse_blocks[gain], bit_blocks[gain])
    samples_to_save[gain][adc][channel] = pulse_blocks[gain]
    bits_to_save[gain][adc][channel] = bit_blocks[gain]

  #pulse_blocks = shift_and_scale_blocks(pulse_blocks, bit_blocks, shifts, gain_const, True)
  #for gain in available_gains:
  #  samples_to_save[gain][adc][channel] = pulse_blocks[gain]
  #  bits_to_save[gain][adc][channel] = bit_blocks[gain]

def save_blocks(f, s, b):
  global n_channels
  for adc in range(len(s)):
    for channel in range(n_channels):
      f.create_dataset("coluta"+str(adc+1)+"/channel"+str(channel+1)+"/samples", data=s[adc][channel])
      f.create_dataset("coluta"+str(adc+1)+"/channel"+str(channel+1)+"/bits", data=b[adc][channel])
      #f["/adc"+str(adc)].create_dataset("samples", data = s[adc])
      #f["/adc"+str(adc)].create_dataset("bits", data = b[adc])
  f.close()

def get_data_from_hdf5(f, n_measurements):
  global n_adcs
  global n_channels
  global n_samples_per_run
  global pulse_length
  data = np.empty(shape=(n_adcs, n_channels, 2, n_measurements, n_samples_per_run))
  extra_sample_offset = 2
  n_pulse_skips = 12
  for adc in range(n_adcs):
    for channel in range(n_channels):
      if(f["Measurement_0/coluta"+str(adc+1)].attrs["channels"][channel] == 'frame'): continue
      for measurement in range(n_measurements):
        #data[adc][channel][0][measurement] = f["Measurement_"+str(measurement)+"/coluta"+str(adc+1)+"/channel"+str(channel+1)+"/samples"][()][3*pulse_length:3*pulse_length+n_samples_per_run]
        #data[adc][channel][1][measurement] = f["Measurement_"+str(measurement)+"/coluta"+str(adc+1)+"/channel"+str(channel+1)+"/bits"][()][3*pulse_length:3*pulse_length+n_samples_per_run].astype(int)
        data[adc][channel][0][measurement] = f["Measurement_"+str(measurement)+"/coluta"+str(adc+1)+"/channel"+str(channel+1)+"/samples"][()][n_pulse_skips*pulse_length - extra_sample_offset:n_pulse_skips*pulse_length - extra_sample_offset + n_samples_per_run]
        data[adc][channel][1][measurement] = f["Measurement_"+str(measurement)+"/coluta"+str(adc+1)+"/channel"+str(channel+1)+"/bits"][()][n_pulse_skips*pulse_length - extra_sample_offset:n_pulse_skips*pulse_length - extra_sample_offset + n_samples_per_run].astype(int)
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
      print("-%s" % short(run))
    continue_script = input("\nContinue running? [y/n] (default=y): ")
    if(continue_script == "n"):
      sys.exit(0)
     

  #Run over data found in each run's directory
  for run in runs:
    amps = getAmps(run)
    #n_files = getNFiles(run+"Pulse_Amp"+amps[0]+"_1x") #TODO Less hard coding
    tmpFile = h5py.File(glob.glob(run+"Pulse_Amp"+amps[0]+"*.hdf5")[0], 'r') #TODO Less hard coding

    #Collect data from attributes
    n_adcs = tmpFile.attrs["n_adcs"]
    #n_channels = tmpFile.attrs["n_channels"]
    n_channels = 2 
    #n_samples_per_run = tmpFile.attrs["n_samples_per_run"]
    n_pulses = tmpFile.attrs["n_pulses"]
    pulse_length = tmpFile.attrs["n_samples_per_pulse"]
    n_samples_per_run = n_pulses*pulse_length
    #n_offset = tmpFile.attrs["n_samples_offset"]
    n_offset = 1
    adc_freq = tmpFile.attrs["adc_freq"]
    awg_freq = tmpFile.attrs["awg_freq"]
    #block_length = int(((awg_freq/np.abs(n_offset))/adc_freq)*pulse_length)
    block_length = int(((awg_freq/n_offset)/adc_freq)*pulse_length)
    n_blocks = int((n_samples_per_run/block_length))
    #n_pulses = int(n_samples_per_run/pulse_length)
    


    n_skips = 1 #Of files skipped (skipping first file for Feb6 data)
    
    print("=================================================================")
    print("= Processing "+str(len(amps))+" amplitudes")
    print("= Amplitudes: "+str(amps).replace("p","."))
    print("= Setup: "+str(n_adcs)+" ADCs per board, "+str(n_channels)+" channels per ADC")
    saveDir = checkDir(getRootDir()+"/Physics/Blocks_Final/"+short(run))

    for amp in amps:      
      gain_exist = {"1x":False, "4x":False, "AG":False} #Dictionary which tracks which gains exist
      if not glob.glob(run+"*Amp"+amp+"*1x*"):
        print("1x Data not found for amplitude "+amp.replace("p", "."))
      else:
        gain_exist["1x"] = True
      if not glob.glob(run+"*Amp"+amp+"*4x*"):
        print("4x Data not found for amplitude "+amp.replace("p", "."))
      else:
        gain_exist["4x"] = True
      if not glob.glob(run+"*Amp"+amp+"*AG*"):
        print("AG Data not found for amplitude "+amp.replace("p", "."))
      else:
        gain_exist["AG"] = True

      data_full = {}
      samples_to_save = {}
      bits_to_save = {}
      save_files = {}
      for gain in ("1x", "4x", "AG"):
        if not gain_exist[gain]: continue 
        tmpfile = h5py.File(run+'Pulse_Amp'+amp+'_'+gain+'.hdf5', 'r')
        n_measurements = tmpfile.attrs["n_measurements"]
        tmpdata = get_data_from_hdf5(tmpfile, n_measurements)
        data_full.update({gain:tmpdata})
        samples_to_save.update({gain:np.zeros(shape=(n_adcs, n_channels, n_measurements*n_blocks, block_length))})
        bits_to_save.update({gain:np.zeros(shape=(n_adcs, n_channels, n_measurements*n_blocks, block_length))})
        save_files.update({gain:makeH5File(saveDir+'Blocks_Amp'+amp+'_'+gain+'.hdf5', 'core')})
        makeAttributes(save_files[gain])
      for adc in range(n_adcs):
        #for channel in range(1,2):
        for channel in range(n_channels):
          data_tmp = {}
          for gain in data_full.keys():
            data_tmp.update({gain:data_full[gain][adc][channel]})
          print("Running analysis for ADC "+str(adc)+", Channel "+str(channel)+", Amplitude "+str(amp)+" ...")
          #The input file has a 1D array for each channel, the output file with have a 2D array of blocks for each channel, but the file structure will remain the same
          runAnalysis(data_tmp, adc, channel) #Run the analysis, and store the pulse blocks
      #Save the aligned blocks
      for gain in data_full.keys():
        save_blocks(save_files[gain], samples_to_save[gain], bits_to_save[gain])

  print("This process took: " + str(np.round((time.time() - t_start), 3))+" seconds to complete")


