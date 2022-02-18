import numpy as np
import csv
import os
import sys
import h5py
"""
ProcessData.py: 

Input: A raw .csv file of data from COLUTAv2
Output: A 2D array of two rows and however many colums that there are measured samples from COLUTAv2
  -First row: Decision bits
  -Second row: ADC samples
The output is saved as a numpy file (.npy)
"""


###GLOBAL### TODO: Make this easier to edit per chip
#chipWeights = [2048.00, 1024.25, 512.50, 256.25, 128.00, 64.00, 130.75, 66.00, 32.75, 16.25, 8.25, 4.00, 2.00, 1.00]
chipWeights = [2046.25, 1023.00, 511.75, 255.75, 128.00, 64.00, 132.75, 67.50, 33.00, 16.75, 8.25, 4.00, 2.00, 1.00]
mean1x = 0
mean4x = 0
#gainConst = 4.01
#gainConst = 1


"""
getModifiedSample(1D Array): Return the modified, weighted sample as described by:

Check if the decision bit is 1 (1x) or 0 (4x):
- If 1(1x): modified sample = raw sample - mean1x
- If 0(4x): modified sample = (raw sample - mean4x)/gain


Input: Row of a raw CSV file, ignoring the first bit. row[0] should be the decision bit with row[1:15] equivalent to bits 3:16 in the csv file
Output: The modified ADC sample

"""
def getModifiedSample(row):
  modifiedSample = 0
  for i in range(14):
    modifiedSample += float(row[i+1])*chipWeights[i] #Get raw sample 
  return modifiedSample



def makeH5File(name, mode, dr = None):
  if os.path.exists(name):
    os.remove(name)
  if dr is not None:
    return h5py.File(name, mode, driver=dr)
  else:
    return h5py.File(name, mode)



def setHDF5Attributes(filename, n_adcs, n_channels, awg_amp, adc_freq, awg_freq, sine_freq, dc_offset, runtype, n_samples_per_run, n_files, n_skips, n_samples_per_pulse, n_samples_offset, notes):
  filename.attrs.create("n_adcs", n_adcs)
  filename.attrs.create("n_channels", n_channels)
  filename.attrs.create("awg_amp", awg_amp)
  filename.attrs.create("adc_freq", adc_freq)
  filename.attrs.create("awg_freq", awg_freq)
  filename.attrs.create("sine_freq", sine_freq)
  filename.attrs.create("dc_offset", dc_offset)
  filename.attrs.create("runtype", runtype, dtype = h5py.special_dtype(vlen=str))
  filename.attrs.create("n_samples_per_run", n_samples_per_run)
  filename.attrs.create("n_samples_per_pulse", n_samples_per_pulse)
  filename.attrs.create("n_samples_offset", n_samples_offset)
  filename.attrs.create("notes", notes, dtype = h5py.special_dtype(vlen=str))
  filename.attrs.create("n_measurements", int(n_files - n_skips))


def processFeb6():

  gitRootDir = os.path.abspath(os.path.join(sys.path[0], os.pardir))+"/" #Root directory of the repo.
  dataDirRaw = gitRootDir+"Data/Raw/" #Directory of Raw Data from COLUTAv2
  dataDirProcessed = gitRootDir+"Data/Processed/" #Output directory for processed data


  gains = ["1x", "4x", "AG"]
  gainsdict = { "AG":0, "1x":1, "4x":2 }
  bLines = ["_BLShift"]
  bLine = 0
  startSample = 125 
  nSamples = 6000
  nFiles = 30
  nSkips = 1
  runNames = [
    "Run_0669" + "/",
    ]
  amps = [
    "0p0001",
    "0p001",
    "0p01",
    "0p025",
    "0p05",
    "0p075",
    "0p1",
    "0p125",
    "0p15",
    "0p175",
    "0p2",
    ]


  for runNum in range(len(runNames)): #Iterate through all runs
  #for amp in range(len(amps)): #Iterate through all runs
    filenames = [] #List of filenames to be processed
    for fname in os.listdir(dataDirRaw+runNames[runNum]):
      if ".csv" not in fname and "CHANNEL1" not in fname:
        continue
      filenames.append(fname) 
    #For each file in filenames, process the file, and save it in Data/Processed
    filenames.sort()
    for amp in range(len(amps)):
      for gain in range(len(gains)):
        f = makeH5File(dataDirProcessed+"Pulses/Feb6/Pulse_Amp"+amps[amp]+"_"+gains[gain]+".hdf5", "w")
        setHDF5Attributes(f, 1, 1, float(amps[amp].replace("p",".")), 40, 1200, 0, 0, "pulse", nSamples, nFiles, nSkips, 40, -1, "")
        for i in range(nSkips, nFiles): #Loop through files, skip the first file
          processedData = np.zeros([2, nSamples]) #Initialize a 2D array of 2 rows and nSamples columns. First row will be the array of the 2nd bits in the raw data (decision bits) second row will be the weighted ADC counts
          #Open the file again and fill the data array
          with open(str(dataDirRaw+runNames[runNum]+filenames[i + nFiles*gain+amp*90])) as csvfile:
            reader = csv.reader(csvfile)
            next(reader, None) #Skip the header
            ct=0 #counter
            for row in reader:
              if ct >= startSample and ct < (nSamples + startSample):
                processedData[0][ct - startSample] = getModifiedSample(row[1::]) #Modified ADC samples
                processedData[1][ct - startSample] = row[1] #Decision bit     
              ct+=1
        
            if not os.path.exists(os.path.dirname(dataDirProcessed+"Pulses/Feb6/")):
              try:
                os.makedirs(os.path.dirname(dataDirProcessed+"Pulses/Feb6/"))
              except OSError as exc: # Guard against race condition
                if exc.errno != errno.EEXIST:
                  raise

            #Save the File
            #In numpy format
            #np.save(dataDirProcessed+"Pulses/Feb6/Pulse_Amp"+amps[amp]+"_"+gains[gain]+bLines[bLine]+"_"+str(i)+".npy", processedData) #Save the file
            #print ("Saved File: "+dataDirProcessed+"Pulses/Feb6/Pulse_Amp"+amps[amp]+"_"+gains[gain]+bLines[bLine]+"_"+str(i)+".npy. File: " + filenames[i + nFiles*gain+amp*90])        

          #In HDF5 Format
          f.create_group("measurement"+str(i - nSkips))
          for adc in range(f.attrs["n_adcs"]):
            f.create_group("measurement"+str(i-nSkips)+"/adc"+str(adc))
            #f["measurement"+str(i-nSkips)+"/adc"+str(adc)].attrs.create("gains", np.full(f.attrs["n_channels"], gainsdict[gains[gain]]), dtype=int)
            #f["measurement"+str(i-nSkips)+"/adc"+str(adc)].attrs.create("sar_weights", chipWeights) 
            f["measurement"+str(i-nSkips)+"/adc"+str(adc)].create_dataset("samples", data=processedData[0].reshape(1, -1)) 
            f["measurement"+str(i-nSkips)+"/adc"+str(adc)].create_dataset("bits", data=processedData[1].reshape(1, -1).astype(bool)) 
        print("Saved File: "+dataDirProcessed+"Pulses/Feb6/Pulse_Amp"+amps[amp]+"_"+gains[gain]+".hdf5")        
            
  runNames = [
    "Run_0670" + "/",
    ]
  amps = [
    "0p225",
    "0p25",
    "0p275",
    "0p3",
    "0p325",
    "0p35",
    "0p375",
    "0p4",
    "0p425",
    "0p45",
    "0p475",
    "0p5",
    "0p525",
    "0p55",
    "0p575",
    "0p6",
    "0p625",
    "0p65",
    "0p675",
    "0p7",
    ]

  skipFor0p6 = 0

  for runNum in range(len(runNames)): #Iterate through all runs
  #for amp in range(len(amps)): #Iterate through all runs
    filenames = [] #List of filenames to be processed
    for fname in os.listdir(dataDirRaw+runNames[runNum]):
      if ".csv" not in fname and "CHANNEL1" not in fname:
        continue
      filenames.append(fname) 
    #For each file in filenames, process the file, and save it in Data/Processed
    filenames.sort()
    for amp in range(len(amps)):
      for gain in range(len(gains)):
        f = makeH5File(dataDirProcessed+"Pulses/Feb6/Pulse_Amp"+amps[amp]+"_"+gains[gain]+".hdf5", "w")
        setHDF5Attributes(f, 1, 1, float(amps[amp].replace("p",".")), 40, 1200, 0, 0, "pulse", nSamples, nFiles, nSkips, 40, -1, "")
        for i in range(1, nFiles): #Loop through files, skip the first file
          processedData = np.zeros([2, nSamples]) #Initialize a 2D array of 2 rows and nSamples columns. First row will be the array of the 2nd bits in the raw data (decision bits) second row will be the weighted ADC counts
          #Open the file again and fill the data array
          if(amps[amp] == "0p6" and gain == 1):
            skipFor0p6 = 30
          elif(amps[amp] == "0p6" and gain == 2):
            continue
          else:
            skipFor0p6 = 0
          with open(str(dataDirRaw+runNames[runNum]+filenames[i + nFiles*gain+amp*90+skipFor0p6])) as csvfile:
            reader = csv.reader(csvfile)
            next(reader, None) #Skip the header
            ct=0 #counter
            for row in reader:
              if ct >= startSample and ct < (nSamples + startSample):
                processedData[0][ct - startSample] = getModifiedSample(row[1::]) #Modified ADC samples
                processedData[1][ct - startSample] = row[1] #Decision bit     
              ct+=1
        
            if not os.path.exists(os.path.dirname(dataDirProcessed+"Pulses/Feb6/")):
              try:
                os.makedirs(os.path.dirname(dataDirProcessed+"Pulses/Feb6/"))
              except OSError as exc: # Guard against race condition
                if exc.errno != errno.EEXIST:
                  raise
            #np.save(dataDirProcessed+"Pulses/Feb6/Pulse_Amp"+amps[amp]+"_"+gains[gain]+bLines[bLine]+"_"+str(i)+".npy", processedData) #Save the file
            #print ("Saved File: "+dataDirProcessed+"Pulses/Feb6/Pulse_Amp"+amps[amp]+"_"+gains[gain]+bLines[bLine]+"_"+str(i)+".npy. File: " + filenames[i + nFiles*gain+amp*90+skipFor0p6])       
            #In HDF5 Format
          f.create_group("measurement"+str(i - nSkips))
          for adc in range(f.attrs["n_adcs"]):
            f.create_group("measurement"+str(i-nSkips)+"/adc"+str(adc))
            #f["measurement"+str(i-nSkips)+"/adc"+str(adc)].attrs.create("gains", np.full(f.attrs["n_channels"], gainsdict[gains[gain]]), dtype=int)
            #f["measurement"+str(i-nSkips)+"/adc"+str(adc)].attrs.create("sar_weights", chipWeights) 
            f["measurement"+str(i-nSkips)+"/adc"+str(adc)].create_dataset("samples", data=processedData[0].reshape(1, -1)) 
            f["measurement"+str(i-nSkips)+"/adc"+str(adc)].create_dataset("bits", data=processedData[1].reshape(1, -1).astype(bool)) 
        print("Saved File: "+dataDirProcessed+"Pulses/Feb6/Pulse_Amp"+amps[amp]+"_"+gains[gain]+".hdf5")        
  runNames = [
    "Run_0671" + "/",
    ]
  amps = [
    "0p725",
    "0p75",
    "0p775",
    "0p8",
    "0p825",
    "0p85",
    "0p875",
    "0p9",
    ]


  for runNum in range(len(runNames)): #Iterate through all runs
  #for amp in range(len(amps)): #Iterate through all runs
    filenames = [] #List of filenames to be processed
    for fname in os.listdir(dataDirRaw+runNames[runNum]):
      if ".csv" not in fname and "CHANNEL1" not in fname:
        continue
      filenames.append(fname) 
    #For each file in filenames, process the file, and save it in Data/Processed
    filenames.sort()
    for amp in range(len(amps)):
      for gain in range(len(gains)):
        f = makeH5File(dataDirProcessed+"Pulses/Feb6/Pulse_Amp"+amps[amp]+"_"+gains[gain]+".hdf5", "w")
        setHDF5Attributes(f, 1, 1, float(amps[amp].replace("p",".")), 40, 1200, 0, 0, "pulse", nSamples, nFiles, nSkips, 40, -1, "")
        for i in range(1, nFiles): #Loop through files, skip the first file
          processedData = np.zeros([2, nSamples]) #Initialize a 2D array of 2 rows and nSamples columns. First row will be the array of the 2nd bits in the raw data (decision bits) second row will be the weighted ADC counts
          #Open the file again and fill the data array
          with open(str(dataDirRaw+runNames[runNum]+filenames[i + nFiles*gain+amp*90])) as csvfile:
            reader = csv.reader(csvfile)
            next(reader, None) #Skip the header
            ct=0 #counter
            for row in reader:
              if ct >= startSample and ct < (nSamples + startSample):
                processedData[0][ct - startSample] = getModifiedSample(row[1::]) #Modified ADC samples
                processedData[1][ct - startSample] = row[1] #Decision bit     
              ct+=1
        
            if not os.path.exists(os.path.dirname(dataDirProcessed+"Pulses/Feb6/")):
              try:
                os.makedirs(os.path.dirname(dataDirProcessed+"Pulses/Feb6/"))
              except OSError as exc: # Guard against race condition
                if exc.errno != errno.EEXIST:
                  raise
            #np.save(dataDirProcessed+"Pulses/Feb6/Pulse_Amp"+amps[amp]+"_"+gains[gain]+bLines[bLine]+"_"+str(i)+".npy", processedData) #Save the file
            #print ("Saved File: "+dataDirProcessed+"Pulses/Feb6/Pulse_Amp"+amps[amp]+"_"+gains[gain]+bLines[bLine]+"_"+str(i)+".npy. File: " + filenames[i + nFiles*gain+amp*90])        
            #In HDF5 Format
          f.create_group("measurement"+str(i - nSkips))
          for adc in range(f.attrs["n_adcs"]):
            f.create_group("measurement"+str(i-nSkips)+"/adc"+str(adc))
            #f["measurement"+str(i-nSkips)+"/adc"+str(adc)].attrs.create("gains", np.full(f.attrs["n_channels"], gainsdict[gains[gain]]), dtype=int)
            #f["measurement"+str(i-nSkips)+"/adc"+str(adc)].attrs.create("sar_weights", chipWeights) 
            f["measurement"+str(i-nSkips)+"/adc"+str(adc)].create_dataset("samples", data=processedData[0].reshape(1, -1)) 
            f["measurement"+str(i-nSkips)+"/adc"+str(adc)].create_dataset("bits", data=processedData[1].reshape(1, -1).astype(bool)) 
        print("Saved File: "+dataDirProcessed+"Pulses/Feb6/Pulse_Amp"+amps[amp]+"_"+gains[gain]+".hdf5")        



  #Remove bad files
  for g in ("1x", "4x", "AG"):
    os.remove(dataDirProcessed+"Pulses/Feb6/Pulse_Amp0p025_"+g+".hdf5")
    os.remove(dataDirProcessed+"Pulses/Feb6/Pulse_Amp0p6_"+g+".hdf5")

















def process40MSPS():

  gitRootDir = os.path.abspath(os.path.join(sys.path[0], os.pardir))+"/" #Root directory of the repo.
  dataDirRaw = gitRootDir+"Data/Raw/" #Directory of Raw Data from COLUTAv2
  dataDirProcessed = gitRootDir+"Data/Processed/" #Output directory for processed data
  #40MSPS
  runNames = [
    "Run_0436" + "/",
    "Run_0437" + "/",
    "Run_0438" + "/",
    "Run_0439" + "/",
    "Run_0440" + "/",
    "Run_0441" + "/",
    "Run_0442" + "/",
    "Run_0443" + "/",
    "Run_0444" + "/",
    "Run_0445" + "/",
    ]
  #40MSPS
  amps = [
    "0p9",
    "0p8",
    "0p5",
    "0p22",
    "0p2",
    "0p18",
    "0p1",
    "0p05",
    "0p01",
    "0p001",
    ]


  gains = ["1x", "4x", "AG"]
  bLines = ["", "_BLShift"]#40MSPS
  #bLines = ["_BLShift"]#10MSPS
  startSample = 130 
  nSamples = 6000
  nFiles = 10
  for runNum in range(len(runNames)): #Iterate through all runs
  #for amp in range(len(amps)): #Iterate through all runs
    filenames = [] #List of filenames to be processed
    for fname in os.listdir(dataDirRaw+runNames[runNum]):
      if ".csv" not in fname and "CHANNEL1" not in fname:
        continue
      filenames.append(fname) 
    #For each file in filenames, process the file, and save it in Data/Processed
    filenames.sort()
    for bLine in range(len(bLines)):
      for gain in range(len(gains)):
        for i in range(1, nFiles): #Loop through files, skip the first file
          processedData = np.zeros([2, nSamples]) #Initialize a 2D array of 2 rows and nSamples columns. First row will be the array of the 2nd bits in the raw data (decision bits) second row will be the weighted ADC counts
          #Open the file again and fill the data array
          with open(str(dataDirRaw+runNames[runNum]+filenames[i + nFiles*gain+bLine*30])) as csvfile:
            reader = csv.reader(csvfile)
            next(reader, None) #Skip the header
            ct=0 #counter
            for row in reader:
              if ct >= startSample and ct < (nSamples + startSample):
                processedData[0][ct - startSample] = getModifiedSample(row[1::]) #Modified ADC samples
                processedData[1][ct - startSample] = row[1] #Decision bit     
              ct+=1
        
            if not os.path.exists(os.path.dirname(dataDirProcessed+"Pulses/40MSPS/")):
              try:
                os.makedirs(os.path.dirname(dataDirProcessed+"Pulses/40MSPS/"))
              except OSError as exc: # Guard against race condition
                if exc.errno != errno.EEXIST:
                  raise
            np.save(dataDirProcessed+"Pulses/40MSPS/Pulse_Amp"+amps[runNum]+"_"+gains[gain]+bLines[bLine]+"_"+str(i)+".npy", processedData) #Save the file
            print("Saved File: "+dataDirProcessed+"Pulses/40MSPS/Pulse_Amp"+amps[runNum]+"_"+gains[gain]+bLines[bLine]+"_"+str(i)+".npy. File: " + filenames[i + nFiles*gain+bLine*30]        )

def process10MSPS():

  gitRootDir = os.path.abspath(os.path.join(sys.path[0], os.pardir))+"/" #Root directory of the repo. TODO: Automate

  dataDirRaw = gitRootDir+"Data/Raw/" #Directory of Raw Data from COLUTAv2
  dataDirProcessed = gitRootDir+"Data/Processed/" #Output directory for processed data
  runNames = ["Run_0469/"] #10MSPS Samples
  #10MSPS
  amps = [
    "0p15",
    "0p125",
    "0p1",
    "0p075",
    "0p05",
    "0p025",
    "0p01",
    ]
  ampsHigh = [
    "0p9",
    "0p8",
    "0p5"
    ]
  



  startSample = 118 
  nSamples = 6000
  nFiles = 30


  #LOW AMP 1-30: 1x, 31-60: 4x, 61-90: AG
  gains = ["1x", "4x", "AG"]
  #bLines = ["", "_BLShift"]#40MSPS
  bLines = ["_BLShift"]#10MSPS
  #for runNum in range(len(runNames)): #Iterate through all runs
  filenames = [] #List of filenames to be processed
  for fname in os.listdir(dataDirRaw+runNames[0]):
    if ".csv" not in fname and "CHANNEL1" not in fname:
      continue
    filenames.append(fname) 
  #For each file in filenames, process the file, and save it in Data/Processed
  filenames.sort()
  for amp in range(len(amps)): #Iterate through all amps
    for bLine in range(len(bLines)):
      for gain in range(len(gains)):
        for i in range(1, nFiles): #Loop through files, skip the first file
          processedData = np.zeros([2, nSamples]) #Initialize a 2D array of 2 rows and nSamples columns. First row will be the array of the 2nd bits in the raw data (decision bits) second row will be the weighted ADC counts
          #Open the file again and fill the data array
          with open(str(dataDirRaw+runNames[0]+filenames[amp*len(gains)*nFiles + nFiles*gain + i])) as csvfile:
            reader = csv.reader(csvfile)
            next(reader, None) #Skip the header
            ct=0 #counter
            for row in reader:
              if ct >= startSample and ct < (nSamples + startSample):
                processedData[0][ct - startSample] = getModifiedSample(row[1::]) #Modified ADC samples
                processedData[1][ct - startSample] = row[1] #Decision bit     
              ct+=1
            
        
            if not os.path.exists(os.path.dirname(dataDirProcessed+"Pulses/10MSPS/")):
              try:
                os.makedirs(os.path.dirname(dataDirProcessed+"Pulses/10MSPS/"))
              except OSError as exc: # Guard against race condition
                if exc.errno != errno.EEXIST:
                  raise
            np.save(dataDirProcessed+"Pulses/10MSPS/Pulse_Amp"+amps[amp]+"_"+gains[gain]+bLines[bLine]+"_"+str(i)+".npy", processedData) #Save the file
            print("Saved File: "+dataDirProcessed+"Pulses/10MSPS/Pulse_Amp"+amps[amp]+"_"+gains[gain]+bLines[bLine]+"_"+str(i)+".npy. File: " + filenames[amp*len(gains)*nFiles + nFiles*gain + i]        )



  #HIGH AMP 691-720: 1x, 721-750: AG
  gains = ["1x", "AG"]
  #bLines = ["", "_BLShift"]#40MSPS
  bLines = ["_BLShift"]#10MSPS
  #for runNum in range(len(runNames)): #Iterate through all runs
  filenames = [] #List of filenames to be processed
  for fname in os.listdir(dataDirRaw+runNames[0]):
    if ".csv" not in fname and "CHANNEL1" not in fname:
      continue
    filenames.append(fname) 
  #For each file in filenames, process the file, and save it in Data/Processed
  filenames.sort()
  for amp in range(1, len(ampsHigh)): #Iterate through all amps
    for bLine in range(len(bLines)):
      for gain in range(len(gains)):
        for i in range(1, nFiles): #Loop through files, skip the first file
          processedData = np.zeros([2, nSamples]) #Initialize a 2D array of 2 rows and nSamples columns. First row will be the array of the 2nd bits in the raw data (decision bits) second row will be the weighted ADC counts
          #Open the file again and fill the data array
          with open(str(dataDirRaw+runNames[0]+filenames[amp*len(gains)*nFiles + nFiles*gain + i + 630])) as csvfile:
            reader = csv.reader(csvfile)
            next(reader, None) #Skip the header
            ct=0 #counter
            for row in reader:
              if ct >= startSample and ct < (nSamples + startSample):
                processedData[0][ct - startSample] = getModifiedSample(row[1::]) #Modified ADC samples
                processedData[1][ct - startSample] = row[1] #Decision bit     
              ct+=1
            
        
            if not os.path.exists(os.path.dirname(dataDirProcessed+"Pulses/10MSPS/")):
              try:
                os.makedirs(os.path.dirname(dataDirProcessed+"Pulses/10MSPS/"))
              except OSError as exc: # Guard against race condition
                if exc.errno != errno.EEXIST:
                  raise
            np.save(dataDirProcessed+"Pulses/10MSPS/Pulse_Amp"+ampsHigh[amp]+"_"+gains[gain]+bLines[bLine]+"_"+str(i)+".npy", processedData) #Save the file
            print("Saved File: "+dataDirProcessed+"Pulses/10MSPS/Pulse_Amp"+ampsHigh[amp]+"_"+gains[gain]+bLines[bLine]+"_"+str(i)+".npy. File: " + filenames[amp*len(gains)*nFiles + nFiles*gain + i + 630]        )







if __name__ == "__main__":
# process10MSPS()
# process40MSPS()
  processFeb6()
