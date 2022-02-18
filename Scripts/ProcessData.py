import numpy as np
import csv
import os
import sys

"""
ProcessData.py: 

Input: A raw .csv file of data from COLUTAv2
Output: A 2D array of two rows and however many colums that there are measured samples from COLUTAv2
	-First row: Decision bits
	-Second row: ADC samples
The output is saved as a numpy file (.npy)
"""


###GLOBAL### TODO: Make this easier to edit per chip
chipWeights = [2048.00, 1024.25, 512.50, 256.25, 128.00, 64.00, 130.75, 66.00, 32.75, 16.25, 8.25, 4.00, 2.00, 1.00]
mean1x = 2149.11
mean4x = 2149.14
gainConst = 4.0655


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
	if int(row[0]) == 1:
		modifiedSample = modifiedSample - mean1x
	elif int(row[0]) == 0:
		modifiedSample = (modifiedSample - mean4x)/gainConst
	return modifiedSample





if __name__ == "__main__":

	gitRootDir = os.path.abspath(os.path.join(sys.path[0], os.pardir))+"/" #Root directory of the repo.
	dataDirRaw = gitRootDir+"Data/Raw/" #Directory of Raw Data from COLUTAv2
	dataDirProcessed = gitRootDir+"Data/Processed/" #Output directory for processed data
	runNames = [
		"Run_0379" + "/",
#		"Run_0418" + "/",
#		"Run_0419" + "/",
#		"Run_0420" + "/",
#		"Run_0421" + "/",
#		"Run_0422" + "/",
#		"Run_0423" + "/",
#		"Run_0424" + "/",
#		"Run_0425" + "/",
#		"Run_0426" + "/",
#		"Run_0427" + "/",
#		"Run_0428" + "/",
#		"Run_0429" + "/",
#		"Run_0430" + "/",
#		"Run_0431" + "/",
#		"Run_0432" + "/",
#		"Run_0433" + "/",
		]

	frequencies = [
		"Sine_0p2MHz" + "/",
#		"Sine_1MHz" + "/",
#		"Sine_2MHz" + "/",
#		"Sine_5MHz" + "/",
#		"Sine_8MHz" + "/",
#		"Sine_10p6MHz" + "/",
#		"Sine_13MHz" + "/",
#		"Sine_18MHz" + "/",
#		"Sine_19p5MHz" + "/",
#		"Sine_20p5MHz" + "/",
#		"Sine_22MHz" + "/",
#		"Sine_29p4MHz" + "/",
#		"Sine_32MHz" + "/",
#		"Sine_35MHz" + "/",
#		"Sine_39MHz" + "/",
#		"Sine_39p8MHz" + "/",
		]



	gains0p22 = ["1x", "4x", "AG"]
	gains = ["1x", "AG"]
	amps = ["0p5", "0p9"]
	#amps = []
	startSample = 58
	nSamples = 4096
	for runNum in range(len(runNames)): #Iterate through all runs
		filenames = [] #List of filenames to be processed
		for fname in os.listdir(dataDirRaw+runNames[runNum]):
			if ".csv" not in fname and "CHANNEL1" not in fname:
				continue
			filenames.append(fname)	
			#print fname		
		#For each file in filenames, process the file, and save it in Data/Processed
		filenames.sort()
		#for i in range(len(filenames)):
		for gain in range(3):
			averagedData = np.zeros([2,nSamples]) 
			for i in range(1, 30): #Loop through 0.22V files, skip the first file
				"""
				#Open the file and count the lines TODO: Make this not require another iterator	
				with open(dataDirRaw+runNames[runNum]+str(filenames[i])) as csvfile:
					lineCounter = csv.reader(csvfile)
					next(lineCounter, None) #Skip the header
					nSamples = sum(1 for row in lineCounter) #Sum up the number of lines
				"""
				processedData = np.zeros([2, nSamples]) #Initialize a 2D array of 2 rows and nSamples columns. First row will be the array of the 2nd bits in the raw data (decision bits) second row will be the weighted ADC counts
				#processedData = np.empty([2, 4096]) #Initialize a 2D array of 2 rows and nSamples columns. First row will be the array of the 2nd bits in the raw data (decision bits) second row will be the weighted ADC counts
				


				#Open the file again and fill the data array
				with open(str(dataDirRaw+runNames[runNum]+filenames[i + 30*gain])) as csvfile:
					reader = csv.reader(csvfile)
					next(reader, None) #Skip the header
					ct=0 #counter
					for row in reader:
						if ct >= startSample and ct < (nSamples + startSample):
							processedData[0][ct - startSample] = getModifiedSample(row[1::]) #Modified ADC samples
							processedData[1][ct - startSample] = row[1]	#Decision bit			
							averagedData[0][ct - startSample] += processedData[0][ct - startSample]/29.0
							averagedData[1][ct - startSample] += processedData[1][ct - startSample]/29.0
						ct+=1
			
					if not os.path.exists(os.path.dirname(dataDirProcessed+frequencies[runNum])):
						try:
							os.makedirs(os.path.dirname(dataDirProcessed+frequencies[runNum]))
						except OSError as exc: # Guard against race condition
							if exc.errno != errno.EEXIST:
								raise
					np.save(dataDirProcessed+frequencies[runNum]+"Amp0p22_"+gains0p22[gain]+"_"+str(i)+".npy", processedData) #Save the file
					print "Saved File: "+dataDirProcessed+frequencies[runNum]+"Amp0p22_"+gains0p22[gain]+"_"+str(i)+".npy: "+filenames[30*gain + i].replace("CHANNEL1_", "").replace("_Binary", "")
				
			if not os.path.exists(os.path.dirname(dataDirProcessed+frequencies[runNum]+"Averaged/")):
				try:
					os.makedirs(os.path.dirname(dataDirProcessed+frequencies[runNum]+"Averaged/"))
				except OSError as exc: # Guard against race condition
					if exc.errno != errno.EEXIST:
						raise
			np.save(dataDirProcessed+frequencies[runNum]+"Averaged/Amp0p22_"+gains0p22[gain]+"_Averaged.npy", averagedData) #Save the file
			print "Saved File: "+dataDirProcessed+frequencies[runNum]+"Averaged/Amp0p22_"+gains0p22[gain]+"_Averaged.npy"



		for amp in range(2):
			for gain in range(2):
				averagedData = np.zeros([2,nSamples]) 
				for i in range(1, 30): #Loop through 0.22V files, skip the first file
					"""
					#Open the file and count the lines TODO: Make this not require another iterator	
					with open(dataDirRaw+runNames[runNum]+str(filenames[i])) as csvfile:
						lineCounter = csv.reader(csvfile)
						next(lineCounter, None) #Skip the header
						nSamples = sum(1 for row in lineCounter) #Sum up the number of lines
					"""
					processedData = np.zeros([2, nSamples]) #Initialize a 2D array of 2 rows and nSamples columns. First row will be the array of the 2nd bits in the raw data (decision bits) second row will be the weighted ADC counts
					#processedData = np.empty([2, 4096]) #Initialize a 2D array of 2 rows and nSamples columns. First row will be the array of the 2nd bits in the raw data (decision bits) second row will be the weighted ADC counts
					


					#Open the file again and fill the data array
					with open(str(dataDirRaw+runNames[runNum]+filenames[90 + 60*amp + 30*gain + i])) as csvfile:
						reader = csv.reader(csvfile)
						next(reader, None) #Skip the header
						ct=0 #counter
						for row in reader:
							if ct >= startSample and ct < (nSamples + startSample):
								processedData[0][ct - startSample] = getModifiedSample(row[1::]) #Modified ADC samples
								processedData[1][ct - startSample] = row[1]	#Decision bit			
								averagedData[0][ct - startSample] += processedData[0][ct - startSample]/29.0
								averagedData[1][ct - startSample] += processedData[1][ct - startSample]/29.0
							ct+=1
				
						if not os.path.exists(os.path.dirname(dataDirProcessed+frequencies[runNum])):
							try:
								os.makedirs(os.path.dirname(dataDirProcessed+frequencies[runNum]))
							except OSError as exc: # Guard against race condition
								if exc.errno != errno.EEXIST:
									raise
						np.save(dataDirProcessed+frequencies[runNum]+"Amp"+amps[amp]+"_"+gains[gain]+"_"+str(i)+".npy", processedData) #Save the file
						print "Saved File: "+dataDirProcessed+frequencies[runNum]+"Amp"+amps[amp]+"_"+gains[gain]+"_"+str(i)+".npy: "+filenames[90 + 60*amp + 30*gain + i].replace("CHANNEL1_", "").replace("_Binary", "")
				if not os.path.exists(os.path.dirname(dataDirProcessed+frequencies[runNum]+"Averaged/")):
					try:
						os.makedirs(os.path.dirname(dataDirProcessed+frequencies[runNum]+"Averaged/"))
					except OSError as exc: # Guard against race condition
						if exc.errno != errno.EEXIST:
							raise
				np.save(dataDirProcessed+frequencies[runNum]+"Averaged/Amp"+amps[amp]+"_"+gains[gain]+"_Averaged.npy", averagedData) #Save the file
				print "Saved File: "+dataDirProcessed+frequencies[runNum]+"Averaged/Amp"+amps[amp]+"_"+gains[gain]+"_Averaged.npy"



