import numpy as np
import csv
import matplotlib.pyplot as plt


"""
AlignAndSubtract.py: 

This script aligns two sin waves from different datasets, and then takes the difference and plots them.

Input: A processed .npy file of COLUTA sine wave data
"""


###GLOBAL### 






def calcSumSquares(list1, list2):
	sumSquares = 0
	for i in range(len(list1)):
	  sumSquares += np.square(list1[i] - list2[i])
	return sumSquares





if __name__ == "__main__":

	gitRootDir = "../" #Root directory of the repo. TODO: Automate

	dataDirProcessed = gitRootDir+"Data/Processed/" #Output directory for processed data
	filenames = ["Sine_0p2MHz/Averaged/379_Amp0p22_1x_Averaged.npy", "Sine_0p2MHz/Averaged/379_Amp0p22_AG_Averaged.npy"] #List of filenames to be processed. Exclusively two are required
	#filenames = ["Sine_0p2MHz/Averaged/Amp0p22_1x_Averaged.npy", "Sine_0p2MHz/Averaged/Amp0p22_AG_Averaged.npy"] #List of filenames to be processed. Exclusively two are required
	#filenames = ["Sine_0p2MHz/Amp0p22_1x_1.npy", "Sine_0p2MHz/Amp0p22_AG_1.npy"] #List of filenames to be processed. Exclusively two are required

		

	data1x = np.load(dataDirProcessed+filenames[0])
	dataAG = np.load(dataDirProcessed+filenames[1])

	startSample = 50
	nSamplesToPlot = 512

	print len(data1x[0])
	"""
	#Find the max within the first oscillation period and roll the array so that it becomes the 10th element
	max_idx = np.argmax(data1x[0][0:50])
	data1x[0] = np.roll(data1x[0], 10 - max_idx)
	data1x[1] = np.roll(data1x[1], 10 - max_idx)
	max_idx = np.argmax(dataAG[0][0:50])
	dataAG[0] = np.roll(dataAG[0], 10 - max_idx)
	dataAG[1] = np.roll(dataAG[1], 10 - max_idx)

	#Compare the two arrays and align them by minimizing the sum of differences squares around the first peak
	sumSquares = calcSumSquares(data1x[0][5:15], dataAG[0][5:15])
	index = 0
	for i in range(0,10):
		if(calcSumSquares(data1x[0][5:15], dataAG[0][i:i+10]) < sumSquares):
			index = i
			sumSquares = calcSumSquares(data1x[0][5:15], dataAG[0][i:i+10]) 
	dataAG[0] = np.roll(dataAG[0], 5-index)
	dataAG[1] = np.roll(dataAG[1], 5-index)
	"""
#	plt.plot(np.linspace(5, 54, 50), data1x[0][5:55])
#	plt.plot(np.linspace(5, 54, 50), dataAG[0][5:55])
#	plt.show()
	
	subtracted = np.subtract(dataAG[0], data1x[0])
	ax = plt.subplot(311) #Top Plot: ADC Samples
	plt.title("ADC Samples")
	plt.plot(np.linspace(startSample, startSample+nSamplesToPlot, len(dataAG[0][startSample:startSample+nSamplesToPlot])), dataAG[0][startSample:startSample+nSamplesToPlot], '.-')
	plt.plot(np.linspace(startSample, startSample+nSamplesToPlot, len(data1x[0][startSample:startSample+nSamplesToPlot])), data1x[0][startSample:startSample+nSamplesToPlot], '.-')
	ax = plt.subplot(312) #Top Plot: ADC Samples
	plt.title("Subtracted ADC Samples")
	plt.plot(np.linspace(startSample, startSample+nSamplesToPlot, len(subtracted[startSample:startSample+nSamplesToPlot])), subtracted[startSample:startSample+nSamplesToPlot], '.-')
	ax = plt.subplot(313) #Bottom Plot: Decision Bits
	plt.title("Decision bits")
	plt.plot(np.linspace(startSample, startSample+nSamplesToPlot, len(dataAG[1][startSample:startSample+nSamplesToPlot])), dataAG[1][startSample:startSample+nSamplesToPlot], '.-')
	plt.ylim(-0.2, 1.2)
	plt.subplots_adjust(hspace=0.4)
	plt.xlabel("Comparing: "+filenames[0].replace("/Averaged/", ", ").replace("_Averaged.npy", " (Averaged)").replace("Sine_","").replace("0p","0.") + " with " + filenames[1].replace("/Averaged/", ", ").replace("_Averaged.npy", " (Averaged)").replace("Sine_","").replace("0p","0."))
	#plt.show()
	plt.savefig("../Plotting/PlotsSaved/comparison_Run379.png")
	print "Saved File: ../Plotting/PlotsSaved/comparison_Run379.png"
	plt.clf()
	


