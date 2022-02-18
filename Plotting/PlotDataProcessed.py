import numpy as np
import matplotlib.pyplot as plt

if __name__ == "__main__":

	gitRootDir = "../" #Root directory of the repo. TODO: Automate

	inputDir = gitRootDir+"Data/Processed/Run_0418/" #Output directory for processed data
	filenames = ["CHANNEL1_0008_Binary_processed.npy"] #List of filenames to be processed

	#Plot samples 3000 to 3500
	startSample = 0
	nSamplesToPlot = 4096
	print nSamplesToPlot
	for i in range(len(filenames)):

		data = np.load(inputDir+filenames[i])		
		ax = plt.subplot(211) #Top Plot: ADC Samples
		plt.title("Modified ADC Samples")
		plt.plot(np.linspace(startSample, startSample+nSamplesToPlot, len(data[0][startSample:startSample+nSamplesToPlot])), data[0][startSample:startSample+nSamplesToPlot])
		ax = plt.subplot(212) #Bottom Plot: Decision Bits
		plt.title("Decision bits")
		plt.xlabel("Filename: " + filenames[i])
		plt.plot(np.linspace(startSample, startSample+nSamplesToPlot, len(data[1][startSample:startSample+nSamplesToPlot])), data[1][startSample:startSample+nSamplesToPlot])
		plt.ylim(-0.2, 1.2)
		plt.show()
		#plt.savefig("PlotsSaved/plots_"+str(filenames[i]).replace(".npy",".png"))
		print "Saved File: PlotsSaved/plots_"+str(filenames[i]).replace(".npy", ".png")
		plt.clf()
	

