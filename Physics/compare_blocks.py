import filters_functions
from filters_functions import *
import numpy as np
import matplotlib.pyplot as plt
import os
import sys

from helperfunctions import *


#### GLOBAL VARIABLES #############################################
_print_ = 0
_show_all_figs_ = 0

exclude_n_points = 800 #initial points to exclude

n_blocks = 10 #18
n_pulses_per_block = 30 #keep
#n_points_per_pulse = 400 #240
n_points_per_pulse = 40 #240
n_points_per_block = n_points_per_pulse*n_pulses_per_block #6000

n_pulses = n_blocks*n_pulses_per_block
n_points = n_pulses*n_points_per_pulse #points read in

colpulse = 2
colnoise = 1

ratio = 0


dataDir = sys.path[0]+"/Blocks_Final/"
plotsDir = os.path.abspath(os.path.join(sys.path[0], os.pardir))+"/Plotting/PlotsSaved/Pulses/Compare/"

amps = []
gains = []

versions = getVersions()

#dataDir = 'Blocks_Final/10MSPS/' 
###################################################################

#Plot Legend
legend_names = []
plot_style = ['b.-', 'c.-', 'k.-', 'r.-', 'm.-', 'g.-', 'b.-']

#Run:
amp_calc = []
amp_func = []
pulse_stdev_master_avg = numpy.empty(n_points_per_block)

n_inputs = len(sys.argv) - 1


#Input: python compare_two.py amp first second shift
def compareTwo(version):
	
	plotsDir = checkDir(os.path.abspath(os.path.join(sys.path[0], os.pardir))+"/Plotting/PlotsSaved/Pulses/Compare/"+version+"/")

	for amp in range(len(amps)):
		has4x = os.path.isfile(dataDir+version+"/Pulse_Amp"+amps[amp]+"_4x.npy")
		if(has4x):
			gains = ["1x", "4x", "AG"]
		else:
			gains = ["1x", "AG"]
		for gain1 in range(len(gains)-1):
			for gain2 in range(gain1+1, len(gains)):
				name1 = dataDir+version+'/Pulse_Amp'+amps[amp]+'_'+gains[gain1]+'.npy'
				name2 = dataDir+version+'/Pulse_Amp'+amps[amp]+'_'+gains[gain2]+'.npy'
				name1Bits = dataDir+version+'/Bits_Amp'+amps[amp]+'_'+gains[gain1]+'.npy'
				name2Bits = dataDir+version+'/Bits_Amp'+amps[amp]+'_'+gains[gain2]+'.npy'
				
				print('Reading in ' + name1 )
				tmp_blocks1 = np.load(name1)
				gc = tmp_blocks1[len(tmp_blocks1)-1][0]
				tmp_blocks1 = np.delete(tmp_blocks1, len(tmp_blocks1)-1, 0)
				#pulses_blocks1 = np.transpose(numpy.load(name1))
				pulses_blocks1 = np.transpose(tmp_blocks1)
				bits_blocks1 = np.transpose(numpy.load(name1Bits))
				print('Reading in ' + name2 )
				tmp_blocks2 = np.load(name2)
				tmp_blocks2 = np.delete(tmp_blocks2, len(tmp_blocks2)-1, 0)
				#pulses_blocks2 = np.transpose(numpy.load(name2))
				pulses_blocks2 = np.transpose(tmp_blocks2)
				bits_blocks2 = np.transpose(numpy.load(name2Bits))
				
				#pulses_blocks = numpy.load(sys.argv[i_f+1])
				
				#Get all pulses overlayed ("Master")
				pulse_overlayed_master1 = numpy.empty(len(pulses_blocks1))
				pulse_overlayed_master1Bits = numpy.empty(len(pulses_blocks1))
				pulse_stdev_master1 = numpy.empty(len(pulses_blocks1))
				pulse_overlayed_master2 = numpy.empty(len(pulses_blocks2))
				pulse_overlayed_master2Bits = numpy.empty(len(pulses_blocks2))
				pulse_stdev_master2 = numpy.empty(len(pulses_blocks2))
				
				
				for i_s2 in xrange(n_points_per_block):
				    pulse_overlayed_master1[i_s2] = numpy.mean(pulses_blocks1[i_s2])                      
				    pulse_overlayed_master1Bits[i_s2] = numpy.round(numpy.mean(bits_blocks1[i_s2]))                      
				    pulse_stdev_master1[i_s2] = numpy.std(pulses_blocks1[i_s2])
				    pulse_overlayed_master2[i_s2] = numpy.mean(pulses_blocks2[i_s2])                      
				    pulse_overlayed_master2Bits[i_s2] = numpy.round(numpy.mean(bits_blocks2[i_s2]))                      
				    pulse_stdev_master2[i_s2] = numpy.std(pulses_blocks2[i_s2])
				
				
				
#				g_overlayed_master1 = numpy.gradient(pulse_overlayed_master1, 1)
#				g_overlayed_master2 = numpy.gradient(pulse_overlayed_master2, 1)
				
				
				#ALIGN PULSES
				nSSQ = 10
				pulse_overlayed_master2 = np.roll(pulse_overlayed_master2, -10)
				pulse_overlayed_master2Bits = np.roll(pulse_overlayed_master2Bits, -10)
				index = 0
				dss = calcSumSquares(pulse_overlayed_master1[np.argmax(np.gradient(pulse_overlayed_master1)) - nSSQ:np.argmax(np.gradient(pulse_overlayed_master1)) + nSSQ], pulse_overlayed_master2[np.argmax(np.gradient(pulse_overlayed_master1)) - nSSQ:np.argmax(np.gradient(pulse_overlayed_master1)) + nSSQ])
				for j in range(1, 21):
					if (calcSumSquares(pulse_overlayed_master1[np.argmax(np.gradient(pulse_overlayed_master1)) - nSSQ:np.argmax(np.gradient(pulse_overlayed_master1)) + nSSQ], np.roll(pulse_overlayed_master2, j)[np.argmax(np.gradient(pulse_overlayed_master1)) - nSSQ:np.argmax(np.gradient(pulse_overlayed_master1)) + nSSQ]) < dss):
#						pulse_overlayed_master2 = np.roll(pulse_overlayed_master2, 1)
#						pulse_overlayed_master2Bits = np.roll(pulse_overlayed_master2Bits, 1)
						index = j
						dss = calcSumSquares(pulse_overlayed_master1[np.argmax(np.gradient(pulse_overlayed_master1)) - nSSQ:np.argmax(np.gradient(pulse_overlayed_master1)) + nSSQ], np.roll(pulse_overlayed_master2, j)[np.argmax(np.gradient(pulse_overlayed_master1)) - nSSQ:np.argmax(np.gradient(pulse_overlayed_master1)) + nSSQ])
				
				pulse_overlayed_master2 = np.roll(pulse_overlayed_master2, index)
				pulse_overlayed_master2Bits = np.roll(pulse_overlayed_master2Bits, index)
				#pulse_overlayed_master2 = np.roll(pulse_overlayed_master2, -1) #0p1 1x AG
				#pulse_overlayed_master2 = np.roll(pulse_overlayed_master2, -2) #0p1 4x AG
				#pulse_overlayed_master2 = np.roll(pulse_overlayed_master2, 1) #0p1 1x 4x
				#pulse_overlayed_master2 = np.roll(pulse_overlayed_master2, 0) #0p01 1x 4x
				#pulse_overlayed_master2 = np.roll(pulse_overlayed_master2, 1) #0p01 1x AG
				#pulse_overlayed_master2 = np.roll(pulse_overlayed_master2, 1) #0p01 4x AG
				#pulse_overlayed_master2 = np.roll(pulse_overlayed_master2, 1) #0p05 1x 4x
				#pulse_overlayed_master2 = np.roll(pulse_overlayed_master2, 2) #0p05 1x AG
				#pulse_overlayed_master2 = np.roll(pulse_overlayed_master2, 1) #0p05 4x AG
				
				#pulse_overlayed_master2 = np.roll(pulse_overlayed_master2, int(sys.argv[4])) 
				
				
				
				g_overlayed_master1 = numpy.gradient(pulse_overlayed_master1, 1)
				g_overlayed_master2 = numpy.gradient(pulse_overlayed_master2, 1)
				
				
				
				
				
				
				x = numpy.arange(len(pulse_overlayed_master1))
				
				difference = np.empty(len(pulse_overlayed_master1))
				bits1 = np.empty(len(pulse_overlayed_master1))
				bits2 = np.empty(len(pulse_overlayed_master1))
				g_difference = np.empty(len(pulse_overlayed_master1))
				for i in range(len(pulse_overlayed_master1)):
					difference[i] = pulse_overlayed_master1[i] - pulse_overlayed_master2[i]
					bits1[i] = pulse_overlayed_master1Bits[i]
					bits2[i] = pulse_overlayed_master2Bits[i]
					g_difference[i] = g_overlayed_master1[i] - g_overlayed_master2[i]
				
				
				#plt.figure(1)
				#plt.plot(x, difference, plot_style[1%len(plot_style)])
				#plt.figure(2)
				#plt.plot(x, ratio, plot_style[1%len(plot_style)])
				#plt.plot(np.linspace(150, 179, 30), ratio[150:180], plot_style[1%len(plot_style)])
				
				
				
				
				ax = plt.subplot(311) #Top Plot: ADC Samples
				plt.title("ADC Samples")
				#xmin, xmax, ymin, ymax = plt.axis()
				#plt.ylim(ymin, ymax*1.1)
				# plt.axis((xmin, 960, ymin, ymax*1.1))
				ax.set_xticks(np.arange(0, 1201, 120))
				ax.set_xticks(np.arange(0, 1201, 30), minor=True)
				plt.grid(which='minor', alpha=0.4)
				plt.grid(which='major', alpha=0.7)
				plt.plot(x, pulse_overlayed_master1, '.-')
				plt.plot(x, pulse_overlayed_master2, '.-')
				xmin, xmax, ymin, ymax = plt.axis()
				plt.ylim(ymin, ymax*1.1)
				ax = plt.subplot(312) #Top Plot: ADC Samples
				plt.title("Difference of ADC Samples (Series 1 - Series 2)")
				ax.set_xticks(np.arange(0, 1201, 120))
				ax.set_xticks(np.arange(0, 1201, 30), minor=True)
				plt.grid(which='minor', alpha=0.4)
				plt.grid(which='major', alpha=0.7)
				plt.plot(x, difference, '.-')
				ax = plt.subplot(313) #Bottom Plot: Decision Bits
				plt.title("Decision Bits")
				plt.ylim(-0.1, 1.1)
				#plt.plot(x, bits1, '.-')
				plt.plot(x, bits2, '.-')
				plt.subplots_adjust(hspace=0.4)
				#plt.xlabel("Comparing: "+filenames[0].replace("/Averaged/", ", ").replace("_Averaged.npy", " (Averaged)").replace("Sine_","").replace("0p","0.") + " with " + filenames[1].replace("/Averaged/", ", ").replace("_Averaged.npy", " (Averaged)").replace("Sine_","").replace("0p","0."))
				plt.xlabel("Gain: "+str(np.round(gc, 4))+", Series 1: Amp"+amps[amp]+", "+gains[gain1]+" --- Series 2: Amp"+amps[amp]+", "+gains[gain2])
				ax.set_xticks(np.arange(0, 1201, 120))
				ax.set_xticks(np.arange(0, 1201, 30), minor=True)
				plt.grid(which='minor', alpha=0.4)
				plt.grid(which='major', alpha=0.7)
				#plt.show()
				plt.savefig(""+checkDir(plotsDir)+"/Pulse_Amp"+amps[amp]+"_"+gains[gain1]+"_"+gains[gain2]+".png")
				print("Saved File: "+checkDir(plotsDir)+"/Pulse_Amp"+amps[amp]+"_"+gains[gain1]+"_"+gains[gain2]+".png")
				plt.clf()
				
				ax = plt.subplot(311) #Top Plot: ADC Samples
				plt.title("Derivative: ADC Samples")
				#xmin, xmax, ymin, ymax = plt.axis()
				#plt.ylim(ymin, ymax*1.1)
				# plt.axis((xmin, 960, ymin, ymax*1.1))
				ax.set_xticks(np.arange(0, 1201, 120))
				ax.set_xticks(np.arange(0, 1201, 30), minor=True)
				plt.grid(which='minor', alpha=0.4)
				plt.grid(which='major', alpha=0.7)
				plt.plot(x, g_overlayed_master1, '.-')
				plt.plot(x, g_overlayed_master2, '.-')
				xmin, xmax, ymin, ymax = plt.axis()
				plt.ylim(ymin, ymax*1.1)
				ax = plt.subplot(312) #Top Plot: ADC Samples
				plt.title("Derivative: Difference of ADC Samples (Series 1 - Series 2)")
				ax.set_xticks(np.arange(0, 1201, 120))
				ax.set_xticks(np.arange(0, 1201, 30), minor=True)
				plt.grid(which='minor', alpha=0.4)
				plt.grid(which='major', alpha=0.7)
				plt.plot(x, g_difference, '.-')
				ax = plt.subplot(313) #Bottom Plot: Decision Bits
				plt.title("Decision Bits")
				plt.ylim(-0.1, 1.1)
				#plt.plot(x, bits1, '.-')
				plt.plot(x, bits2, '.-')
				plt.subplots_adjust(hspace=0.4)
				#plt.xlabel("Comparing: "+filenames[0].replace("/Averaged/", ", ").replace("_Averaged.npy", " (Averaged)").replace("Sine_","").replace("0p","0.") + " with " + filenames[1].replace("/Averaged/", ", ").replace("_Averaged.npy", " (Averaged)").replace("Sine_","").replace("0p","0."))
				plt.xlabel("Gain: "+str(np.round(gc, 4))+", Series 1: Amp"+amps[amp]+", "+gains[gain1]+" --- Series 2: Amp"+amps[amp]+", "+gains[gain2])
				ax.set_xticks(np.arange(0, 1201, 120))
				ax.set_xticks(np.arange(0, 1201, 30), minor=True)
				plt.grid(which='minor', alpha=0.4)
				plt.grid(which='major', alpha=0.7)
				#plt.show()
				plt.savefig(""+checkDir(plotsDir)+"/Derivative_Amp"+amps[amp]+"_"+gains[gain1]+"_"+gains[gain2]+".png")
				print("Saved File: "+checkDir(plotsDir)+"/Derivative_Amp"+amps[amp]+"_"+gains[gain1]+"_"+gains[gain2]+".png")
				plt.clf()

if __name__ == "__main__":

	for version in range(len(versions)):
		amps = getAmps(dataDir+versions[version])
		compareTwo(versions[version])
