import numpy as np
import matplotlib.pyplot as plt
import sys

if __name__ == "__main__":


	filename = sys.argv[1]
	
	data = np.load(filename)
	
	ax = plt.subplot(211) #Top Plot: ADC Samples
	plt.title("Modified ADC Samples")
	plt.plot(np.linspace(0, len(data[0])-1, len(data[0])), data[0], '.-')
	ax = plt.subplot(212) #Bottom Plot: Decision Bits
	plt.title("Decision bits")
	plt.xlabel("Filename: " + filename)
	plt.plot(np.linspace(0, len(data[1])-1, len(data[1])), data[1], '.-')
	plt.ylim(-0.2, 1.2)
	plt.show()
	plt.clf()
	

