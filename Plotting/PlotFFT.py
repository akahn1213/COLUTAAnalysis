##################################################################
#
# arguments = 1x or AG (for now)
#
##################################################################
import numpy as np
import matplotlib.pyplot as plt
import sys
from scipy.fftpack import fft
from matplotlib.colors import LogNorm

if __name__ == "__main__":

	datapath = "../Data/Processed/data"+sys.argv[1]+"_processed.npy"
	plotpath = "PlotsSaved/fft_"+sys.argv[1]+".png"

	print('Reading in ' + datapath )
	data = np.load( datapath )
	gain_bit = data[0]
	adc_counts = data[1]

	print('Calculating FFT...')
	fft = fft( adc_counts[500:1000] )

        #Plot FFT
        N = len(fft)
	print(N)
        T = 25 #ns
        x = np.linspace(0.0, N*T, N)
        xf = np.linspace(0.0, 1.0/(2.0*T), N/2)
        x = np.arange(len(fft)/2.)

	print('Plotting FFT...')
	plt.plot(xf, 2.0/N * np.abs(fft[0:N/2]), "b.-")
	plt.semilogy()
	plt.title('FFT') 
	plt.savefig( plotpath )
	plt.clf()

