import h5py
import numpy as np
import matplotlib.pyplot as plt



filter_out = h5py.File('Filter_Output/Run_1147/OF/OF_Output.hdf5', 'r') 
#energy = filter_out["/adc"+str(adc)+"/channel"+str(channel)+"/"+str(n_samples)+"S/"+gain+"/Amp"+amps[amp]+"/energy"][()]
energy = filter_out["/adc0/channel0/5S/1x/Amp1p0/energy"][()]
print(energy)
#plt.plot(energy)
#plt.show()
print(np.mean(energy))
print(np.std(energy))
