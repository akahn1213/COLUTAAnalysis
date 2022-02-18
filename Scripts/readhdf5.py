import h5py
import numpy as np
import os
import time
import matplotlib.pyplot as plt
import helperfunctions as hf


start = time.time()

raw_dir = hf.getRootDir()+"/Data/Raw/Pulses/"
filename = raw_dir+"Run_1232_Output.hdf5"

f = h5py.File(filename, "r")
print("\n\n\n##########################")
print("###### HDF5 Example ######")
print("##########################\n\n\n")


print("------------------------")
print("------ Attributes ------")
print("------------------------\n\n\n")


for key, val in f.attrs.items():
  print ("%s: %s\n" % (key, val))


print(f["Measurement_400"].attrs["pulse_amplitude"])


print("\n\n-----------------------------------------")
print("------ Group and Dataset Structure ------")
print("-----------------------------------------\n\n\n")
"""
for key in f.keys():
  print("\n"+f[key].name)
  print("   -Attributes-")
  for k, v, in f[key].attrs.items():
    print("   - %s: %s" % (k, v))
  for key2 in f[key].keys():
    n = f[key+"/"+key2].name
    print("\n"+n)
    print("   -Attributes-")
    for k, v, in f[key+"/"+key2].attrs.items():
      print("   - %s: %s" % (k, v))

    #print(n+": dataset shape = "+str(np.shape(f[n][()])))
print("\n\n")
"""



"""
if(os.path.getsize(filename)/1048576 > 1.):
  print("File Size: "+str(np.round(os.path.getsize(filename)/1048576, 2))+"MB\n")
else:
  print("File Size: "+str(np.round(os.path.getsize(filename)/1024, 2))+"KB\n")

print("Loading and Processing Time: "+str(np.round((time.time() - start)*1000., 3))+"ms\n\n") 
"""



samples = f["Measurement_205/coluta1/channel2/samples"][()]
bits = f["Measurement_205/coluta1/channel2/bits"][()]
#bits = f["adc0/bits"][()]
ax = plt.subplot(211) #Top Plot: ADC Samples
plt.title("Modified ADC Samples")
plt.plot(np.linspace(0, len(samples)-1, len(samples)), samples, '.-')
ax = plt.subplot(212) #Bottom Plot: Decision Bits
plt.title("Decision bits")
plt.xlabel("Filename: " + filename)
plt.plot(np.linspace(0, len(bits)-1, len(bits)), bits, '.-')
plt.ylim(-0.2, 1.2)
plt.show()
plt.clf()


