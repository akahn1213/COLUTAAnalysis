#/usr/local/bin/python3

import numpy as np

import os
import glob
import sys
import h5py

#Versions: The list of directory names in Data/Processed/Pulses that the analysis will run on
#versions = ["10MSPS", "40MSPS"]
versions = ["Feb6"]

#getVersions: return the list of versions
def getVersions():
  return versions

#getAmps: Find the list of amplitudes, to be used to find and make the lists of files. i.e. ["0p1", "0p25", "0p9", ....] etc
def getAmps(dir): 
  amps = []
  for filename in os.listdir(dir):
    amps.append((filename.split("Amp")[1]).split("_")[0])
  return np.unique(amps)


def get_rms(data):
  data = data.astype(float)
  return np.sqrt( np.sum( np.square(data) )/len(data))



#getAmpNums replacing the p in the list of amps with a period
def getAmpNums(amps):
  ampNums = np.empty(len(amps))
  for i in range(len(amps)):
    #ampNum = float(str(amps[i].decode("utf-8")).replace("p", "."))
    ampNums[i] = np.asarray(float(amps[i].replace("p", "."))).astype(np.float)
    #ampNums[i] = ampNum
  return ampNums

#checkDir: Make sure a directory exists. If not, make it. Then, return the directory name so that it can be used in save commands i.e. plt.savefig(checkDir("Plots/subdirectory1/"))
def checkDir(directory):
  if not os.path.exists(os.path.dirname(directory)):
    try:
      os.makedirs(os.path.dirname(directory))
    except OSError as exc: # Guard against race condition
      if exc.errno != errno.EEXIST:
        raise
  return directory

def makeH5File(name, dr = None):
  if os.path.exists(name):
    os.remove(name)
  if dr is not None:
    return h5py.File(name, 'w', driver=dr)
  else:
    return h5py.File(name, 'w')

#getNFiles: Get the number of files that is used for each amplitude. Uses the names of the processed data files in Data/Processed/Pulses/version to find how many files correspond to the same amplitude and gain
def getNFiles(filename):
  return len(glob.glob(filename+"*"))

def getRootDir():
  return os.path.abspath(os.path.join(sys.path[0], os.pardir))


"""
CalcSumSquares: Returns the sum of differences squared between elements of two lists (must be of the same size)
"""
def calcSumSquares(list1, list2):

    sumSquares = 0
    for i in range(len(list1)):
      sumSquares += np.square(list1[i] - list2[i])  
    return sumSquares


def short(directory):
  return os.path.basename(os.path.normpath(directory)) + "/"

def getMatrixAvg(m, xi, xf, yi, yf):
  return np.mean(m[np.ix_(np.arange(yi, yf), np.arange(xi, xf))].mean(0)) 



def getDivisionUncertainty(a, b, sig_a, sig_b):
  t1 = np.square(np.divide(sig_a,a))
  t2 = np.square(np.divide(sig_b,b))
  t3 = np.multiply(2, np.divide(np.multiply(sig_a,sig_b)  ,  np.multiply(a, b)))
  return np.multiply(np.divide(a, b), np.sqrt(np.add(t1, np.subtract(t2, t3))))
