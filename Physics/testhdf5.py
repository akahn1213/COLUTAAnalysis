import h5py
import numpy as np
from math import *
import matplotlib.pyplot as plt
import argparse
import glob
from scipy import stats
from helperfunctions import *
import pickle


dataDir = checkDir(getRootDir()+"/Data/Processed/Pulses/")
#Get directories to run in
parser = argparse.ArgumentParser()
parser.add_argument("-r", "--runs", default = [], type=str, nargs='+',
                   help="list of runs (directories in Physics/Blocks_Final/) to include")
#Allow user input for choice of runs
args = parser.parse_args()
#runs = args.runs
runs = [dataDir + run + "/" for run in args.runs]
directories = []

#Search for runs if none are provided
if runs == []:
  directories = glob.glob(dataDir+"*/")
  if directories == []:
    print("No runs found in "+dataDir+" exiting...")
    sys.exit(0)
  for directory in directories:
    if glob.glob(directory+"*.hdf5"):
      runs.append(directory)
  print("Found blocks in the following directories:\n")
  for run in runs:
    print("-%s" % short(run))
  continue_script = input("\nContinue running? [y/n] (default=y): ")
  if(continue_script == "n"):
    sys.exit(0)

#amps = getAmps(runs[0])

amps = ["1p0", "0p75", "0p5", "0p25", "0p1"]
dat = [ [] for i in range(len(amps)*3)]
#dat- list for saving all amplitudes
#amps = ["0p75"]
j = 0 #for indexing the data array
for run in runs:
  saveDir = "Physics/Shape_Plots/"+short(run)
  for amp in amps:
    f = h5py.File("Blocks_Final/"+short(run)+"Blocks_Amp"+amp+"_1x.hdf5", "r")
  
    #Sample Pulse
    plt.plot(np.gradient(f["/coluta1/channel1/samples"][()][0][0:1200]), '.-')
    plt.xticks(np.arange(0,1300,100)); plt.xlabel('ns', horizontalalignment='right', x=1.0)
    plt.title(short(run)+" COLUTA 1, Channel 1, Example Pulse Derivative, Amp "+amp)
    plt.savefig(saveDir+"derivative_Co1Ch1_Amp"+amp)
    plt.clf()
    #plt.show()

    plt.plot(np.gradient(f["/coluta1/channel2/samples"][()][0][0:1200]), '.-')
    plt.xticks(np.arange(0,1300,100)); plt.xlabel('ns', horizontalalignment='right', x=1.0)
    plt.title(short(run)+" COLUTA 1, Channel 2, Example Pulse Derivative, Amp "+amp)
    plt.savefig(saveDir+"derivative_Co1Ch2_Amp"+amp)
    plt.clf()

    plt.plot(np.gradient(f["/coluta2/channel1/samples"][()][0][0:1200]), '.-')
    plt.xticks(np.arange(0,1300,100)); plt.xlabel('ns', horizontalalignment='right', x=1.0)
    plt.title(short(run)+" COLUTA 2, Channel 1, Example Pulse Derivative, Amp "+amp)
    plt.savefig(saveDir+"derivative_Co2Ch1_Amp"+amp)
    plt.clf()


    plt.plot(f["/coluta1/channel1/samples"][()][0][0:1200], '.-')
    plt.xticks(np.arange(0,1300,100)); plt.xlabel('ns', horizontalalignment='right', x=1.0)
    plt.title(short(run)+" COLUTA 1, Channel 1, Example Pulse, Amp "+amp)
    plt.savefig(saveDir+"pulse_Co1Ch1_Amp"+amp)
    #plt.show()
    plt.clf()

    plt.plot(f["/coluta1/channel2/samples"][()][0][0:1200], '.-')
    plt.xticks(np.arange(0,1300,100)); plt.xlabel('ns', horizontalalignment='right', x=1.0)
    plt.title(short(run)+" COLUTA 1, Channel 2, Example Pulse, Amp "+amp)
    plt.savefig(saveDir+"pulse_Co1Ch2_Amp"+amp)
    #plt.show()
    plt.clf()

    plt.plot(f["/coluta2/channel1/samples"][()][0][0:1200], '.-')
    plt.xticks(np.arange(0,1300,100)); plt.xlabel('ns', horizontalalignment='right', x=1.0)
    plt.title(short(run)+" COLUTA 2, Channel 1, Example Pulse, Amp "+amp)
    plt.savefig(saveDir+"pulse_Co2Ch1_Amp"+amp)
    #plt.show()
    plt.clf()
  
    blocks11 = f["/coluta1/channel1/samples"][()]
    blocks21 = f["/coluta2/channel1/samples"][()]
    blocks12 = f["/coluta1/channel2/samples"][()]
    #plt.plot(np.mean(blocks, axis=0))
    #plt.title(short(run)+"COLUTA 1, Channel 1, Amp " + str(amp))
    #plt.show()
    #plt.savefig("Run1238_COLUTA 1_Channel 1_Amp" + str(amp))
  
    #print np.argmax(np.mean(blocks, axis=0))
  

    #Amplitude pulse shape depenedence tests
    norm_fact11 = max(np.mean(blocks11, axis=0))
    norm_fact21 = max(np.mean(blocks21, axis=0))
    norm_fact12 = max(np.mean(blocks12, axis=0))
    #print(str(norm_fact) + "Max of pulse for amp " + amp)
    #plt.plot((np.mean(blocks, axis=0))/norm_fact)
    #plt.title(short(run)+"COLUTA 1, Channel 1, Amp Compare")

    dat[j] = np.mean(blocks11, axis=0)/norm_fact11
    dat[j+len(amps)] = np.mean(blocks21, axis=0)/norm_fact21
    dat[j+len(amps)*2] = np.mean(blocks12, axis=0)/norm_fact12
    j+=1
    #Histograms along pulse
    points_dict = {
      "20":"Baseline Before Pulse",
      "142":"Rising Edge",
      "165":"Peak",
      "195":"Falling Edge",
      "500":"Negative Lobe",
      "1000":"Baseline After Pulse"
      }


  
    for channel in range(2): 
      n_samp = len(f["/coluta1/channel1/samples"][()][0])
      t_data = np.transpose(f["/coluta1/channel"+str(channel+1)+"/samples"][()])
      rms_data = np.empty(n_samp)
      #for i in [20, 142, 165, 195, 500, 1000]: 
      for point, desc in points_dict.items(): 
  
        data = t_data[int(point)]
        ax=plt.subplot(111)
        plt.hist(data, bins = range(int(min(data)), int(max(data))+1), normed=False)
        xt = plt.xticks()[0]  
        xmin, xmax = min(xt), max(xt)  
        lnspc = np.linspace(xmin, xmax, len(data))
        m, s = stats.norm.fit(data) # get mean and standard deviation  
        #pdf_g = stats.norm.pdf(lnspc, m, s) # now get theoretical values in our interval  
        #plt.plot(lnspc, pdf_g, label="Gaussian Fit") # plot it
        plt.title("Sample "+str(point)+": "+desc+", Channel "+str(channel+1))
        plt.text(0.95, 0.85, r'$\mu = $'+str(np.round(np.mean(data), 3)), horizontalalignment='right', verticalalignment='center', transform=ax.transAxes) 
        plt.text(0.95, 0.75, r'$\sigma = $'+str(np.round(np.std(data), 3)), horizontalalignment='right', verticalalignment='center', transform=ax.transAxes) 
        #plt.show()
        #plt.savefig("dist_s"+str(point)+"_c"+str(channel+1)+".png")
        plt.clf()
  
  
  
  
  
      #RMS
      for i in range(n_samp): 
        #rms_data[i] = get_rms([x - np.mean(t_data[i]) for x in t_data[i]])
        rms_data[i] = get_rms(t_data[i] - np.mean(t_data[i]))
        #rms_data[i] = np.std(t_data[i])
   
  
  
  
      ax = plt.subplot(111)
      plt.plot(rms_data[0:1200], '.-', label="Channel "+str(channel+1))
      plt.title('RMS, Amplitude '+amp.replace("p", "."))
      plt.savefig("rmsedit"+str(channel+1)+".png")
      plt.clf()
  
  #Plot normalized pulses and deviations
  colorscheme = ['r', 'b', 'g', 'm', 'y']
  fig,ax = plt.subplots(3)
  fig1,ax1 = plt.subplots(3)
  for i in range(3):
    ax[i].plot(dat[0+len(amps)*i][0:(floor(len(dat[0])/2))], 'r', linewidth =1, alpha = 0.7, label='1p0')
    ax[i].plot(dat[1+len(amps)*i][0:(floor(len(dat[0])/2))], 'b', linewidth =1, alpha = 0.7, label='0p75')
    ax[i].plot(dat[2+len(amps)*i][0:(floor(len(dat[0])/2))], 'g', linewidth =1, alpha = 0.7, label='0p5')
    ax[i].plot(dat[3+len(amps)*i][0:(floor(len(dat[0])/2))], 'm', linewidth =1, alpha = 0.7, label='0p25')
    ax[i].plot(dat[4+len(amps)*i][0:(floor(len(dat[0])/2))], 'y', linewidth =1, alpha = 0.7, label='0p1')
    if i==0:
      ax[i].legend(loc='upper right')
      ax[i].set_title("Normalized Pulses - Run1260/ COLUTA1 Channel1")
      #pickle.dump(fig,open(saveDir+'normalizedPulses_Co1Ch1.fig.pickle', 'wb'))
    elif i==1:
      ax[i].set_title("Normalized Pulses - Run1260/ COLUTA2 Channel1")
      #pickle.dump(fig,open(saveDir+'normalizedPulses_Co2Ch1.fig.pickle', 'wb'))
    else:
      ax[i].set_title("Normalized Pulses - Run1260/ COLUTA1 Channel2")
      #pickle.dump(fig,open(saveDir+'normalizedPulses_Co1Ch2.fig.pickle', 'wb'))
    #plt.show()
    #plt.clf()
    ax1[i].plot(dat[1+len(amps)*i][0:(floor(len(dat[0])/2))]-dat[0+len(amps)*i][0:(floor(len(dat[0])/2))], 'b', alpha = 0.7, label='0p75')
    ax1[i].plot(dat[2+len(amps)*i][0:(floor(len(dat[0])/2))]-dat[0+len(amps)*i][0:(floor(len(dat[0])/2))], 'g', alpha = 0.7, label='0p5')
    ax1[i].plot(dat[3+len(amps)*i][0:(floor(len(dat[0])/2))]-dat[0+len(amps)*i][0:(floor(len(dat[0])/2))], 'm', alpha = 0.7, label='0p25')
    ax1[i].plot(dat[4+len(amps)*i][0:(floor(len(dat[0])/2))]-dat[0+len(amps)*i][0:(floor(len(dat[0])/2))], 'y', alpha = 0.7, label='0p1')
    if i==0:
      ax1[i].legend(loc='upper right')
      ax1[i].set_title("Deviation from Amp=1.0 - Run1260/ COLUTA1 Channel1")
      #pickle.dump(fig1,open(saveDir+'deviation_Co1Ch1.fig.pickle', 'wb'))
    elif i==1:
      ax1[i].set_title("Deviation from Amp=1.0 - Run1260/ COLUTA2 Channel1")
      #pickle.dump(fig1,open(saveDir+'deviation_Co2Ch1.fig.pickle', 'wb'))
    else:
      ax1[i].set_title("Deviation from Amp=1.0 - Run1260/ COLUTA1 Channel2")
      #pickle.dump(fig1,open(saveDir+'deviation_Co1Ch2.fig.pickle', 'wb'))
    #plt.show()
    #plt.clf()
  for a in fig.get_axes():
      a.set_xlabel('ns')
      a.label_outer(); a.xaxis.set_major_locator(plt.MultipleLocator(100))
  for a1 in fig1.get_axes():
      a1.set_xlabel('ns')
      a1.label_outer()
      a1.xaxis.set_major_locator(plt.MultipleLocator(100))


  pickle.dump(fig,open(saveDir+'normalizedPulses.fig.pickle', 'wb'))
  pickle.dump(fig1,open(saveDir+'deviationPulses.fig.pickle', 'wb'))
  plt.close(fig)
  plt.close(fig1)

  #High vs low gain -one plot
  plt.plot(dat[0][0:900]-dat[0+len(amps)*2][0:900],'r', alpha = 0.7, label='1p0')
  plt.plot(dat[1][0:900]-dat[1+len(amps)*2][0:900],'b', alpha = 0.7, label='0p75')
  plt.plot(dat[2][0:900]-dat[2+len(amps)*2][0:900],'g', alpha = 0.7, label='0p5')
  plt.plot(dat[3][0:900]-dat[3+len(amps)*2][0:900],'m', alpha = 0.7, label='0p25')
  plt.plot(dat[4][0:900]-dat[4+len(amps)*2][0:900],'y', alpha = 0.7, label='0p1')
  plt.legend()
  plt.title("Diff between high and low gain (CH1 -CH2)")
  plt.savefig(saveDir+"highVSlowGain_allamps.png")
  #plt.show()
  plt.clf()

  #high vs low gain - pulse shapes
  fig2,ax2 = plt.subplots(2,3)
  j=0; k=0; flag = False;
  for i in range(len(amps)):
    if i==2 or i==4:
      k+=1
    if flag:
      j=1
      flag = not(flag)
    elif not(flag):
      j=0
      flag = not(flag)
    ax2[j,k].plot(dat[i][0:900], colorscheme[i], label='low')
    ax2[j,k].plot(dat[i+len(amps)*2][0:900], colorscheme[i]+':', label='high')
    if i==0:
      ax2[j,k].legend()
    ax2[j,k].set_title('Amp' + amps[i])

  for a2 in fig.get_axes():
    a2.label_outer()
    a2.label_outer(); a2.xaxis.set_major_locator(plt.MultipleLocator(100))
  pickle.dump(fig2,open(saveDir+'highGainVSlowGain.fig.pickle', 'wb'))
  plt.close(fig2)

  for i in range(len(amps)):
      fig,axs = plt.subplots(2)
      axs[0].plot(dat[i][0:900], colorscheme[i], label='low')
      axs[0].plot(dat[i+len(amps)*2][0:900], colorscheme[i]+':', label='high')
      axs[1].plot(dat[i][0:900]-dat[i+len(amps)*2][0:900], 'black')
      axs[0].legend()
      axs[0].set_title('Amp' + amps[i])
      axs[1].set_title('Vertical Deviation')
      for a in fig.get_axes():
          a.set_xlabel('ns')
          a.label_outer(); a.xaxis.set_major_locator(plt.MultipleLocator(100))
      pickle.dump(fig,open(saveDir+'highGainVSlowGain_Amp'+amps[i]+'.fig.pickle', 'wb'))
      plt.close(fig)

  #find avg of amplitdes
  avg = np.mean(dat, axis=0)
  #print(avg)

  #plot deviations
  plt.plot(dat[0][:]-avg, 'r', alpha = 0.7, label='1p0')
  plt.plot(dat[1][:]-avg, 'b', alpha = 0.7, label='0p75')
  plt.plot(dat[2][:]-avg, 'g', alpha = 0.7, label='0p5')
  plt.plot(dat[3][:]-avg, 'm', alpha = 0.7, label='0p25')
  plt.plot(dat[4][:]-avg, 'y', alpha = 0.7, label='0p1')
  plt.legend()
  plt.title("Deviation from Average")
  #plt.show()
  plt.clf()
  
  #plot deivations

  
    #plt.title('RMS - Mean, Amplitude '+amp.replace("p", "."))
    #plt.title('RMS, Amplitude '+amp.replace("p", "."))
    #plt.legend()
    #plt.savefig("rmsedit.png")
    #plt.show()
    #for key, val in f.attrs.items():
    #  print ("%s: %s\n" % (key, val))
    
    
  """  
    for key in f.keys():
      print("\n"+f[key].name)
      print("   -Attributes-")
      for k, v, in f[key].attrs.items():
        print("   - %s: shape = %s" % (k, np.shape(v)))
      for key2 in f[key].keys():
        n = f[key+"/"+key2].name
        print(n+": dataset shape = "+str(np.shape(f[n][()])))
  """
  """
    for i in range(0,10):
      #Check interlacing
      s = f["/adc0/samples"][()][0][i]
      b = f["/adc0/bits"][()][0][i]
  
  
      #print np.std(s[0:20])
      
      ax = plt.subplot(111)
      plt.plot(np.linspace(0, len(s)-1, len(s)), s, '.-')
      #plt.title(r'Amp'+amp+': Pulse Baseline: $\sigma$ of first 40 samples = '+str(np.round(np.std(s[0:40]), 3)))
      plt.title('Run 1230, COLUTA 1, Channel 2, Amplitude '+amp.replace("p", "."))
      plt.text(0.55, 0.8, r'$\sigma$ of first 40 samples = '+str(np.round(np.std(s[0:40]), 3)), transform=ax.transAxes)
  #    plt.savefig(getRootDir()+"/Plotting/PlotsSaved/Pulses/OFCs/Run_1147/Pulse_Example_Channel1_Amp"+amp+".png")
  #    plt.clf()
      plt.show()
  """
  #  plt.plot(np.linspace(0, len(b)-1, len(b)), b, '.-')
  #  plt.show()
    #Check alignment
  """
    s = []*10
    for j in range(20):
      for i in range(5):
        s = np.asarray(f["/adc0/samples"][()][1])
        plt.plot(np.linspace(145, 164, 20), s[i+5*j][145:165], '.-')
        #plt.plot(np.linspace(0, 1199, 1200), s[i*5 + j], '.-')
      plt.title("Amp"+amp)
      plt.show()
      plt.clf()
  """
