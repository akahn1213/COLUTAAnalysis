import h5py
import numpy as np
from math import *
import matplotlib.pyplot as plt
import argparse
import glob
from scipy import stats
from helperfunctions import *
import pickle
from itertools import product

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

# amps = getAmps(runs[0])
amps = ["0p2"]

#amps = ["1p0", "0p75", "0p5", "0p25", "0p1"]
for run in runs:
  saveDir = checkDir(getRootDir()+"/Physics/Shape_Plots/"+short(run))

  # Point for histograms along pulse
  points_dict = {
    "20":"Baseline Before Pulse",
    "142":"Rising Edge",
    "165":"Peak",
    "195":"Falling Edge",
    "500":"Negative Lobe",
    "1000":"Baseline After Pulse"
  }

  for amp in amps:
    f = h5py.File(getRootDir()+"/Physics/Blocks_Final/"+short(run)+"Blocks_Amp"+amp+"_1x.hdf5", "r")

    #Identify variable attributes
    default_attrs = ['adc_freq', 'awg_freq', 'block_length', 'n_adcs', 'n_channels', 'attribute_combos']
    cut_attrs = [attr for attr in f.attrs if attr not in default_attrs]
    #print(cut_attrs)
    cut_combos = [combo.split('-') for combo in f.attrs["attribute_combos"]]
    print(cut_combos)
    # print(cut_combos)

    #Dictionaries to save data
    dat11 = {}; dat12 = {};
    norm_dat11 = {}; norm_dat12 = {};
    shape_time11 = {}; shape_time12 = {};
    height11 = {}; height12 = {}

    #Loop over all data combinations
    for prod in cut_combos:
      #print(prod)
      #Select and name data
      # https://stackoverflow.com/a/41922008
      #print(np.array( [f.attrs[cut_attrs[i]] == prod[i] for i in range(len(prod))] ))
      cut_index = np.logical_and.reduce(np.array( [f.attrs[cut_attrs[i]] == prod[i] for i in range(len(prod))] ))
      if not any(cut_index): continue # don't make plots if we don't have that combo
      cut_attrs_str_name = "_".join([cut_attrs[i]+'_'+str(prod[i]) for i in range(len(prod))])
      cut_attrs_str_title = "; ".join([cut_attrs[i]+' '+str(prod[i]) for i in range(len(prod))])

      for channel in range(1,3):

        pulse_data = f["/coluta1/channel"+str(channel)+"/samples"][()][cut_index][0:1200]

        #Plot Average Pulse Derivative and Average Pulse
        plt.plot(np.gradient(np.mean(pulse_data,axis=0)), '.-')
        plt.xticks(np.arange(0,1300,100))
        plt.xlabel('ns', horizontalalignment='right', x=1.0)
        plt.title(short(run)+" COLUTA 1, Channel "+str(channel)+", Average Pulse Derivative, Amp "+amp+"\n"+cut_attrs_str_title)
        plt.savefig(saveDir+"derivative_Co1Ch"+str(channel)+"_Amp"+amp+'_'+cut_attrs_str_name)
        # plt.show()
        plt.clf()

        plt.plot(np.mean(pulse_data,axis=0), '.-')
        plt.xticks(np.arange(0,1300,100))
        plt.xlabel('ns', horizontalalignment='right', x=1.0)
        plt.title(short(run)+" COLUTA 1, Channel "+str(channel)+", Example Pulse, Amp "+amp+"\n "+cut_attrs_str_title)
        plt.savefig(saveDir+"pulse_Co1Ch"+str(channel)+"_Amp"+amp+"_"+cut_attrs_str_name)
        #plt.show()
        plt.clf()

        for point, desc in points_dict.items():

          data = pulse_data[:,int(point)]
          ax=plt.subplot(111)
          plt.hist(data, bins = range(int(min(data)), int(max(data))+1), normed=False)
          m, s = stats.norm.fit(data) # get mean and standard deviation
          plt.title("Sample "+str(point)+": "+desc+", Channel "+str(channel)+", Amp "+amp+"\n"+cut_attrs_str_title)
          plt.text(0.95, 0.85, r'$\mu = $'+str(np.round(np.mean(data), 3)), horizontalalignment='right', verticalalignment='center', transform=ax.transAxes)
          plt.text(0.95, 0.75, r'$\sigma = $'+str(np.round(np.std(data), 3)), horizontalalignment='right', verticalalignment='center', transform=ax.transAxes)
          #plt.show()
          plt.savefig(checkDir(saveDir+"/Noise_Hists/")+"dist_s"+str(point)+"_c"+str(channel)+"_Amp"+amp+"_"+cut_attrs_str_name+".png")
          plt.clf()


      #Data to save
      blocks11 = np.mean(f["/coluta1/channel1/samples"][()][cut_index], axis = 0)
      blocks12 = np.mean(f["/coluta1/channel2/samples"][()][cut_index], axis = 0)

      #Normalization factors
      norm_fact11 = max(blocks11)
      norm_fact12 = max(blocks12)

      #5% - 100% shape time calculution
      time100_11 = np.argmax(blocks11); time100_12 = np.argmax(blocks12)
      baseline_11 = np.mean(blocks11[:5]); baseline_12 = np.mean(blocks12[:5])
      max11 = max(blocks11); max12 = max(blocks12);
      blocks11_shifted = blocks11 - baseline_11
      blocks12_shifted = blocks12 - baseline_12
      time05_11 = 9999; time05_12 = 9999;
      for val1 in blocks11_shifted:
        if val1 >= 0.05*max(blocks11_shifted):
          time05_11 = np.where(blocks11_shifted == val1)[0][0]
          break
      for val2 in blocks12_shifted:
        if val2 >= 0.05*max(blocks12_shifted):
          time05_12 = np.where(blocks12_shifted==val2)[0][0]
          break
      if time05_11 == 9999:
        print('Uh oh. No value above 5 percent found for COLUTA1 Channel 1')
      if time05_12 == 9999:
        print('Uh oh. No value above 5 percent found for COLUTA1 Channel 2')
      shape_time11.update({cut_attrs_str_title:time100_11-time05_11})
      shape_time12.update({cut_attrs_str_title:time100_12-time05_12})
      height11.update({cut_attrs_str_title:max11 - baseline_11})
      height12.update({cut_attrs_str_title:max12 - baseline_12})

      print("0.05 - 1.0 Shaping Time for COLUTA1 Channel 1, with attributes " + cut_attrs_str_title + '  ->  ' + str(shape_time11[cut_attrs_str_title]) + ' ns')
      print("0.05 - 1.0 Shaping Time for COLUTA1 Channel 2, with attributes " + cut_attrs_str_title + '  ->  ' + str(shape_time12[cut_attrs_str_title]) + ' ns')

      #Save the data
      norm_dat11.update({cut_attrs_str_title:blocks11/norm_fact11})
      dat11.update({cut_attrs_str_title:blocks11})
      norm_dat12.update({cut_attrs_str_title:blocks12/norm_fact12})
      dat12.update({cut_attrs_str_title:blocks12})

      points_dict = {
        "20":"Baseline Before Pulse",
        "142":"Rising Edge",
        "165":"Peak",
        "195":"Falling Edge",
        "500":"Negative Lobe",
        "1000":"Baseline After Pulse"
        }
    if amp == '0p2':
      for channel in range(2):
        if channel == 0:
          dat = dat11
          st = shape_time11
          hgt = height11
        else:
          dat = dat12
          st = shape_time12
          hgt = height12
        #print(dat11.keys())
        on_g20_true = 'dac_ibi_g20 63; hg_lg_c2 220; on_g20 True; shaper_constants 14_6_15; sw_ibo_g20 False'
        on_g20_false = 'dac_ibi_g20 63; hg_lg_c2 220; on_g20 False; shaper_constants 14_6_15; sw_ibo_g20 False'
        l1, = plt.plot(dat[on_g20_true])
        l2, = plt.plot(dat[on_g20_false])
        plt.title('Channel' + str(channel+1) +': on_g20 Toggle')
        plt.xlabel('ns')
        legend1 = plt.legend([l1,l2], ["on_g20 True (default)", 'on_g20 False'], loc = 'upper right')
        legend2 = plt.legend([l1,l2], ["Rise time: " + str(st[on_g20_true]) + " ns", 'Rise time: ' + str(st[on_g20_false]) + " ns"], loc = 'upper right', bbox_to_anchor= (1.0,0.7))
        legend3 = plt.legend([l1,l2], ["Height: " + str(int(hgt[on_g20_true])) + " counts", 'Height: ' + str(int(hgt[on_g20_false])) + " counts"], loc = 'upper right', bbox_to_anchor= (1.0,0.85))
        ax1 = plt.gca().add_artist(legend1)
        ax2 = plt.gca().add_artist(legend2)
        ax3 = plt.gca().add_artist(legend3)
        plt.savefig(saveDir + 'channel' + str(channel+1) + '_on_g20_toggle.png')
        plt.clf()

        dac_ibi_g20_63 = 'dac_ibi_g20 63; hg_lg_c2 220; on_g20 True; shaper_constants 14_6_15; sw_ibo_g20 False'
        dac_ibi_g20_30 = 'dac_ibi_g20 30; hg_lg_c2 220; on_g20 True; shaper_constants 14_6_15; sw_ibo_g20 False'
        dac_ibi_g20_12 = 'dac_ibi_g20 12; hg_lg_c2 220; on_g20 True; shaper_constants 14_6_15; sw_ibo_g20 False'
        dac_ibi_g20_0 = 'dac_ibi_g20 0; hg_lg_c2 220; on_g20 True; shaper_constants 14_6_15; sw_ibo_g20 False'
        plts = [dac_ibi_g20_63, dac_ibi_g20_30, dac_ibi_g20_12, dac_ibi_g20_0]
        l1, = plt.plot(dat[dac_ibi_g20_63])
        l2, = plt.plot(dat[dac_ibi_g20_30])
        l3, = plt.plot(dat[dac_ibi_g20_12])
        l4, = plt.plot(dat[dac_ibi_g20_0])
        plt.title('Channel' + str(channel+1) +': Dac_ibi_g20: 63, 30, 12, and 0')
        plt.xlabel('ns')
        legend1 = plt.legend([l1,l2,l3,l4], ['dac_ibi_g20 63 (default)', 'dac_ibi_g20 30', 'dac_ibi_g20 12', 'dac_ibi_g20 0'], loc = 'upper right')
        legend2 = plt.legend([l1,l2,l3,l4], ["Rise time: " + str(st[plts[i]]) + " ns" for i in range(len(plts))], loc = 'upper right', bbox_to_anchor= (1.0,0.5))
        legend3 = plt.legend([l1,l2,l3,l4], ["Height: " + str(int(hgt[plts[i]])) + " counts" for i in range(len(plts))], loc = 'upper right', bbox_to_anchor= (1.0,0.75))
        ax1 = plt.gca().add_artist(legend1)
        ax2 = plt.gca().add_artist(legend2)
        ax3 = plt.gca().add_artist(legend3)
        plt.savefig(saveDir + 'channel' + str(channel+1) +'_dac_ibi_g20_range.png')
        plt.clf()  

        sw_ibo_g20_false = 'dac_ibi_g20 63; hg_lg_c2 220; on_g20 True; shaper_constants 14_6_15; sw_ibo_g20 False'
        sw_ibo_g20_true = 'dac_ibi_g20 63; hg_lg_c2 220; on_g20 True; shaper_constants 14_6_15; sw_ibo_g20 True'
        l1, = plt.plot(dat[sw_ibo_g20_false])
        l2, = plt.plot(dat[sw_ibo_g20_true])
        plt.title('Channel' + str(channel+1) +': sw_ibo_g20 Toggle')
        plt.xlabel('ns')
        legend1 = plt.legend([l1,l2], ["sw_ibi_g20 False (default)", 'sw_ibi_g20 True'], loc = 'upper right')
        legend2 = plt.legend([l1,l2], ["Rise time: " + str(st[sw_ibo_g20_false]) + " ns", 'Rise time: ' + str(st[sw_ibo_g20_true]) + " ns"], loc = 'upper right', bbox_to_anchor= (1.0,0.7))
        legend3 = plt.legend([l1,l2], ["Height: " + str(int(hgt[sw_ibo_g20_false])) + " counts", 'Height: ' + str(int(hgt[sw_ibo_g20_true])) + " counts"], loc = 'upper right', bbox_to_anchor= (1.0,0.85))
        ax1 = plt.gca().add_artist(legend1)
        ax2 = plt.gca().add_artist(legend2)
        ax3 = plt.gca().add_artist(legend3)
        plt.savefig(saveDir + 'channel' + str(channel+1) + '_sw_ibo_g20_toggle.png')
        plt.clf()

    #plt.plot(dat11['dac_ibi_g20 63; hg_lg_c2 220; on_g20 True; shaper_constants 14_6_15; sw_ibo_g20 False'])
    #plt.plot(dat11['dac_ibi_g20 63; hg_lg_c2 220; on_g20 True; shaper_constants 14_6_15; sw_ibo_g20 False'])
    #combos = dat11.keys()
    #print(combos)
    # print(combos)
    # for attr in cut_attrs:
    #   param_vals = np.unique(f.attrs[attr])
    #   for val in param_vals:
    #     for setting in [prod for prod in combos if val in prod]:
    #       plt.plot(dat11[setting][0:1000], label = setting.replace(attr+" "+val,'').replace("; ",''))
    #     plt.legend()
    #     plt.xlabel('ns')
    #     plt.title("Attribute "+attr+" fixed at "+val)
    #     print('should be saving plots to ' + saveDir)
    #     plt.savefig(saveDir + attr + "_" + val + ".png")
    #     plt.clf()
    #plt.plot(dat)

