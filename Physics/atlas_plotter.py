import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
from scipy import stats
from helperfunctions import *



#plotsDir = checkDir(getRootDir()+"/Plotting/PlotsSaved/Pulses/")
plotsDir = checkDir(getRootDir()+"/Plotting/PlotsSaved/Pulses/OFCs/")
pedestalDir = checkDir(getRootDir()+"/Plotting/PlotsSaved/Pedestal/")

colors = 'bcgrymk'
c_iter = 0

enable_latex = True

#ATLAS Style
if(enable_latex): mpl.rcParams['text.usetex'] = True
if(enable_latex): mpl.rcParams['font.size'] = 12
if(enable_latex): mpl.rcParams['text.latex.preamble'] = [r'\usepackage{helvet}',
                                        r'\usepackage{sansmath}',
                                        r'\setlength{\parindent}{0pt}',
                                        r'\sansmath',
                                        r'\boldmath']

##########################

def new_log_formatter(x, p): #Formats log plot scales with smaller negaative signs
    """ Own formatting function """
    return r"$10$\textsuperscript{%i}" % np.log10(x)  #  raw string to avoid "\\"


def next_color():
  global colors
  global c_iter
  color = colors[c_iter]
  c_iter = (c_iter + 1)%len(colors)
  return color

def reset_color():
  global c_iter
  c_iter = 0 

def atlabel(adc, channel, x = 0.05, y = 0.95):
  atstyle.use_atlas_style()
#  atstyle.draw_atlas_label(x, y, status="Upgrade", desc="COLUTAv2\nADC "+str(adc)+", Channel "+str(channel))


def make_xlabel(plt, label):
    if(enable_latex): plt.xlabel(r'\textbf{'+label+'}', x=1.0, ha='right', size=12, labelpad=0)
    else: plt.xlabel(r''+label+'', x=1.0, ha='right', size=12, labelpad=0)

def make_ylabel(plt, label):
    #plt.ylabel(r'\textbf{'+label+'}', y=1.0, ha='right', size=12, labelpad=0)
    if(enable_latex): plt.ylabel(r'\textbf{'+label+'}', y=1.0, ha='right', size=12, labelpad=0)
    else: plt.ylabel(r''+label+'', y=1.0, ha='right', size=12, labelpad=0)

def drawLabel(ax, adc, channel, filter_name, filter_info, run, pos):
  global _atlas_label
  if(enable_latex): label = (r'\textbf{\textit{ATLAS}} Upgrade')
  else: label = (r'ATLAS Upgrade')
  label = label+('\n'+r'COLUTAv2' 
                +'\n'+r''+run.replace("/","").replace("_"," ") 
                +'\n'+r'ADC '+str(adc)+', Channel '+str(channel) 
                #+'\n'+r''+str(filter_name)+': 5 Samples' 
                +'\n'+r''+str(filter_name)+': '+str(filter_info) 
                #+'\n'+r'Averaged') 
                #+'\n'+r'Single Pulse') 
                )

  x = 0.05
  y = 0.05
  va = 'bottom'
  ha = 'left'
  if 'b' not in pos: #Top
    y = 0.95
    va = 'top'
  if 'l' not in pos: #Right
    x = 0.95
    ha = 'right'
  plt.text(x, y,label, size=14, transform=ax.transAxes, multialignment=ha, verticalalignment=va, horizontalalignment=ha)



def drawLabel_pedestal(ax, adc, channel, run, pos):
  global _atlas_label
  if(enable_latex): label = (r'\textbf{\textit{ATLAS}} Upgrade')
  else: label = (r'ATLAS Upgrade')
  label = label+('\n'+r'COLUTAv2'
                +'\n'+r'ADC '+str(adc)+', Channel '+str(channel) 
                #+'\n'+r'Run 0994'
                +'\n'+run
                +'\n'+r'Pedestal Run')
                #+'\n'+r'Single Pulse')

  x = 0.05
  y = 0.05
  va = 'bottom'
  ha = 'left'
  if 'b' not in pos: #Top
    y = 0.95
    va = 'top'
  if 'l' not in pos: #Right
    x = 0.95
    ha = 'right'
  plt.text(x, y,label, size=14, transform=ax.transAxes, multialignment=ha, verticalalignment=va, horizontalalignment=ha)








def get_ref_noise(ref_energies):
  ref_noise = np.zeros(len(ref_energies))
  for i in range(len(ref_energies)):
    ref_noise[i] = np.sqrt( np.square(0.1/np.sqrt(ref_energies[i])) + np.square(0.0025) + np.square(0.0546/ref_energies[i]))
  return ref_noise


#Data = 2D Array of shape (n_measurements, n_samples)
def make_pedestal_plots(data, adc, channel, run):
  
  data = data.flatten()
  data_rms = get_rms(data)


  #Histo
  ax=plt.subplot(111)
#  plt.grid(linestyle='--')
  plt.minorticks_on()
  #plt.hist(data, bins = np.arange(min(data), max(data)+1, 1), normed=True)
  plt.hist(data, bins = np.arange(min(data), max(data)+1, 1), normed=False)
  xt = plt.xticks()[0]  
  xmin, xmax = min(xt), max(xt)  
  lnspc = np.linspace(xmin, xmax, len(data))
  m, s = stats.norm.fit(data) # get mean and standard deviation  
  pdf_g = stats.norm.pdf(lnspc, m, s) # now get theoretical values in our interval  
#  plt.plot(lnspc, pdf_g, label="Gaussian Fit") # plot it
  make_xlabel(plt, r'ADC Counts')
  #make_ylabel(plt, r'Normalized Entries')
  make_ylabel(plt, r'Entries')
  #ymin, ymax = plt.ylim()
  #plt.ylim(ymin-0.3*(ymax-ymin), ymax+0.5*(ymax-ymin))
  #plt.ylim(0.5, 4.5)
  drawLabel_pedestal(ax, adc+1, channel+1, short(run).replace("/","").replace("_"," "), 'tl')
  fmt=mpl.ticker.StrMethodFormatter("{x:,.0f}")
  fmt_x=mpl.ticker.StrMethodFormatter("{x:,.0f}")
  ax.yaxis.set_major_formatter(fmt)
  ax.xaxis.set_major_formatter(fmt_x)
  ax.tick_params(direction='in', which='both')
  ax.get_xaxis().set_ticks_position('both')
  ax.get_yaxis().set_ticks_position('both')
  locs, labels = plt.xticks()
  plt.legend(loc=1, frameon=False)
  plt.text(0.95, 0.85, r'$\mu = $'+str(np.round(np.mean(data), 3)), horizontalalignment='right', verticalalignment='center', transform=ax.transAxes) 
  plt.text(0.95, 0.75, r'$\sigma = $'+str(np.round(np.std(data), 3)), horizontalalignment='right', verticalalignment='center', transform=ax.transAxes) 
#  plt.text(0.95, 0.65, 'rms = '+str(np.round(data_rms, 3)), horizontalalignment='right', verticalalignment='center', transform=ax.transAxes) 
#  plt.show()
  plt.savefig(checkDir(pedestalDir+short(run)+"/Pedestal_ADC"+str(adc+1)+"_Channel"+str(channel+1)+".png"))
  print("File Saved: "+pedestalDir+short(run)+"/Pedestal_ADC"+str(adc+1)+"_Channel"+str(channel+1)+".png")
  plt.clf()
  plt.cla()
  plt.close()
  reset_color()

  """
  ax = plt.subplot(111) 
  xt = plt.xticks()[0]  
  xmin, xmax = min(xt), max(xt)  
  lnspc = np.linspace(xmin, xmax, len(amplitudes[amp][nSamples][phase]))
  m, s = stats.norm.fit(amplitudes[amp][nSamples][phase]) # get mean and standard deviation  
  pdf_g = stats.norm.pdf(lnspc, m, s) # now get theoretical values in our interval  
  ax.plot(lnspc, pdf_g, label="Norm") # plot it
  ax.text(0.8, 0.9, r'$\mu$: '+str(np.round(m, 4)), verticalalignment='top', horizontalalignment='left', transform=ax.transAxes)
  ax.text(0.8, 0.8, r'$\sigma$: '+str(np.round(s, 4)), verticalalignment='top', horizontalalignment='left', transform=ax.transAxes)
  plt.title('Energy, Samples: '+str(nSamples)+' Phase: '+str(phase))
  plt.ylabel('Normalized Entries')
  plt.xlabel('Energy')
  plt.savefig(''+checkDir(plotsDir+'Histograms/')+'Energy_S'+str(nSamples)+'_P'+str(phase)+'.png')
  print('Created File: '+checkDir(plotsDir+'Histograms/')+'Energy_S'+str(nSamples)+'_P'+str(phase)+'.png')
  plt.clf()
  """




#def makePlots(adc, channel, n_samples, gain, amps, awg_amps, run, adc_freq, awg_freq):
def makePlots(data, adc, channel, setup_params, run):
  adc += 1
  channel += 1
  amps = setup_params["amps"]
  awg_amps = setup_params["awg_amps"]
  run = setup_params["run"]
  adc_freq = setup_params["adc_freq"]
  awg_freq = setup_params["awg_freq"]
  filter_name = setup_params["filter_name"]
  filter_info = setup_params["filter_info"]

  amps = amps.astype("U13")

  gains = ""
  for g in data:
    gains += "_"+g["gain"]
  gains = gains[1:]

  filter_name_subdir = filter_name.replace(" ","")
  filter_info_subdir = filter_info.replace(" ","")


  #Energy Resolution
  #Sigma E
  ax=plt.subplot(111)
  plt.grid(linestyle='--')
  plt.minorticks_on()
  for points in data:
    plt.errorbar(points["energy_mean"], points["energy_std"], yerr=points["energy_std"]/np.sqrt(145), fmt='o-', color=next_color(), linewidth=0.4, markersize = 5, capsize=2, label=points["gain"])
  #make_xlabel(plt, r'$\textbf{E}$ from OFCs [ADC Counts]')
  #make_ylabel(plt, r'$\sigma_{\textbf{\tiny E}}$ from OFCs [ADC Counts]')
  #TEST: Single sample (peak)
  if(enable_latex):
    make_xlabel(plt, r'$\textbf{E}$ from OFCs [ADC Counts]')
    make_ylabel(plt, r'$\sigma_{\textbf{\tiny E}}$ from OFCs [ADC Counts]')
  else:
    make_xlabel(plt, r'$E$ from OFCs [ADC Counts]')
    make_ylabel(plt, r'$\sigma E$ from OFCs [ADC Counts]')
  ymin, ymax = plt.ylim()
  #plt.ylim(ymin-0.3*(ymax-ymin), ymax+0.5*(ymax-ymin))
  plt.ylim(0., 8.5)
  drawLabel(ax, adc, channel, filter_name, filter_info, run, 'tl')
  fmt=mpl.ticker.StrMethodFormatter("{x}")
  ax.yaxis.set_major_formatter(fmt)
  ax.xaxis.set_major_formatter(fmt)
  ax.tick_params(direction='in', which='both')
  ax.get_xaxis().set_ticks_position('both')
  ax.get_yaxis().set_ticks_position('both')
  locs, labels = plt.xticks()
  plt.legend(loc=1, frameon=False)
  #plt.show(block=True)
  plt.savefig(checkDir(plotsDir+short(run)+"/Performance/"+filter_name_subdir+"/"+gains+"/"+filter_info_subdir+"/Sigma_E_ADC"+str(adc)+"_Channel"+str(channel)+"_"+gains+".png"))
  print("File Saved: "+plotsDir+short(run)+"/Performance/"+filter_name_subdir+"/"+gains+"/"+filter_info_subdir+"/Sigma_E_ADC"+str(adc)+"_Channel"+str(channel)+"_"+gains+".png")
  plt.clf()
  plt.cla()
  plt.close()
  reset_color()



  #Energy Resolution From Peak
  #Sigma E Peak
  ax=plt.subplot(111)
  plt.grid(linestyle='--')
  plt.minorticks_on()
  for points in data:
    plt.errorbar(points["energy_mean_peak"], points["energy_std_peak"], yerr=points["energy_std_peak"]/np.sqrt(145), fmt='o-', color=next_color(), linewidth=0.4, markersize = 5, capsize=2, label=points["gain"])
  #make_xlabel(plt, r'$\textbf{E}$ from OFCs [ADC Counts]')
  #make_ylabel(plt, r'$\sigma_{\textbf{\tiny E}}$ from OFCs [ADC Counts]')
  #TEST: Single sample (peak)
  #make_xlabel(plt, r'$\textbf{E}$ from Peak [ADC Counts]')
  #make_ylabel(plt, r'$\sigma_{\textbf{\tiny E}}$ from Peak [ADC Counts]')
  if(enable_latex):
    make_xlabel(plt, r'$\textbf{E}$ from Peak [ADC Counts]')
    make_ylabel(plt, r'$\sigma_{\textbf{\tiny E}}$ from Peak [ADC Counts]')
  else:
    make_xlabel(plt, r'$E$ from Peak [ADC Counts]')
    make_ylabel(plt, r'$\sigma E$ from Peak [ADC Counts]')
  ymin, ymax = plt.ylim()
  #plt.ylim(ymin-0.3*(ymax-ymin), ymax+0.5*(ymax-ymin))
  plt.ylim(0., 8.5)
  #drawLabel(ax, adc, channel, filter_name, filter_info, 'tl')
  drawLabel(ax, adc, channel, "Single Point", "Peak", run, 'tl')
  fmt=mpl.ticker.StrMethodFormatter("{x}")
  ax.yaxis.set_major_formatter(fmt)
  ax.xaxis.set_major_formatter(fmt)
  ax.tick_params(direction='in', which='both')
  ax.get_xaxis().set_ticks_position('both')
  ax.get_yaxis().set_ticks_position('both')
  locs, labels = plt.xticks()
  plt.legend(loc=1, frameon=False)
  #plt.show(block=True)
  plt.savefig(checkDir(plotsDir+short(run)+"/Performance/"+filter_name_subdir+"/"+gains+"/"+filter_info_subdir+"/Sigma_E_Peak_ADC"+str(adc)+"_Channel"+str(channel)+"_"+gains+".png"))
  print("File Saved: "+plotsDir+short(run)+"/Performance/"+filter_name_subdir+"/"+gains+"/"+filter_info_subdir+"/Sigma_E_Peak_ADC"+str(adc)+"_Channel"+str(channel)+"_"+gains+".png")
  plt.clf()
  plt.cla()
  plt.close()
  reset_color()





  #Energy Reconstruction
  #OFC E vs Peak E
  ax=plt.subplot(111)
  plt.grid(linestyle='--')
  plt.minorticks_on()
  for points in data:
    plt.errorbar(points["energy_mean_peak"], np.subtract(points["energy_mean"], points["energy_mean_peak"]), yerr=points["energy_std"], fmt='o-', color=next_color(), linewidth=0.4, markersize = 5, capsize=2, label=points["gain"])
  #make_xlabel(plt, r'$\textbf{E}$ from OFCs [ADC Counts]')
  #make_ylabel(plt, r'$\sigma_{\textbf{\tiny E}}$ from OFCs [ADC Counts]')
  #TEST: Single sample (peak)
  #make_xlabel(plt, r'$\textbf{E}$ from peak [ADC Counts]')
  #make_ylabel(plt, r'$\textbf{E}$ from OFCs, deviation from y=x [ADC Counts]')
  if(enable_latex):
    make_xlabel(plt, r'$\textbf{E}$ from peak [ADC Counts]')
    make_ylabel(plt, r'$\textbf{E}$ from OFCs, deviation from y=x [ADC Counts]')
  else:
    make_xlabel(plt, r'$E$ from peak [ADC Counts]')
    make_ylabel(plt, r'$E$ from OFCs, deviation from y=x [ADC Counts]')
  ymin, ymax = plt.ylim()
  plt.ylim(ymin-0.3*(ymax-ymin), ymax+0.5*(ymax-ymin))
  drawLabel(ax, adc, channel, "Single Point", "Peak", run, 'tl')
  fmt=mpl.ticker.StrMethodFormatter("{x}")
  ax.yaxis.set_major_formatter(fmt)
  ax.xaxis.set_major_formatter(fmt)
  ax.tick_params(direction='in', which='both')
  ax.get_xaxis().set_ticks_position('both')
  ax.get_yaxis().set_ticks_position('both')
  locs, labels = plt.xticks()
  plt.legend(loc=1, frameon=False)
  #plt.show(block=True)
  plt.savefig(checkDir(plotsDir+short(run)+"/Performance/"+filter_name_subdir+"/"+gains+"/"+filter_info_subdir+"/E_Reco_ADC"+str(adc)+"_Channel"+str(channel)+"_"+gains+".png"))
  print("File Saved: "+plotsDir+short(run)+"/Performance/"+filter_name_subdir+"/"+gains+"/"+filter_info_subdir+"/E_Reco_ADC"+str(adc)+"_Channel"+str(channel)+"_"+gains+".png")
  plt.clf()
  plt.cla()
  plt.close()
  reset_color()


  #Energy Resolution Ratio
  #Sigma E / E
  ax=plt.subplot(111)
  plt.grid(linestyle='--')
  plt.minorticks_on()
  for points in data:
    #plt.errorbar(points["energy_mean"], np.divide(points["energy_std"], points["energy_mean"]), yerr=np.divide(points["energy_std"]/np.sqrt(145), points["energy_mean"]), fmt='o-', color=next_color(), linewidth=0.4, markersize = 2, capsize=2, label=points["gain"])
    plt.errorbar(points["energy_mean"], np.divide(points["energy_std"], points["energy_mean"]), yerr=getDivisionUncertainty(points["energy_std"], points["energy_mean"], np.divide(points["energy_std"], np.sqrt(145)), points["energy_std"]), xerr=points["energy_std"], fmt='o-', color=next_color(), linewidth=0.4, markersize = 2, capsize=2, label=points["gain"])
  #make_ylabel(plt, r'$\sigma_{\textbf{\tiny E}/E}$, from OFCs')
  #make_xlabel(plt, r'$\textbf{E}$ from OFCs [ADC Counts]')
  if(enable_latex):
    make_ylabel(plt, r'$\sigma_{\textbf{\tiny E}/E}$, from OFCs')
    make_xlabel(plt, r'$\textbf{E}$ from OFCs [ADC Counts]')
  else:
    make_ylabel(plt, r'$\sigma E/E$, from OFCs')
    make_xlabel(plt, r'$E$ from OFCs [ADC Counts]')
  #TEST: Single sample (peak)
  #make_ylabel(plt, r'$\sigma_{\textbf{\tiny E}/E}$, from peak')
  #make_xlabel(plt, r'$\textbf{E}$ from peak [ADC Counts]')
  plt.yscale("log", nonposy='clip')
  plt.xscale("log", nonposx='clip')
  ymin, ymax = plt.ylim()
  xmin, xmax = plt.xlim()
  plt.ylim(ymin/10, ymax*20000)
  plt.xlim(xmin, xmax*5)
  #Calorimeter Expected Error 
  ref_energies = np.linspace(xmin, xmax*5, 100)
  ref_noise = get_ref_noise(ref_energies)
#  plt.plot(ref_energies, ref_noise, '-', color='k', label="Calorimeter Expected Error") 
  drawLabel(ax, adc, channel, filter_name, filter_info, run, 'tl')
  fmt=mpl.ticker.FuncFormatter(new_log_formatter)
  ax.yaxis.set_major_formatter(fmt)
  ax.xaxis.set_major_formatter(fmt)
  ax.tick_params(direction='in', which='both')
  ax.get_xaxis().set_ticks_position('both')
  ax.get_yaxis().set_ticks_position('both')
  locs, labels = plt.xticks()
  plt.legend(loc=1, frameon=False)
  #plt.show(block=True)
  plt.savefig(checkDir(plotsDir+short(run)+"/Performance/"+filter_name_subdir+"/"+gains+"/"+filter_info_subdir+"/Sigma_E_Ratio_ADC"+str(adc)+"_Channel"+str(channel)+"_"+gains+".png"))
  print("File Saved: "+plotsDir+short(run)+"/Performance/"+filter_name_subdir+"/"+gains+"/"+filter_info_subdir+"/Sigma_E_Ratio_ADC"+str(adc)+"_Channel"+str(channel)+"_"+gains+".png")
  plt.clf()
  plt.cla()
  plt.close()
  reset_color()


  #Timing Resolution
  #Sigma t
  
  ax=plt.subplot(111)
  plt.grid(linestyle='--')
  plt.minorticks_on()
  for points in data:
    plt.errorbar(points["energy_mean"], points["timing_std"], yerr=points["timing_std"]/np.sqrt(145), fmt='o-', color=next_color(), linewidth=0.4, markersize = 2, capsize=2, label=points["gain"])
  #make_ylabel(plt, r'$\sigma_{\textbf{t}}$, from OFCs [ns]')
  #make_xlabel(plt, r'$\textbf{E}$ from OFCs [ADC Counts]')
  if(enable_latex):
    make_ylabel(plt, r'$\sigma_{\textbf{t}}$, from OFCs [ns]')
    make_xlabel(plt, r'$\textbf{E}$ from OFCs [ADC Counts]')
  else:
    make_ylabel(plt, r'$\sigma t}$, from OFCs [ns]')
    make_xlabel(plt, r'$E$ from OFCs [ADC Counts]')
  plt.yscale("log", nonposy='clip')
  plt.xscale("log", nonposx='clip')
  ymin, ymax = plt.ylim()
  xmin, xmax = plt.xlim()
  plt.ylim(ymin/10, ymax*50)
  plt.xlim(xmin, xmax*5)
  drawLabel(ax, adc, channel, filter_name, filter_info, run, 'tl')
  fmt=mpl.ticker.FuncFormatter(new_log_formatter)
  ax.yaxis.set_major_formatter(fmt)
  ax.xaxis.set_major_formatter(fmt)
  ax.tick_params(direction='in', which='both')
  ax.get_xaxis().set_ticks_position('both')
  ax.get_yaxis().set_ticks_position('both')
  locs, labels = plt.xticks()
  plt.legend(loc=1, frameon=False)
  #plt.show(block=True)
  plt.savefig(checkDir(plotsDir+short(run)+"/Performance/"+filter_name_subdir+"/"+gains+"/"+filter_info_subdir+"/Sigma_T_ADC"+str(adc)+"_Channel"+str(channel)+"_"+gains+".png"))
  print("File Saved: "+plotsDir+short(run)+"/Performance/"+filter_name_subdir+"/"+gains+"/"+filter_info_subdir+"/Sigma_T_ADC"+str(adc)+"_Channel"+str(channel)+"_"+gains+".png")
  plt.clf()
  plt.cla()
  plt.close()
  reset_color()

  #Timing Correlation
  #t16 vs t15
  ax=plt.subplot(111)
  plt.grid(linestyle='--')
  plt.minorticks_on()
  for points in data:
    plt.scatter(points["timing"], points["timing_p16"], c=next_color(), label=points["gain"])
  if(enable_latex):
    make_ylabel(plt, r'$t_{16}$, from OFCs [ns]')
    make_xlabel(plt, r'$t_{15}$, from OFCs [ns]')
  else:
    make_ylabel(plt, r'$t 16$, from OFCs [ns]')
    make_xlabel(plt, r'$t 15$ from OFCs [ns]')
  ymin, ymax = plt.ylim()
  plt.ylim(ymin - 0.15*(ymax-ymin), ymax+0.5*(ymax-ymin))
  drawLabel(ax, adc, channel, filter_name, filter_info, run, 'tl')
  fmtx=mpl.ticker.StrMethodFormatter("{x:,.2f}")
  fmty=mpl.ticker.StrMethodFormatter("{x:,.2f}")
  ax.yaxis.set_major_formatter(fmty)
  ax.xaxis.set_major_formatter(fmtx)
  ax.tick_params(direction='in', which='both')
  ax.get_xaxis().set_ticks_position('both')
  ax.get_yaxis().set_ticks_position('both')
  locs, labels = plt.xticks()
  plt.legend(loc=1, frameon=False)
  #plt.show(block=True)
  plt.savefig(checkDir(plotsDir+short(run)+"/Performance/"+filter_name_subdir+"/"+gains+"/"+filter_info_subdir+"/T_2D_Correlation_ADC"+str(adc)+"_Channel"+str(channel)+"_"+gains+".png"))
  print("File Saved: "+plotsDir+short(run)+"/Performance/"+filter_name_subdir+"/"+gains+"/"+filter_info_subdir+"/T_2D_Correlation_ADC"+str(adc)+"_Channel"+str(channel)+"_"+gains+".png")
  plt.clf()
  plt.cla()
  plt.close()
  reset_color()

  #Distribution of t15
  ax=plt.subplot(111)
  plt.minorticks_on()
  try:
    plt.hist(points["timing"], bins = 25, normed=False)
  except ValueError:
    pass
  make_xlabel(plt, r'$t_{15}$, from OFCs [ns]')
  make_ylabel(plt, r'Entries')
  drawLabel(ax, adc, channel, filter_name, filter_info, run, 'tl')
  fmt=mpl.ticker.StrMethodFormatter("{x:,.2f}")
  fmt_x=mpl.ticker.StrMethodFormatter("{x:,.2f}")
  ax.yaxis.set_major_formatter(fmt)
  ax.xaxis.set_major_formatter(fmt_x)
  ax.tick_params(direction='in', which='both')
  ax.get_xaxis().set_ticks_position('both')
  ax.get_yaxis().set_ticks_position('both')
  locs, labels = plt.xticks()
  plt.legend(loc=1, frameon=False)
  plt.text(0.95, 0.85, r'$\mu = $'+str(np.round(np.mean(points["timing"]), 5)), horizontalalignment='right', verticalalignment='center', transform=ax.transAxes) 
  plt.text(0.95, 0.75, r'$\sigma = $'+str(np.round(np.std(points["timing"]), 5)), horizontalalignment='right', verticalalignment='center', transform=ax.transAxes) 
#  plt.text(0.95, 0.65, 'rms = '+str(np.round(points["timing"]_rms, 3)), horizontalalignment='right', verticalalignment='center', transform=ax.transAxes) 
#  plt.show()
  plt.savefig(checkDir(plotsDir+short(run)+"/Performance/"+filter_name_subdir+"/"+gains+"/"+filter_info_subdir+"/Hist_T15_ADC"+str(adc)+"_Channel"+str(channel)+"_"+gains+".png"))
  print("File Saved: "+plotsDir+short(run)+"/Performance/"+filter_name_subdir+"/"+gains+"/"+filter_info_subdir+"/Hist_T15_ADC"+str(adc)+"_Channel"+str(channel)+"_"+gains+".png")
  plt.clf()
  plt.cla()
  plt.close()
  reset_color()

  #Distribution of t16
  ax=plt.subplot(111)
  plt.minorticks_on()
  try:
    plt.hist(points["timing_p16"], bins = 25, normed=False)
  except ValueError:
    pass
  make_xlabel(plt, r'$t_{16}$, from OFCs [ns]')
  make_ylabel(plt, r'Entries')
  drawLabel(ax, adc, channel, filter_name, filter_info, run, 'tl')
  fmt=mpl.ticker.StrMethodFormatter("{x:,.2f}")
  fmt_x=mpl.ticker.StrMethodFormatter("{x:,.2f}")
  ax.yaxis.set_major_formatter(fmt)
  ax.xaxis.set_major_formatter(fmt_x)
  ax.tick_params(direction='in', which='both')
  ax.get_xaxis().set_ticks_position('both')
  ax.get_yaxis().set_ticks_position('both')
  locs, labels = plt.xticks()
  plt.legend(loc=1, frameon=False)
  plt.text(0.95, 0.85, r'$\mu = $'+str(np.round(np.mean(points["timing_p16"]), 5)), horizontalalignment='right', verticalalignment='center', transform=ax.transAxes) 
  plt.text(0.95, 0.75, r'$\sigma = $'+str(np.round(np.std(points["timing_p16"]), 5)), horizontalalignment='right', verticalalignment='center', transform=ax.transAxes) 
#  plt.text(0.95, 0.65, 'rms = '+str(np.round(points["timing_p16"]_rms, 3)), horizontalalignment='right', verticalalignment='center', transform=ax.transAxes) 
#  plt.show()
  plt.savefig(checkDir(plotsDir+short(run)+"/Performance/"+filter_name_subdir+"/"+gains+"/"+filter_info_subdir+"/Hist_T16_ADC"+str(adc)+"_Channel"+str(channel)+"_"+gains+".png"))
  print("File Saved: "+plotsDir+short(run)+"/Performance/"+filter_name_subdir+"/"+gains+"/"+filter_info_subdir+"/Hist_T16_ADC"+str(adc)+"_Channel"+str(channel)+"_"+gains+".png")
  plt.clf()
  plt.cla()
  plt.close()
  reset_color()

  #Distribution of t16-t15
  ax=plt.subplot(111)
  plt.minorticks_on()
  try:
    plt.hist(np.subtract(points["timing_p16"], points["timing"]), bins = 25, normed=False)
  except ValueError:
    pass
  make_xlabel(plt, r'$t_{16}-t_{15}$, from OFCs [ns]')
  make_ylabel(plt, r'Entries')
  drawLabel(ax, adc, channel, filter_name, filter_info, run, 'tl')
  fmt=mpl.ticker.StrMethodFormatter("{x:,.2f}")
  fmt_x=mpl.ticker.StrMethodFormatter("{x:,.2f}")
  ax.yaxis.set_major_formatter(fmt)
  ax.xaxis.set_major_formatter(fmt_x)
  ax.tick_params(direction='in', which='both')
  ax.get_xaxis().set_ticks_position('both')
  ax.get_yaxis().set_ticks_position('both')
  locs, labels = plt.xticks()
  plt.legend(loc=1, frameon=False)
  plt.text(0.95, 0.85, r'$\mu = $'+str(np.round(np.mean(np.subtract(points["timing_p16"], points["timing"])), 5)), horizontalalignment='right', verticalalignment='center', transform=ax.transAxes) 
  plt.text(0.95, 0.75, r'$\sigma = $'+str(np.round(np.std(np.subtract(points["timing_p16"], points["timing"])), 5)), horizontalalignment='right', verticalalignment='center', transform=ax.transAxes) 
#  plt.text(0.95, 0.65, 'rms = '+str(np.round(np.subtract(points["timing_p16"], points["timing"])_rms, 3)), horizontalalignment='right', verticalalignment='center', transform=ax.transAxes) 
#  plt.show()
  plt.savefig(checkDir(plotsDir+short(run)+"/Performance/"+filter_name_subdir+"/"+gains+"/"+filter_info_subdir+"/Hist_Timing_Difference_ADC"+str(adc)+"_Channel"+str(channel)+"_"+gains+".png"))
  print("File Saved: "+plotsDir+short(run)+"/Performance/"+filter_name_subdir+"/"+gains+"/"+filter_info_subdir+"/Hist_Timing_Difference_ADC"+str(adc)+"_Channel"+str(channel)+"_"+gains+".png")
  plt.clf()
  plt.cla()
  plt.close()
  reset_color()


  #Distribution of t16+t15
  ax=plt.subplot(111)
  plt.minorticks_on()
  try:
    plt.hist(np.add(points["timing_p16"], points["timing"]), bins = 25, normed=False)
  except ValueError:
    pass
  make_xlabel(plt, r'$t_{16}+t_{15}$, from OFCs [ns]')
  make_ylabel(plt, r'Entries')
  drawLabel(ax, adc, channel, filter_name, filter_info, run, 'tl')
  fmt=mpl.ticker.StrMethodFormatter("{x:,.2f}")
  fmt_x=mpl.ticker.StrMethodFormatter("{x:,.2f}")
  ax.yaxis.set_major_formatter(fmt)
  ax.xaxis.set_major_formatter(fmt_x)
  ax.tick_params(direction='in', which='both')
  ax.get_xaxis().set_ticks_position('both')
  ax.get_yaxis().set_ticks_position('both')
  locs, labels = plt.xticks()
  plt.legend(loc=1, frameon=False)
  plt.text(0.95, 0.85, r'$\mu = $'+str(np.round(np.mean(np.add(points["timing_p16"], points["timing"])), 5)), horizontalalignment='right', verticalalignment='center', transform=ax.transAxes) 
  plt.text(0.95, 0.75, r'$\sigma = $'+str(np.round(np.std(np.add(points["timing_p16"], points["timing"])), 5)), horizontalalignment='right', verticalalignment='center', transform=ax.transAxes) 
#  plt.text(0.95, 0.65, 'rms = '+str(np.round(np.add(points["timing_p16"], points["timing"])_rms, 3)), horizontalalignment='right', verticalalignment='center', transform=ax.transAxes) 
#  plt.show()
  plt.savefig(checkDir(plotsDir+short(run)+"/Performance/"+filter_name_subdir+"/"+gains+"/"+filter_info_subdir+"/Hist_Timing_Sum_ADC"+str(adc)+"_Channel"+str(channel)+"_"+gains+".png"))
  print("File Saved: "+plotsDir+short(run)+"/Performance/"+filter_name_subdir+"/"+gains+"/"+filter_info_subdir+"/Hist_Timing_Sum_ADC"+str(adc)+"_Channel"+str(channel)+"_"+gains+".png")
  plt.clf()
  plt.cla()
  plt.close()
  reset_color()



  #Linearity
  #E OFC vs E AWG 
  ax = plt.subplot(111)
  plt.grid(linestyle='--')
  plt.minorticks_on()
  for points in data:
    plt.errorbar(awg_amps, points["energy_mean"], yerr=points["energy_std"], fmt='o-', color=next_color(), linewidth=0.4, markersize = 2, capsize=2, label=points["gain"])
  #make_ylabel(plt, r'$\textbf{E}_{\textbf{OFC}}$ [ADC Counts]')
  #make_xlabel(plt, r'$\textbf{E}_{\textbf{AWG}}$ [A.U.]')
  if(enable_latex):
    make_ylabel(plt, r'$\textbf{E}_{\textbf{OFC}}$ [ADC Counts]')
    make_xlabel(plt, r'$\textbf{E}_{\textbf{AWG}}$ [A.U.]')
  else:
    make_ylabel(plt, r'$E OFC$ [ADC Counts]')
    make_xlabel(plt, r'$E AWG$ [A.U.]')
  plt.yscale("log", nonposy='clip')
  plt.xscale("log", nonposx='clip')
  ymin, ymax = plt.ylim()
  xmin, xmax = plt.xlim()
  plt.ylim(ymin/10, ymax*50)
  plt.xlim(xmin, xmax*5)
  drawLabel(ax, adc, channel, filter_name, filter_info, run, 'tl')
  fmt=mpl.ticker.FuncFormatter(new_log_formatter)
  ax.yaxis.set_major_formatter(fmt)
  ax.xaxis.set_major_formatter(fmt)
  ax.tick_params(direction='in', which='both')
  ax.get_xaxis().set_ticks_position('both')
  ax.get_yaxis().set_ticks_position('both')
  locs, labels = plt.xticks()
  plt.legend(loc=1, frameon=False)
  #plt.show(block=True)
  plt.savefig(checkDir(plotsDir+short(run)+"/Performance/"+filter_name_subdir+"/"+gains+"/"+filter_info_subdir+"/Linearity_OFC_AWG_ADC"+str(adc)+"_Channel"+str(channel)+"_"+gains+".png"))
  print("File Saved: "+plotsDir+short(run)+"/Performance/"+filter_name_subdir+"/"+gains+"/"+filter_info_subdir+"/Linearity_OFC_AWG_ADC"+str(adc)+"_Channel"+str(channel)+"_"+gains+".png")
  plt.clf()
  plt.cla()
  plt.close()
  reset_color()

  #Nonlinearity
  #E OFC vs E AWG, Deviation from linear fit 
  ax = plt.subplot(111)
  plt.grid(linestyle='--')
  plt.minorticks_on()
  for points in data:
    diffAmps = np.empty(len(awg_amps))
    for i in range(len(awg_amps)):
      diffAmps[i] = points["energy_mean"][i] - (  ((points["energy_mean"][len(points["energy_mean"])-1]-points["energy_mean"][0])/(awg_amps[len(awg_amps)-1]-awg_amps[0]))*(awg_amps[i]-awg_amps[0])  + points["energy_mean"][0] )#Not skipping first point 
    plt.errorbar(awg_amps, diffAmps, yerr=points["energy_std"], fmt='o-', color=next_color(), linewidth=0.4, markersize = 5, capsize=2, label=points["gain"])
  #make_ylabel(plt, r'$\textbf{E}_{\textbf{OFC}}$, deviation from fit [ADC Counts]')
  #make_xlabel(plt, r'$\textbf{E}_{\textbf{AWG}}$ [A.U.]')
  if(enable_latex):
    make_ylabel(plt, r'$\textbf{E}_{\textbf{OFC}}$, deviation from fit [ADC Counts]')
    make_xlabel(plt, r'$\textbf{E}_{\textbf{AWG}}$ [A.U.]')
  else:
    make_ylabel(plt, r'$E OFC$, deviation from fit [ADC Counts]')
    make_xlabel(plt, r'$E AWG$ [A.U.]')
  ymin, ymax = plt.ylim()
  plt.ylim(ymin - 0.15*(ymax-ymin), ymax+0.5*(ymax-ymin))
  drawLabel(ax, adc, channel, filter_name, filter_info, run, 'tl')
  fmtx=mpl.ticker.StrMethodFormatter("{x:,.2f}")
  fmty=mpl.ticker.StrMethodFormatter("{x:,.0f}")
  ax.yaxis.set_major_formatter(fmty)
  ax.xaxis.set_major_formatter(fmtx)
  ax.tick_params(direction='in', which='both')
  ax.get_xaxis().set_ticks_position('both')
  ax.get_yaxis().set_ticks_position('both')
  locs, labels = plt.xticks()
  plt.legend(loc=1, frameon=False)
  #plt.show(block=True)
  plt.savefig(checkDir(plotsDir+short(run)+"/Performance/"+filter_name_subdir+"/"+gains+"/"+filter_info_subdir+"/Nonlinearity_ADC"+str(adc)+"_Channel"+str(channel)+"_"+gains+".png"))
  print("File Saved: "+plotsDir+short(run)+"/Performance/"+filter_name_subdir+"/"+gains+"/"+filter_info_subdir+"/Nonlinearity_ADC"+str(adc)+"_Channel"+str(channel)+"_"+gains+".png")
  plt.clf()
  plt.cla()
  plt.close()
  reset_color()



  #Nonlinearity
  #E Peak vs E AWG, Deviation from linear fit 
  ax = plt.subplot(111)
  plt.grid(linestyle='--')
  plt.minorticks_on()
  for points in data:
    diffAmps = np.empty(len(awg_amps))
    for i in range(len(awg_amps)):
      diffAmps[i] = points["energy_mean_peak"][i] - (  ((points["energy_mean_peak"][len(points["energy_mean_peak"])-1]-points["energy_mean_peak"][0])/(awg_amps[len(awg_amps)-1]-awg_amps[0]))*(awg_amps[i]-awg_amps[0])  + points["energy_mean_peak"][0] )#Not skipping first point 
    plt.errorbar(awg_amps, diffAmps, yerr=points["energy_std_peak"], fmt='o-', color=next_color(), linewidth=0.4, markersize = 5, capsize=2, label=points["gain"])
  #make_ylabel(plt, r'$\textbf{E}_{\textbf{OFC}}$, deviation from fit [ADC Counts]')
  #make_xlabel(plt, r'$\textbf{E}_{\textbf{AWG}}$ [A.U.]')
  if(enable_latex):
    make_ylabel(plt, r'$\textbf{E}_{\textbf{Peak}}$, deviation from fit [ADC Counts]')
    make_xlabel(plt, r'$\textbf{E}_{\textbf{AWG}}$ [A.U.]')
  else:
    make_ylabel(plt, r'$E Peak$, deviation from fit [ADC Counts]')
    make_xlabel(plt, r'$E AWG$ [A.U.]')
  ymin, ymax = plt.ylim()
  plt.ylim(ymin - 0.15*(ymax-ymin), ymax+0.5*(ymax-ymin))
  drawLabel(ax, adc, channel, "Single Point", "Peak", run, 'tl')
  fmtx=mpl.ticker.StrMethodFormatter("{x:,.2f}")
  fmty=mpl.ticker.StrMethodFormatter("{x:,.0f}")
  ax.yaxis.set_major_formatter(fmty)
  ax.xaxis.set_major_formatter(fmtx)
  ax.tick_params(direction='in', which='both')
  ax.get_xaxis().set_ticks_position('both')
  ax.get_yaxis().set_ticks_position('both')
  locs, labels = plt.xticks()
  plt.legend(loc=1, frameon=False)
  #plt.show(block=True)
  plt.savefig(checkDir(plotsDir+short(run)+"/Performance/"+filter_name_subdir+"/"+gains+"/"+filter_info_subdir+"/Nonlinearity_Peak_ADC"+str(adc)+"_Channel"+str(channel)+"_"+gains+".png"))
  print("File Saved: "+plotsDir+short(run)+"/Performance/"+filter_name_subdir+"/"+gains+"/"+filter_info_subdir+"/Nonlinearity_Peak_ADC"+str(adc)+"_Channel"+str(channel)+"_"+gains+".png")
  plt.clf()
  plt.cla()
  plt.close()
  reset_color()


  #Nonlinearity
  #E OFC vs E AWG, Percent deviation from linear fit 
  ax = plt.subplot(111)
  plt.grid(linestyle='--')
  plt.minorticks_on()
  for points in data:
    diffAmps = np.empty(len(awg_amps))
    for i in range(len(awg_amps)):
      diffAmps[i] = points["energy_mean"][i] - (  ((points["energy_mean"][len(points["energy_mean"])-1]-points["energy_mean"][0])/(awg_amps[len(awg_amps)-1]-awg_amps[0]))*(awg_amps[i]-awg_amps[0])  + points["energy_mean"][0] )#Not skipping first point 
    plt.errorbar(awg_amps, np.divide(diffAmps,points["energy_mean"])*100, yerr=np.divide(points["energy_std"], points["energy_mean"])*100, fmt='o-', color=next_color(), linewidth=0.4, markersize = 5, capsize=2, label=points["gain"])
  #make_ylabel(plt, r'$\textbf{E}_{\textbf{OFC}}$, \% deviation from fit')
  #make_xlabel(plt, r'$\textbf{E}_{\textbf{AWG}}$ [A.U.]')
  if(enable_latex):
    make_ylabel(plt, r'$\textbf{E}_{\textbf{OFC}}$, \% deviation from fit')
    make_xlabel(plt, r'$\textbf{E}_{\textbf{AWG}}$ [A.U.]')
  else:
    make_ylabel(plt, r'$E OFC$, \% deviation from fit')
    make_xlabel(plt, r'$E AWG$ [A.U.]')
  ymin, ymax = plt.ylim()
  plt.ylim(-3.5, 5)
  drawLabel(ax, adc, channel, filter_name, filter_info, run, 'tl')
  fmtx=mpl.ticker.StrMethodFormatter("{x:,.2f}")
  fmty=mpl.ticker.StrMethodFormatter("{x:,.0f}")
  ax.yaxis.set_major_formatter(fmty)
  ax.xaxis.set_major_formatter(fmtx)
  ax.tick_params(direction='in', which='both')
  ax.get_xaxis().set_ticks_position('both')
  ax.get_yaxis().set_ticks_position('both')
  locs, labels = plt.xticks()
  plt.legend(loc=1, frameon=False)
  #plt.show(block=True)
  plt.savefig(checkDir(plotsDir+short(run)+"/Performance/"+filter_name_subdir+"/"+gains+"/"+filter_info_subdir+"/Nonlinearity_Percent_ADC"+str(adc)+"_Channel"+str(channel)+"_"+gains+".png"))
  print("File Saved: "+plotsDir+short(run)+"/Performance/"+filter_name_subdir+"/"+gains+"/"+filter_info_subdir+"/Nonlinearity_Percent_ADC"+str(adc)+"_Channel"+str(channel)+"_"+gains+".png")
  plt.clf()
  plt.cla()
  plt.close()
  reset_color()

