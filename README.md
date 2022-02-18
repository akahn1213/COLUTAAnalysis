# COLUTAAnalysis

Framework for analyzing the performance of the COLUTAv2 ADC

download with: git clone -b hdf5_processing https://gitlab.cern.ch/alkahn/COLUTAAnalysis.git

This package is used to analyze the hdf5 outputs of the COLUTAv2 GUI

It is split into a number of independent components. A more in-depth summary of each component can be found after the processing examples:

1. Physics/pulseAnalysis.py - Takes input hdf5 files (raw data), interlaces the pulses, and saves the interlaced output in Physics/Blocks_Final
2. Physics/calcOFCs.py - Uses the blocks in Physics/Blocks_Final as input, and calculates Optimal Filter Coefficients using them, which are saved in Physics/OFCs
3. Physics/ApplyFilter.py - Uses the blocks in Physics/Blocks_Final and the OFCs in Physics/OFCs and applies the Optimal Filter, saving its output in Physics/Filter_Output
4. Physics/PerformancePlots.py - Makes performance plots of the filter output saved in Physics/Filter_Output

Processing Examples
------------------
Processing Raw Pulse Data from COLUTAv2 Testing Board

1. Place output HDF5 files into Data/Raw/Pulses

2. Edit Scripts/convert_hdf5.py
  - Change the input file name
  - Edit the measurement ranges and the amplitudes that they correspond to
  - Edit any pulse parameters, such as samples per pulse, number of pulses, etc in the output attributes
  - Run the script

3. Run Physics/pulseAnalysis.py, giving a run number with the -r option
  - Ex: python pulseAnalysis.py -r Run_1174

4. Investigate the finely sampled pulse-blocks using testhdf5.py, giving a run number
  - Ex: python testhdf5.py -r Run_1174

5. Calculate OFCs using the finely sampled pulses with calcOFCs.py, giving a run number, amplitude, gain, adc number, and channel number
  - Ex: python calcOFCs.py -r Run_1174 -a 1p0 -g 1x -d 1 -c 1

6. Apply the OFCs to the pulse-blocks with ApplyFilter.py, giving a run number. OFCs are determined by giving an amplitude and a gain as well
  - Ex: python ApplyFilter.py -r Run_1174 -a 1p0 -g 1x

7. Plot the filter output using PerformancePlots.py, giving a run number. This script must be run with python3
  - Ex: python3 performancePlots.py -r Run_1174



Processing Pedestal Data from COLUTAv2 Testing Board

1. Rename output HDF5 files as Data/Processed/Pedestal/[Run_Name]/Pedestal_Data.hdf5
2. Run Physics/pedestalAnalysis.py with python3, giving a run number with the -r option
  - Ex: python3 pedestalAnalysis.py -r Run_1206



Summary of Components
---------------------

### 1. pulseAnalysis.py: 
#### -Options: -r <Run_Name> 

This script looks for runs in Data/Processed/Pulses/. If no run, or list of runs, is given with the option -r, the script will run automatically on all runs that it finds. 

A run should contain .hdf5 files of amplitude and gain combinations, such as Pulse_Amp0p5_1x.hdf5

The script uses the file names in the run's directory to see what kinds of amplitudes and gains are available. The following analysis depends on the availability of some of these amplitude/gain combinations. For example, if there is no 1x or 4x data available, the code will not try to calculate a gain constant and will just use 4.0 instead.

The script then stores all of the data it finds into a python dictionary, one for each amplitude, whose keys are the gain types ("1x", "4x", or "AG") and whose values are a hyper-rectangular array of shape (#ADCs, #Channels, #Blocks, #Samples_per_block). One block is one interlaced pulse. This dictionary is then analyzed all at once. 

Analysis loop (One per amplitude):
  - Interlace the data by taking each 32nd sample in the original data (in the case of pulses that are 32 samples long), and saving them into a new array.
  - Find the 1x, 4x, and AG baseline shifts by taking the average of the first 5 samples of the interlaced blocks for each respective gain type
  - Determine the gain constant. If 1x or 4x data is not available, or if the 4x data saturates, then a gain constant of 4.0 is used instead. Otherwise, the gain constant is determined by comparing the ratio of samples on the negative lobe of the 1x and 4x pulses, after they have been shifted by their respective baselines found in step 2
  - Shift all pulses by the baseline shifts found in step 2, and multiply the 1x samples by the gain constant found in step 3. This gain constant is applied to all samples whose decision bit == 1. 
  - Align the blocks amongst themselves by minimizing the sum of differences squared along the rising edge, and ensure that the peak of the averaged pulse (over all blocks of the same gain type) occurs at sample 165. 

The blocks are then saved as .hdf5 files in Physics/Blocks_Final/Run_XXXX, where the run number is the same as the input run number. The group structure is as follows: file_name["coluta#/channel#/(samples,bits)"] will return a 2-D array with (#Blocks_per_measurement * #Measurements) rows and (#Samples_per_block) columns. Each row in this 2-D array is a finely sampled pulsewhich was found using the steps described previously. One block will be saved for each combination of amplitude and gain, just as the input data was stored. An example output file name would be Physics/Blocks_Final/Run_XXXX/Block_Amp0p5_1x.hdf5. 


### 2. calcOFCs.py: 
#### -Options -r <Run_Name> -a <Amplitude [e.g. 0p75]> -g <Gain [e.g. 1x]> 

This script uses the output of pulseAnalysis.py, stored in Physics/Blocks_Final/, to calculate a set of Optimal Filter Coefficients (OFCs)

The options -a and -g are used to choose an amplitude and gain to use to calculate OFCs. Given these inputs, along with a Run Number to search for blocks .hdf5 files in, the calculation proceeds as follows:

- Generate an average pulse, by taking the mean of the pulses in the 2-D blocks array according to the run, amplitude, gain, adc, and channel given. Call this averaged interlaced pulse "g"
- Calculate its derivative, "dg"
- Within a predefined range of samples (from 2 to 9 by default), perform the following loop. Note: For ease of discussion, assume each 25ns sample is 30 samples apart in the finely sampled pulse. This corresponds to the "default" configuration, with the ADC clock at 40MHz, and the AWG clock at 1200MHz

OFC Loop:

Loop through #of samples (from 2 to 9)
  - Loop through # of phases (30)
    - Save the 25ns samples and derivatives of which we want to calculate the corresponding OFCs for, with the following method:
      - Start at sample 120 (Remember the peak is at sample 165)
      - Add the current phase # 
      - Add every 30th sample to an array, starting at sample 120+phase#, until #samples samples have been added
        -Ex: Phase 15, 5 samples will save the following samples of g and dg: [135, 165, 195, 225, 255]. Note that the peak, sample 165, is included in phase 15

  - Using these new arrays of g and dg samples, the amplitude of the pulse (max value of g), and the autocorrelation function, calculate the a and b OFCs using the method described in the paper by Cleland & Stern: "Signal processing considerations for liquid ionization calorimeters in high rate environment"
  4- Save the OFCs as a new set of .hdf5 files in Physics/OFCs/[Run]/[Amp]/[Gain]/. The OFC .hdf5 file is saved under the name of COLUTA_OFCs.hdf5 and has a group structure of file["/coluta#/channel#/(nSamples)S/(a,b)"]. For example, file["/coluta1/channel1/5S/a"] will return a 2-D array of the "a-coefficient" OFCs, with 30 rows, corresponding to each phase starting at 0, and 5 columns, corresponding to the 5 samples, for COLUTA 1 Channel 1. OFCs for all channels + ADCs are created with this script 


### 3. ApplyFilter.py: 
#### -Options -r <Blocks Run_Name> -o <OFCs Run_Name> -a <OFCs Amplitude [e.g. 0p75]> -g <OFCs Gain [e.g. 1x]>

  This script uses the blocks in Physics/Blocks_Final/[Blocks Run_Name]/ and the Optimal Filter Coefficients (OFCs) in Physics/OFCs/[OFCs Run_Name]/[OFCs Amplitude]/[OFCs Gain]/COLUTA_OFCs.hdf5 to apply the Optimal Filter, and save the resulting Energy and Time calculations in an output .hdf5 file

  If the -o option is not specified, the OFCs run name will be set to the Blocks run name automatically.

  Using the same loop method as in calcOFCs.py, this script determines E (in ADC counts) and t (in ns) by multiplying the a and b coefficients to the respective pulse samples in the blocks .hdf5 file. A value of E and t is determined for each block, i.e. each row in the blocks file's 2-D array of blocks.

  The resulting E and t measurements are saved in an output .hdf5 file, stored in Physics/Filter_Output/OF/. Two files will be saved here, one named OF_Output.hdf5, which uses the OFCs to calculate E and t, and OF_Output_Peak, which uses the pulse's peak to calculate E, and the OFCs to calculate t. 

  The output .hdf5 files are saved with a group structure of file["/coluta#/channel#/(gain)/Amp(amp)/(#Samples)S/P(Phase#)/(energy,timing)"]. For example, file["/coluta1/channel1/1x/Amp0p75/5S/P15/energy"] will return a 1-D array of energies, with one element for each block in the original blocks .hdf5 file, where each energy value was calculated using OFCs with 5 samples, at phase 15, and the pulses are taken using coluta 1 channel 1, at 1x gain, with a pulse of amplitude 0.75.

### 4. PerformancePlots.py: 
#### -Options -r <Run_Name> 

  This script uses the filter output from ApplyFilter.py to produce performance plots of the OFCs.
  
  All of the following plots are made using all of the data from one run, with as much data is available:

  1. Energy Resolution: Sigma E, using OFCs, vs Amplitude
  2. Energy Resolution: Sigma E, using the pulse peak, vs Amplitude
  3. Energy Reconstruction: E from OFCs vs E from the pulse peak
  4. Energy Resolution Ratio: Sigma E using OFCs / E using OFCs, vs Amplitude
  5. Linearity: E from OFCs, vs input AWG Amplitude
  6. Nonlinearity: E from OFCs, vs input AWG Amplitude, deviation in counts from a linear fit between the lowest and highest amplitudes
  7. Nonlinearity: E from OFCs, vs input AWG Amplitude, percent deviation from a linear fit between the lowest and highest amplitudes 
  8. Timing Resolution: Sigma t using OFCs, vs Amplitude
  9. Timing Correlation: 2D scatter plot of timing in one phase vs timing in the next phase (i.e. phase 15 vs phase 16)
  10. Timing Histogram: Distribution of timing for one phase (i.e. 15)
  11. Timing Histogram: Distribution of timing for the following phase (i.e. 16)
  12. Timing Difference: Distribution of the difference in timing between two adjacent phases (i.e. phase 16 - phase 15)
  13. Timing Sum: Distribution of the sum in timing between two adjacent phases (i.e. phase 15 + phase 16)
  


  These plots are defined in the plotting script atlas_plotter.py. This plotting script uses matplotlib to emulate the standard ATLAS plotting style without the need of ROOT. 


  Outputs are saved in Plotting/PlotsSaved/Pulses/OFCs/[Run_Name]

  TO RUN THIS SCRIPT WITHOUT LATEX: Edit the atlas_plotter.py script and make sure enable_latex is set to False.
