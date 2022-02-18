============================
====== Pulse Analysis ======
============================


Automatically finds processed data and runs the entire alignment+OFC analysis on it


There are four main scripts in this directory

1. alignSections.py
	- Takes processed pulse data in COLUTA_Analysis/Data/Processed/Pulses and produces a .npy file of "blocks"
	- The output files are 2D objects, where each row is a finely-sampled pulse, and the number of rows is equal 
		to the number of finely sampled pulses that can be made from each amplitude and gain. For one file of 6000 
		samples at 40MHz, we get 5 finely-sampled pulses. The total number of rows for this kind of data is then
		5*nFiles. A blocks file is made for BOTH pulse samples and decision bits

	This is performed on the processed data by:
		1. Finding the 1x, 4x, and AG baselines by finely sampling the 1x, 4x, and AG data without shifting or 
			 scaling any samples, and then determining their vertical offsets. These values are saved as shift1x, 
			 shift4x, and shiftAG
		2. Finding the gain constant by shifting the 1x and 4x finely-sampled pulses by their respective baselines, 
			 and then determining the ratio between the mean of samples 400-600, along the negative lobe, of these
			 pulses
			- If there is no 4x data, or the 4x negative lobes saturate because of an amplitude that is too high,
				the gain constant of 4.0 is chosen instead
		3. Applying both the shift and gain constant, and then saving the blocks

	This script automatically finds the processed data, and only needs to know the versions to run on. Versions
	are the names of directories in COLUTA_Analysis/Data/Processed/Pulses/ in which the processed data is stored.
	The list of versions is written and accessed in the helperfunctions.py script

2. SaveAmplitudes.py
	- Takes blocks files and applies OFCs to the finely sampled pulses. It saves the resulting amplitudes, and 
		amplitude*time values in a many-dimensional array, saved in OFCs/Amplitudes/version/. It also produces
		results for when the pulse and OFCs are out-of-phase, either with the pulse at phase 15 and the OFCs 
		changing (suffix "_offset") or with the OFCs at 15 and the pulse changing (suffix "_offset_pulse")
	- The OFCs used are saved in OFCs/Saved_OFCs/

3. OFCAnalysis.py
	- Takes the amplitudes and timing files made from SaveAmplitudes.py and produces plots of these values. 
		The plots are saved in COLUTA_Analysis/Plottings/PlotsSaved/Pulses/OFCs/version/

4. compare_blocks.py
	- Takes blocks files and makes triple-plots comparing two gains of the same amplitude, i.e. comparing 
		Amp0.5 1x to Amp0.5 AG. 
	- The top plot is overlayed finely-sampled pulses for both gains
	- The middle plot is the difference between the two finely-sampled pulses
	- The bottom plot is the decision bits for the second gain in the comparison. For the example above, this would
		be the decision bits for the AG data
	- Plots are saved in COLUTA_Analysis/Plotting/PlotsSaved/Pulses/Compare/version/

5. calcOFCs.py
	- Using the blocks made from alignSections.py, generates sets of OFCs for a given version and amplitude
	- Run python calcOFCs.py -h for information on usage
	- Uses the identity matrix as an autocorrelation matrix




Notes:
	-	There is currently no code to produce OFCs. I will add this soon. If the need is urgent, let me know and I will
		do so faster.
	- This code automatically finds all of the processed data it needs in order to run. This may break with later 
		sets of data.
	- The bash script runAll.sh will run the first 4 of these scripts in order
		
