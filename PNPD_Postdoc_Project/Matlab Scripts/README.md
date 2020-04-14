# Matlab Scripts

Project: NEURAL DYNAMICS DURING AUDITORY FEAR CONDITIONING

Main.m<br />
- Call the function/scripts<br />
- Make some plots<br />

Extracting_LFPs_and_Events.m<br />
- Extract and save data from Intan/Open Ephys:  *.continuous and  *.events<br />
- Required function: load_open_ephys_data.m (https://github.com/open-ephys/analysis-tools)<br />

F_filter.m<br />
- It filters the signal by two options: <br />
   _ Filter with parameters defined manually. Matlab buil function: filtfilt.m<br />
   _ Define parameters to 'fun_myfilters.m'.  Filter made by VRCarva (https://github.com/vrcarva) based on EEG_lab: <br />
   
Pre_processing.m<br />
- Define sound epochs<br />
- Organize channels according to the electrodes map<br />
- Estimate the CS modulating signal from digital pulses<br />
- Concatenate the modulator signal as channel 1<br />
- Organize data by trials <br />
- Make some plots<br /> 

sFFT_spectrogram.m
- Short-time FFT by matlab built function spectrogram <br />
- Organizing data trials from the spectrogram<br />  
- Make some plots<br /> 

sFFT_stats.m
- Performs descriptive analysis from the spectrogram <br />
- Make some plots<br /> 

CorCov.m
- Performs Correlation and Covariance Matrices betwwen channels <br />
- Make some plots<br /> 

Hilbert_phase.m
- Phase analyses based on Hilbert Transform. Phase Coherence / Phase lock value (PLV) <br />
- Required circular-statistics-toolbox<br />
 (https://www.mathworks.com/matlabcentral/fileexchange/10676-circular-statistics-toolbox-directional-statistics)<br />
- Make some plots<br /> 

Track_3d.m<br /> 
- Plot electrophysiological parameters as a heat map on the animal's track<br /> 

Video_Analyses.m<br /> 
- Behavioral Analysis from *.csv files imported from Bonsai software and beginning and end of trials index identified through the Video Guide<br /> 

    Outputs:

    . x and y coordinates in pixels<br /> 
    . x and y coordinates in cm<br /> 
    . x / y derivative in cm<br /> 
    . displacement in cm <br /> 
    . accumulated distance over time in cm <br /> 
    . total distance covered in cm <br /> 
    . time vector<br />  


- Make some plots<br /> 

All code by Flavio Mourao. Nucleo de Neurociencias - NNC<br />
email: mourao.fg@gmail.com<br />
Universidade Federal de Minas Gerais<br />
Brazil<br />
