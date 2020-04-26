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
- Define sound and behavior epochs<br />
- Organize channels according to the electrodes map<br />
- Estimate the CS modulating signal from digital pulses<br />
- Concatenate the modulator signal as channel 1<br />
- Organize data by trials and by behavior events<br />
- Make some plots<br /> 

sFFT_spectrogram.m
- Short-time FFT by matlab built function spectrogram <br />

sFFT_spectrogram_Full_Trials.m
- Organizing data from the spectrogram considering the trial periods - sound on/off<br />  
- Make some plots<br /> 

sFFT_spectrogram_behavior.m
- Organizing data from the spectrogram considering only the behavior events<br />  
- Make some plots<br /> 

sFFT_stats_Full_Trials.m
- Performs descriptive analysis from the spectrogram considering the trial periods - sound on/off<br />
- Make some plots<br /> 

sFFT_stats_behavior.m
- Performs descriptive analysis from the spectrogram considering only the behavior events<br />
- Make some plots<br />

CorCov.m
- Performs Correlation and Covariance Matrices betwwen channels considering the trial periods - sound on/off<br />
- Make some plots<br /> 

Hilbert_phase_Full_Trials.m
- Phase analyses based on Hilbert Transform. Phase Coherence / Phase lock value (PLV) <br />
- Performs  analysis  considering the trial periods - sound on/off<br />
- Required circular-statistics-toolbox<br />
- Make some plots<br /> 

Hilbert_phase_behavior
- Organizing data from the spectrogram considering only the behavior events<br />  
- Performs  analysis  considering only the behavior events<br />
- Required circular-statistics-toolbox<br />
- Make some plots<br /> 

Track_3d.m<br /> 
- Plot electrophysiological parameters as a heat map on the animal's track<br /> 

Video_Analyses.m<br /> 
- Behavioral Analysis from *.csv files imported from Bonsai software and from frame trials (sound on/off) index identified through the Video Guide<br /> 

    Outputs:

    . x and y coordinates in pixels<br /> 
    . x and y coordinates in cm<br /> 
    . x / y derivative in cm<br /> 
    . Displacement in cm <br /> 
    . Accumulated distance over time in cm <br /> 
    . Total distance covered in cm <br /> 
    . Total distance covered in cm in each trial (sound on) <br />   
    . Total distance covered in all trials (sound on) <br /> 
    . Total distance covered out of trials (sound off) <br /> 
    . time vector<br />  

- Make some plots<br /> 

All code by Flavio Mourao. Nucleo de Neurociencias - NNC<br />
email: mourao.fg@gmail.com<br />
Universidade Federal de Minas Gerais<br />
Brazil<br />
