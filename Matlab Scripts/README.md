# Matlab Scripts / Functions

Project: NEURAL DYNAMICS DURING AUDITORY FEAR CONDITIONING

Main.m<br />
- Call functions/scripts<br />
- Make some plots<br />

Extracting_LFPs_and_Events.m<br />
- Extract and save data from Intan/Open Ephys:  *.continuous and  *.events<br />
- The code relies on the following functions : load_open_ephys_data.m (https://github.com/open-ephys/analysis-tools)<br />

Filter_mod.m<br />
- Built filter to avoid deformations at 53.71 modulated frequency.  It guarantees smooth edges in the CS modulating transitions. <br />

fun_myfilters.m<br />
- Set of filters by VRCarva (https://github.com/vrcarva)<br />
   
Pre_processing.m<br />
- Define sound and behavior epochs<br />
- Organize channels according to the electrodes map<br />
- Estimate the modulated signal from digital pulses recorded <br />
- Concatenate the modulated signal as channel 1<br />
- Organize data by trials and by behavior events<br />
- Filter the data<br />
  . Relies on the following functions: filter_mod.m<br />
                                                          fun_myfilters.m by VRCarva (https://github.com/vrcarva)<br />
- Make some plots<br /> 

sFFT_spectrogram.m
- Short-time FFT by matlab built function spectrogram <br />

sFFT_spectrogram_Full_Trials.m
- Organize and plot data from spectrogram considering the trial periods - sound on/off<br />  
- Make some plots<br /> 

sFFT_spectrogram_behavior.m
- Organize and plot data from spectrogram considering only the behavior events<br />  
- Make some plots<br /> 

sFFT_stats_Full_Trials.m
- Performs descriptive analysis from the spectrogram considering the trial periods - sound on/off<br />
- Make some plots<br /> 

sFFT_stats_behavior.m
- Performs descriptive analysis from the spectrogram considering only the behavior events<br />
- Make some plots<br />

p_welch_Full_Trials.m
- Welch power spectral density estimate by matlab built function pwelch <br />
- Performs descriptive analysis considering the trial periods - sound on/off<br />
- Make some plots<br />

p_welch_behavior.m
- Welch power spectral density estimate by matlab built function pwelch <br />
- Performs descriptive analysis considering only the behavior events<br />
- Make some plots<br />

CorCov.m
- Performs Correlation and Covariance Matrices between channels considering the trial periods - sound on/off<br />
- Performs Power Spearman's correlation betwwen channels considering the trial periods - sound on/off<br />
- Make some plots<br /> 

Hilbert_phase_Full_Trials.m
- Phase analyses based on Hilbert Transform. Phase Coherence / Phase lock value (PLV) <br />
- Performs analysis considering the trial periods - sound on/off<br />
- The code relies on the following package: circular-statistics-toolbox  (https://github.com/circstat/circstat-matlab) <br />
- Make some plots<br /> 

Hilbert_phase_behavior.m
- Phase analyses based on Hilbert Transform. Phase Coherence / Phase lock value (PLV) <br />
- Performs analysis considering only the behavior events<br />
- The code relies on the following package: circular-statistics-toolbox  (https://github.com/circstat/circstat-matlab) <br />
- Make some plots<br /> 

MI_behavior.m
- Phase-amplitude Cross-frequency coupling measure<br />
- Performs analysis  considering only the behavior events<br />
- Performs analysis with raw and surrogate values<br />
- The code relies on the following functions: ModIndex.m and shuffle_esc.m<br />
- Make some plots<br /> 

Comodulation_MI_behavior.m
- Phase-amplitude Cross-frequency coupling measure<br />
- Performs comodulation considering only the behavior events<br />
- Code based on the Hindiael Belchior & Adriano Tort script<br />
  Instituto do Cerebro - Universidade Federal do Rio Grande do Norte<br />
- The code relies on the following functions: ModIndex.m <br />
- Make some plots<br /> 

Granger_behavior.m
- Established measure of directed functional connectivity <br />
- Performs analysis considering only the behavior events<br />
- The code relies on the following package: BSMART: A Matlab/C Toolbox for Analyzing Brain Circuits (https://brain-smart.org/)<br />
- Code based on the Analyzing Neural Time Series Data: Theory and Practice by Mike X Cohen<br />

Track_3d.m<br /> 
- Plot electrophysiological parameters as a heat map on the animal's track<br /> 

Video_Analyses.m<br /> 
- Behavioral analysis from *.csv files imported from Bonsai software (https://open-ephys.org/bonsai) and from frame trials (sound on/off) index identified through the Video Guide<br /> 
- Make some plots<br /> 

All code by Flavio Mourao. Nucleo de Neurociencias - NNC<br />
email: mourao.fg@gmail.com<br />
Universidade Federal de Minas Gerais<br />
Brazil<br />
