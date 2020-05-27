%% Phase-amplitude Cross-frequency coupling measure
%  Comodulation considering the behavior events

% Based on the Hindiael Belchior & Adriano Tort script
% Instituto do Cerebro - Universidade Federal do Rio Grande do Norte


% The code relies on the following functions:
% --> ModIndex.m - by Adriano Tort, Instituto do Cerebro - Universidade Federal do Rio Grande do Norte

% See Tort et al, 2010 -> 10.1152/jn.00106.2010 
%     Tort et al, 2008 -> 10.1073/pnas.0810524105


% Adapted by:
% Flavio Mourao. Nucleo de Neurociencias NNC.
% email: mourao.fg@gmail.com
% Universidade Federal de Minas Gerais
% Started in:  05/2020
% Last update: 05/2020

%%
% First extract the data with: Extracting_LFPs_and_events.m
% ... then organize the data with the script: Pre_processing.m

%% Run each session sequentially

%% Define the Amplitude and Phase Frequencies

% choose channel
MI_comod.parameters.ch_phase = 12;
MI_comod.parameters.ch_amp   = 16;

MI_comod.lfp_phase = data.data{1,1}(MI_comod.parameters.ch_phase,:);
MI_comod.lfp_amp   = data.data{1,1}(MI_comod.parameters.ch_amp,:);

% Default by author
% MI_comod.parameters.MI_comod.PhaseFreqVector = 2:2:50;
% MI_comod.parameters.MI_comod.AmpFreqVector   = 10:5:200;

% MI_comod.parameters.MI_comod.PhaseFreq_BandWidth = 4;
% MI_comod.parameters.MI_comod.AmpFreq_BandWidth   = 10;

MI_comod.parameters.PhaseFreqVector = 2:0.5:14;
MI_comod.parameters.AmpFreqVector   = 30:0.5:100;

MI_comod.parameters.PhaseFreq_BandWidth = 2;
MI_comod.parameters.AmpFreq_BandWidth   = 12;

%% Do filtering and Hilbert transform on CPU

'CPU filtering';

tic
MI_comod.AmpFreqTransformed   = zeros(length(MI_comod.parameters.AmpFreqVector), length(MI_comod.lfp_amp));
MI_comod.PhaseFreqTransformed = zeros(length(MI_comod.parameters.PhaseFreqVector), length(MI_comod.lfp_phase));

for ii=1:length(MI_comod.parameters.AmpFreqVector)
    Af1 = MI_comod.parameters.AmpFreqVector(ii);
    Af2 = Af1+MI_comod.parameters.AmpFreq_BandWidth;
    AmpFreq = eegfilt(MI_comod.lfp_amp,parameters.srate,Af1,Af2); % just filtering
    MI_comod.AmpFreqTransformed(ii, :) = abs(hilbert(AmpFreq)); % getting the amplitude envelope
end

for jj=1:length(MI_comod.parameters.PhaseFreqVector)
    Pf1 = MI_comod.parameters.PhaseFreqVector(jj);
    Pf2 = Pf1 + MI_comod.parameters.PhaseFreq_BandWidth;
    PhaseFreq = eegfilt(MI_comod.lfp_phase,parameters.srate,Pf1,Pf2); % this is just filtering 
    MI_comod.PhaseFreqTransformed(jj, :) = angle(hilbert(PhaseFreq)); % this is getting the phase time series
end
toc

clear('ans', 'ii', 'jj', 'Af1', 'Af2', 'Pf1', 'Pf2', 'AmpFreq', 'PhaseFreq');

%% Set the parameters for each behavior epoch

% Trial period in seconds. Considering the time of the biggest event
MI_comod.parameters.trialperiod = ceil(max(data.events.behavior.TS_LFPsec(:,2)-data.events.behavior.TS_LFPsec(:,1))); 
% Pre behavior
MI_comod.parameters.Tpre        = 5; % (seconds)
% Pos behavior
MI_comod.parameters.Tpos        = 5; % (seconds)
% Number of behavior events
MI_comod.parameters.NTrials     = length(data.events.behavior.TS_LFPsec); 


% Cell events. Full period - pre/during/pos
% - rowns - > behavioral events

% in each cell
% - rows        - > phase/amplitude
% - columns     - > time


% Initializing behavior periods
MI_comod.data_amp    = cell(MI_comod.parameters.NTrials,1);
MI_comod.data_phase  = cell(MI_comod.parameters.NTrials,1);

% Cutting time windows...
for ii = 1:size(data.data,2)
    for jj = 1:size(data.events.behavior.TS_LFPindex,1)       
         MI_comod.data_amp{jj,ii}   = MI_comod.AmpFreqTransformed(:,data.events.behavior.TS_LFPindex(jj,1) - MI_comod.parameters.Tpre * parameters.srate : data.events.behavior.TS_LFPindex(jj,2) + MI_comod.parameters.Tpos * parameters.srate);  
         MI_comod.data_phase{jj,ii} = MI_comod.PhaseFreqTransformed(:,data.events.behavior.TS_LFPindex(jj,1) - MI_comod.parameters.Tpre * parameters.srate : data.events.behavior.TS_LFPindex(jj,2) + MI_comod.parameters.Tpos * parameters.srate);  
    end   
end

% Normalizing total samples with "not a numbers (nan)" in each trial 
% to the exact time period according to the sample rate.

totalsamples  = max(max(cellfun(@length,MI_comod.data_amp)));

for ii = 1:size(MI_comod.data_amp,2)
    for jj = 1:size(MI_comod.data_amp,1)        
        if length(MI_comod.data_amp{jj,ii}) < totalsamples
           MI_comod.data_amp{jj,ii}(:,end:totalsamples) = nan;
           MI_comod.data_phase{jj,ii}(:,end:totalsamples) = nan;
        else
           continue
        end
    end
end

clear ('totalsamples')


% Setting Time

% Initializing time trial vectors
MI_comod.parameters.data_time = linspace(-MI_comod.parameters.Tpre,MI_comod.parameters.trialperiod + MI_comod.parameters.Tpos,length(MI_comod.data_amp{1,1}));


clear ('totalsamples','ii', 'jj', 'FileLoaded', 'Path', 'MinTime', 'toRemove')

%% Cut Time periods - pre behavior and during behavior events. Ignore "pos" period

% Time index -->  behavior epoch begins / ends

MI_comod.parameters.data_time_zero_idx = dsearchn(MI_comod.parameters.data_time',0); % time zero index. Behavior Start.

for ii = 1:length(data.events.behavior.TS_LFPsec)

    temp  = data.events.behavior.TS_LFPsec(ii,2) - data.events.behavior.TS_LFPsec(ii,1); % time end index.  Behavior End.
    MI_comod.parameters.data_time_end_idx(ii) = dsearchn(MI_comod.parameters.data_time',temp');

end

% from...all behavior epochs over time
% Initialize variables

% Cell events. Before and during
% - rowns - > behavioral events
% - first column -> phase
% - second column -> amplitude

% in each cell
% - rows        - > phase/amplitude
% - columns     - > time

MI_comod.data_before      = cell(MI_comod.parameters.NTrials,2);
MI_comod.data_before(:,1) = {nan(length(MI_comod.parameters.PhaseFreqVector),max(MI_comod.parameters.data_time_zero_idx)-1)};
MI_comod.data_before(:,2) = {nan(length(MI_comod.parameters.AmpFreqVector),max(MI_comod.parameters.data_time_zero_idx)-1)};

MI_comod.data_during      = cell(MI_comod.parameters.NTrials,2);
MI_comod.data_during(:,1) = {nan(length(MI_comod.parameters.PhaseFreqVector),max(MI_comod.parameters.data_time_end_idx - MI_comod.parameters.data_time_zero_idx)+1)};
MI_comod.data_during(:,2) = {nan(length(MI_comod.parameters.AmpFreqVector),max(MI_comod.parameters.data_time_end_idx - MI_comod.parameters.data_time_zero_idx)+1)};


for jj = 1:MI_comod.parameters.NTrials
        
    MI_comod.data_before{jj,1}(:,1:(MI_comod.parameters.data_time_zero_idx-1)) = MI_comod.data_phase{jj,1}(:,1:MI_comod.parameters.data_time_zero_idx-1);              
    MI_comod.data_before{jj,2}(:,1:(MI_comod.parameters.data_time_zero_idx-1)) = MI_comod.data_amp{jj,1}(:,1:MI_comod.parameters.data_time_zero_idx-1);              

    MI_comod.data_during{jj,1}(:,1:(MI_comod.parameters.data_time_end_idx(jj) - MI_comod.parameters.data_time_zero_idx)+1) = MI_comod.data_phase{jj,1}(:,MI_comod.parameters.data_time_zero_idx:MI_comod.parameters.data_time_end_idx(jj)); 
    MI_comod.data_during{jj,2}(:,1:(MI_comod.parameters.data_time_end_idx(jj) - MI_comod.parameters.data_time_zero_idx)+1) = MI_comod.data_amp{jj,1}(:,MI_comod.parameters.data_time_zero_idx:MI_comod.parameters.data_time_end_idx(jj)); 

end

clear( 'ii','jj','temp')
%% For comodulation calculation (only has to be calculated once)

MI_comod.parameters.nbin     = 18;
MI_comod.parameters.position = zeros(1,MI_comod.parameters.nbin); % this variable will get the beginning (not the center) of each phase bin (in rads)
MI_comod.parameters.winsize  = 2*pi/MI_comod.parameters.nbin;

for j=1:MI_comod.parameters.nbin 
    MI_comod.parameters.position(j) = -pi+(j-1)*MI_comod.parameters.winsize; 
end


% Define time to analyse
MI_comod.parameters.time2analise = 5;
time2samples = MI_comod.parameters.time2analise * parameters.srate; % in samples


% Comodulogram before and during behavior events
% - rows            -> phase
% - columns         -> amplitude
% - third dimension -> events

MI_comod.Comodulogram_before = single(zeros(length(MI_comod.parameters.PhaseFreqVector),length(MI_comod.parameters.AmpFreqVector),MI_comod.parameters.NTrials));
MI_comod.Comodulogram_during = single(zeros(length(MI_comod.parameters.PhaseFreqVector),length(MI_comod.parameters.AmpFreqVector),MI_comod.parameters.NTrials));

tic
disp('Comodulation loop')

for nn = 1:MI_comod.parameters.NTrials
    for ii = 1:length(MI_comod.parameters.PhaseFreqVector)
        for jj = 1:length(MI_comod.parameters.AmpFreqVector)        
            
            [MI_before,~] = ModIndex(MI_comod.data_before{nn,1}(ii,:), MI_comod.data_before{nn,2}(jj,:), MI_comod.parameters.position);
            [MI_during,~] = ModIndex(MI_comod.data_during{nn,1}(ii, 1:time2samples), MI_comod.data_during{nn,2}(jj, 1:time2samples), MI_comod.parameters.position);
            
            MI_comod.Comodulogram_before(ii,jj,nn) = MI_before;
            MI_comod.Comodulogram_during(ii,jj,nn) = MI_during;
        
        end
    end
end

toc

clear ('ii','jj','nn','time2samples','MI_before','MI_during')

%% Graph comodulogram
%  Events pars - Before and during behavior

event = 6;

figure
set(gcf,'color','w');
suptitle({['Behavior event: ' num2str(event)];' ';' '})

subplot(1,2,1)
toplot = MI_comod.Comodulogram_before(:,:,event);
contourf(MI_comod.parameters.PhaseFreqVector+MI_comod.parameters.PhaseFreq_BandWidth/2, MI_comod.parameters.AmpFreqVector+MI_comod.parameters.AmpFreq_BandWidth/2, toplot',30,'lines','none')
set(gca,'fontsize',14)
ylabel('Amplitude Frequency (Hz)')
xlabel('Phase Frequency (Hz)')
xlim([4 7])
ylim([45 80])
%caxis([0 4*10^-3])
colorbar

title ('Pre behavior')

subplot(1,2,2)
toplot = MI_comod.Comodulogram_during(:,:,event);
contourf(MI_comod.parameters.PhaseFreqVector+MI_comod.parameters.PhaseFreq_BandWidth/2, MI_comod.parameters.AmpFreqVector+MI_comod.parameters.AmpFreq_BandWidth/2, toplot',30,'lines','none')
set(gca,'fontsize',14)
ylabel('Amplitude Frequency (Hz)')
xlabel('Phase Frequency (Hz)')
xlim([4 7])
ylim([45 80])
%caxis([0 4*10^-3])
colorbar
title ('Behavior')
    
clear('event','toplot')

%% Graph comodulogram
%  All events - Before and during behavior

figure
set(gcf,'color','w');

count = 1;

for ii = 1:MI_comod.parameters.NTrials
    
    subplot(2,6,ii)
    toplot = MI_comod.Comodulogram_before(:,:,count);
    contourf(MI_comod.parameters.PhaseFreqVector+MI_comod.parameters.PhaseFreq_BandWidth/2, MI_comod.parameters.AmpFreqVector+MI_comod.parameters.AmpFreq_BandWidth/2, toplot',30,'lines','none')
    set(gca,'fontsize',14)
    ylabel('Amplitude Frequency (Hz)')
    xlabel('Phase Frequency (Hz)')
    xlim([4 7])
    ylim([45 80])
    
    caxis([0 4*10^-3])
    
    title ({['Pre behavior event: ' num2str(count)];' '})
    
    colorbar
    
    count = count + 1;
    
end

count = 1;

for ii = 7:MI_comod.parameters.NTrials*2
    
    subplot(2,6,ii)
    toplot = MI_comod.Comodulogram_during(:,:,count);
    contourf(MI_comod.parameters.PhaseFreqVector+MI_comod.parameters.PhaseFreq_BandWidth/2, MI_comod.parameters.AmpFreqVector+MI_comod.parameters.AmpFreq_BandWidth/2, toplot',30,'lines','none')
    set(gca,'fontsize',14)
    ylabel('Amplitude Frequency (Hz)')
    xlabel('Phase Frequency (Hz)')
    xlim([4 7])
    ylim([45 80])
    
    caxis([0 4*10^-3])
    
    title ({['Behavior event: ' num2str(count)];' '})
    
    colorbar
    
    count = count + 1;
    
end

clear('toplot','count','ii')

%% last update 26/05/2020 - 01:31am
%  listening: Mogwai - Hugh Dallas
