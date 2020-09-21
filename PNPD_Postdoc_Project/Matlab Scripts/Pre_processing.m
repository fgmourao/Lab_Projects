
%% Data pre-processing

% (First extract the data using the "Extracting raw LFPs and Events" script)

% - Define sound and behavior epochs
% - Organize channels according to the electrodes map
% - Estimate the modulated signal from digital pulses recorded 
% - Concatenate the modulated signal as channel 1
% - Organize data by trials and by behavior events
% - Filter the data<br />
%   . Relies on the following functions: filter_mod.m
%                                        fun_myfilters.m


%                               - CHANNELS MAP - 

% CS modulating signal
% .Row 1

% mPFC 
% .Row 2,3 -> pre limbic
% .Row 4,5 -> infra limbic

% Hippocampus
% .Row 6   -> CA1
% .Row 7   -> MOL layer
% .Row 8,9 -> GD

% Amygdala
% .Row 10,11 -> lateral
% .Row 12,13 -> basolateral

% Inferior colliculus
% .Row 14,15,16,17 -> Dorsol -> ventral, respectively

% by Flavio Mourao. Nucleo de Neurociencias - NNC.
% email: mourao.fg@gmail.com
% Universidade Federal de Minas Gerais. 
% Started in:  02/2020
% Last update: 09/2020

%% First extract the data with : Extracting_LFPs_and_events.m

%% Organizing channels according to the channels map

% obs.
% During my experiment I set channels 17 to 32 to record but here I corrected the numbers to 1 to 16

% Channels order - correted channels map
% mPFC -> ch/lines: 1 2 3 4
ch_order1 = [20 18 21 19] - 16;

% Hippocampus -> ch/lines: 5 6 7 8
ch_order2 = [22 17 23 24] - 16; 

% Amygdala -> ch/lines: 9 10 11 12
ch_order3 = [25 26 27 28] - 16;

% Inferior colliculus -> ch/lines: 13 14 15 16
ch_order4 = [29 30 31 32] - 16;

% Rearrange  channels
data.raw{1,1}   = data.raw{1,1}([ch_order1 ch_order2 ch_order3 ch_order4],:);       % Raw data - All channels  
data.data{1,1}  = data.data{1,1}([ch_order1 ch_order2 ch_order3 ch_order4],:); % Downsampled data - All channels 

% Clear trash
clear ('ch_order1','ch_order2','ch_order3','ch_order4'); 

%% Define sound epochs. Index and time (sec.)

% lines       : trials
% odd columns : trial start
% even columns: trial stop

% data.events.ts_sort{1,1} --> Video Frames
% Video frames were counted by Bonsai software then I send a sequence of numbers were sent to an Arduino via serial port.
% So,Arduino created packages at every 10 frames and each package sent a digital signal that was recorded by openephys.
data.events.ts_vid = data.events.ts_sort{1,1}; % in seconds 

% data.events.ts_sort{1,2} --> CS modulating frequency
% Timestamps locked to the peaks and valleys of the CS modulating frequency 
data.events.ts_mod = data.events.ts_sort{1,2}; % in seconds

% Time in seconds
data.events.idx_t(:,1) = data.events.ts_mod(abs(diff([0; data.events.ts_mod]))>1);
data.events.idx_t(:,2) = data.events.ts_mod(abs(diff([data.events.ts_mod; 0]))>1);

% Index from decimate data
data.events.idx(:,1) = dsearchn(data.timev,data.events.idx_t(:,1));
data.events.idx(:,2) = dsearchn(data.timev,data.events.idx_t(:,2));

% Index from raw data
data.events.raw_idx(:,1) = dsearchn(data.timev_raw,data.events.idx_t(:,1));
data.events.raw_idx(:,2) = dsearchn(data.timev_raw,data.events.idx_t(:,2));

%% Estimating the CS modulating signal from decimate data

% The timestamps locked to the peaks and valleys of the CS modulating frequency
% were generated in parallel through one of the Arduino Digital I/O pins and then recorded
% by a digital input port of the Open-ephys. Through linear interpolation, these
% time values were used to obtain an instantaneous phase time series, which in turn were
% used to reconstruct/estimate the CS modulating signal itself. In this sense the time-
% frequency analysis could keep engaged with the stimulus presentation.

% Define index for peaks and valleys (these two lines take time ...)
picosTSidxs = dsearchn(data.timev,data.events.ts_mod(1:2:end));
valesTSidxs = dsearchn(data.timev,data.events.ts_mod(2:2:end));

% Interpolated instantaneous phase time series
phiRec = nan(length(data.timev),1);      % initialize the vector
phiRec(picosTSidxs) = 0;                 % peak phase   = 0
phiRec(valesTSidxs) = pi;                % valley phase = pi

xphi = find(~isnan(phiRec)); % Sample points
yphi = phiRec(xphi);         % Sample values    
qp   = 1:length(phiRec);     % Query points

% 1-D data linear interpolation
yphi = interp1(xphi,yphi,(qp)');

% Extrated values from Euler representation of angles
yrec = real(exp(1i.*yphi));

% Reconstructed signal with modulator
data.mod = zeros(length(data.timev),1);  

for ii = 1:size(data.events.idx,1)
    data.mod(data.events.idx(ii,1):data.events.idx(ii,2),1) = yrec(data.events.idx(ii,1):data.events.idx(ii,2),1);
end

clear ('picosTSidxs', 'valesTSidxs', 'phiRec', 'xphi', 'yphi', 'qp', 'yrec','ii')                                                 

% Concatenate the modulator signal as channel 1 in the decimated data variable
data.data{1,1} = [data.mod';data.data{1, 1}];

%% Filter Decimate Data

% Filtering by hand (especially the modulated envelope): 'filter_mod.m'
% and/or
% by the Johnzinho's function: 'fun_myfilters.m'
% Set of filters by VRCarva (https://github.com/vrcarva)


% The 53.71 modulated frequency it will always be positioned 
% in "data.data" cell column 2. The other filters will follow 
% the subsequent columns


% Define frequencies cutoff

% filter_mod function
parameters.filter.modulator       = [51.71 55.71]; % 2

% fun_myfilters function
parameters.filter.deltacutoff     = [1 3];         % 3
parameters.filter.thetacutoff1    = [4 7];         % 4
parameters.filter.thetacutoff2    = [8 12];        % 5
parameters.filter.alphacutoff     = [13 15];       % 6
parameters.filter.betacutoff      = [16 31];       % 7
parameters.filter.lowgammacutoff  = [30 50];       % 8
parameters.filter.highgammacutoff = [62 100];      % 9
parameters.filter.extracutoff1    = [80 140];     % 10
%parameters.filter.extracutoff2    = [1 100];       % 11
%parameters.filter.extracutoff3    = [300 3000];    % 12

% Each cell --> columns: parameters.filters according to the above order 

for jj = 1:size(data.data{1,1},1)
    
    % filter_mod function
    data.data{1,2}(jj,:) = Filter_mod(data.data{1,1}(jj,:),parameters.filter.modulator,parameters.srate);
    
    % fun_myfilters function
    data.data{1,3}(jj,:)  = fun_myfilters(data.data{1,1}(jj,:),parameters.srate,parameters.filter.deltacutoff,'iir','0');
    data.data{1,4}(jj,:)  = fun_myfilters(data.data{1,1}(jj,:),parameters.srate,parameters.filter.thetacutoff1,'iir','0');
    data.data{1,5}(jj,:)  = fun_myfilters(data.data{1,1}(jj,:),parameters.srate,parameters.filter.thetacutoff2,'iir','0');    
    data.data{1,6}(jj,:)  = fun_myfilters(data.data{1,1}(jj,:),parameters.srate,parameters.filter.alphacutoff,'iir','0');
    data.data{1,7}(jj,:)  = fun_myfilters(data.data{1,1}(jj,:),parameters.srate,parameters.filter.betacutoff,'iir','0');
    data.data{1,8}(jj,:)  = fun_myfilters(data.data{1,1}(jj,:),parameters.srate,parameters.filter.lowgammacutoff,'iir','0');
    data.data{1,9}(jj,:)  = fun_myfilters(data.data{1,1}(jj,:),parameters.srate,parameters.filter.highgammacutoff,'iir','0');
    data.data{1,10}(jj,:) = fun_myfilters(data.data{1,1}(jj,:),parameters.srate,parameters.filter.extracutoff1,'iir','0');
%     data.data{1,11}(jj,:) = fun_myfilters(data.data{1,1}(jj,:),parameters.srate,parameters.filter.extracutoff2,'iir','0');
%     data.data{1,12}(jj,:) = fun_myfilters(data.data{1,1}(jj,:),parameters.srate,parameters.filter.modulator,'iir','0');   

end

%% Organizing Trial Data - Considering only behavior periods

% Load idx files (*.mat. behavior events. analyzed in the Guide_Video_Track)
[FileLoaded,Path] = uigetfile({'*.mat'}); % Define file type *.*
data.events.behavior = load(fullfile(Path,FileLoaded));

% Minimum behavior period to be considered
MinTime = 5; % (seconds)

toRemove = find((data.events.behavior.TS_LFPsec(:,2) - data.events.behavior.TS_LFPsec(:,1))<MinTime);

data.events.behavior.TS_LFPindex(toRemove,:) = [];
data.events.behavior.TS_LFPsec(toRemove,:)   = [];

data.events.behavior.TSframes(toRemove,:)    = [];
data.events.behavior.TSseconds(toRemove,:)   = [];

% correcting TS_LFPsec  to zero if the record has started after a viewing time
data.events.behavior.TS_LFPsec = data.events.behavior.TS_LFPsec - min(data.timev_raw);

clear('FileLoaded', 'MinTime','Path','toRemove')

%%
% Correcting timestamps according to the downsamplig
% Check if it is necessary !!!
% The sampling frequency during the video analysis must be the same as the current analysis

data.events.behavior.TS_LFPindex = round(data.events.behavior.TS_LFPindex./parameters.downsampling);

%% Organizing data - Considering only the behavioral events

% Set the parameters for each behavior epoch
% Trial period in seconds. Considering the time of the biggest event
parameters.behavior.trialperiod = ceil(max(data.events.behavior.TS_LFPsec(:,2)-data.events.behavior.TS_LFPsec(:,1))); 
% Pre behavior
parameters.behavior.Tpre        = 10; % (seconds)
% Pos behavior
parameters.behavior.Tpos        = 5; % (seconds)
% Number of behavior events
parameters.behavior.NTrials     = length(data.events.behavior.TS_LFPsec); 


% Cell columns --> frequencies cutoff according filter data
% in each cell
% - rows        - > channels
% - columns     - > time
% - 3 dimension - > behavioral events

% Decimate data
% Initializing behavior periods
data.data_behavior    = cell(parameters.behavior.NTrials,length(data.data));

% Cutting time windows...
for ii = 1:size(data.data,2)
    for jj = 1:size(data.events.behavior.TS_LFPindex,1)       
         data.data_behavior{jj,ii} = data.data{1,ii}(:,data.events.behavior.TS_LFPsec(jj,1) * parameters.srate - parameters.behavior.Tpre * parameters.srate : data.events.behavior.TS_LFPsec(jj,2) * parameters.srate + parameters.behavior.Tpos * parameters.srate);      
    end   
end

% Raw data
% Initializing behavior periods
data.raw_behavior     = cell(parameters.behavior.NTrials,length(data.raw));

% Cutting time windows...
for ii = 1:size(data.raw,2)
    for jj = 1:size(data.events.behavior.TS_LFPindex,1)       
         data.raw_behavior{jj,ii}  = data.raw{1,ii}(:,data.events.behavior.TS_LFPsec(jj,1) * parameters.header.sampleRate - parameters.behavior.Tpre * parameters.header.sampleRate : data.events.behavior.TS_LFPsec(jj,2) * parameters.header.sampleRate + parameters.behavior.Tpos * parameters.header.sampleRate);  
    
    end   
end

% Normalizing total samples with "not a numbers (nan)" in each trial 
% to the exact time period according to the sample rate.

% Decimate data
totalsamples1  = max(max(cellfun(@length,data.data_behavior)));
for ii = 1:size(data.data_behavior,2)
    for jj = 1:size(data.data_behavior,1)        
        if length(data.data_behavior{jj,ii}) < totalsamples1
           data.data_behavior{jj,ii}(:,end:totalsamples1) = nan;
        else
           continue
        end
    end
end

% Raw data
totalsamples2  = max(max(cellfun(@length,data.raw_behavior)));

for ii = 1:size(data.raw_behavior,2)
    for jj = 1:size(data.raw_behavior,1)        
        if length(data.raw_behavior{jj,ii}) < totalsamples2
           data.raw_behavior{jj,ii}(:,end:totalsamples2) = nan;
        else
           continue
        end
    end
end

clear ('totalsamples1','totalsamples2')

% Concatenate trials in 3 dimentions

% Decimate data
for ii = 1:size(data.data_behavior,2)
    data.data_behavior{1,ii} = cat(3,data.data_behavior{:,ii});
end

data.data_behavior(2:end,:) = [];

% Raw Data
for ii = 1:size(data.raw_behavior,2)
    data.raw_behavior{1,ii} = cat(3,data.raw_behavior{:,ii});
end

data.raw_behavior(2:end,:) = [];

% Setting Time

% Initializing time trial vectors
data.time_behavior = linspace(-parameters.behavior.Tpre,parameters.behavior.trialperiod + parameters.behavior.Tpos,length(data.data_behavior{1,1}));
data.time_raw_behavior = linspace(-parameters.behavior.Tpre,parameters.behavior.trialperiod + parameters.behavior.Tpos,length(data.raw_behavior{1,1}));


clear ('totalsamples','ii', 'jj', 'FileLoaded', 'Path', 'MinTime', 'toRemove')

%% Organizing Trial Data - Considering the entire trial period

% Set the parameters for each epoch

parameters.trialperiod = 30; % trial period in seconds
parameters.Tpre        = 30; % pre trial in seconds
parameters.Tpos        = 10; % pos trial in seconds
parameters.NTrials     =  5; % number of trials


% Cell columns --> frequencies cutoff according " F_filter.m "
% in each cell
% - rows        - channels
% - columns     - time
% - 3 dimension - trials

% Decimate Data
% Initializing trial periods
data.data_trials     = cell(parameters.NTrials,length(data.data));

% Cutting time windows...
for ii = 1:size(data.data,2)
    for jj = 1:size(data.events.idx,1)       
         data.data_trials{jj,ii} = data.data{1,ii}(:,data.events.idx(jj,1) - parameters.Tpre * parameters.srate : data.events.idx (jj,2) + parameters.Tpos * parameters.srate);  
    end   
end

% Raw Data
% Initializing trial periods
data.data_raw_trials     = cell(parameters.NTrials,length(data.raw));

% Cutting time windows...
for ii = 1:size(data.raw,2)
    for jj = 1:size(data.events.raw_idx,1)       
         data.data_raw_trials{jj,ii} = data.raw{1,ii}(:,data.events.raw_idx(jj,1) - parameters.Tpre * parameters.header.sampleRate : data.events.raw_idx (jj,2) + parameters.Tpos * parameters.header.sampleRate);  
    end   
end

% Normalizing total samples with "not a numbers (nan)" in each trial 
% to the exact time period according to the sample rate.

% Decimate data
totalsamples1  = max(max(cellfun(@length,data.data_trials)));

for ii = 1:size(data.data_trials,2)
    for jj = 1:size(data.data_trials,1)        
        if length(data.data_trials{jj,ii}) < totalsamples1
           data.data_trials{jj,ii}(:,end:totalsamples1) = nan;
        else
           continue
        end
    end
end

clear ('totalsamples1')

% Raw data
totalsamples2  = max(max(cellfun(@length,data.data_raw_trials)));

for ii = 1:size(data.data_raw_trials,2)
    for jj = 1:size(data.data_raw_trials,1)        
        if length(data.data_raw_trials{jj,ii}) < totalsamples2
           data.data_raw_trials{jj,ii}(:,end:totalsamples2) = nan;
        else
           continue
        end
    end
end

clear ('totalsamples2')

% Concatenate trials in 3 dimentions

% Decimate data
for ii = 1:size(data.data,2)
    data.data_trials{1,ii} = cat(3,data.data_trials{:,ii});
end

data.data_trials(2:end,:) = [];

% Raw data
for ii = 1:size(data.raw,2)
    data.data_raw_trials{1,ii} = cat(3,data.data_raw_trials{:,ii});
end

data.data_raw_trials(2:end,:) = [];

% Setting Time

% Initializing time trial vectors
data.time_trials = linspace(-parameters.Tpre,parameters.trialperiod+parameters.Tpos,length(data.data_trials{1,1}));
data.time_raw_trials = linspace(-parameters.Tpre,parameters.trialperiod+parameters.Tpos,length(data.data_raw_trials{1,1}));


clear ('ii', 'jj')

%% last update 20/09/2020 - 01:04am
%  listening: Hope Sandoval - Isn`t it true
