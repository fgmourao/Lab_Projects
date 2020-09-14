%%  Welch power spectral density estimate 

% Flavio Mourao. Nucleo de Neurociencias NNC.
% email: mourao.fg@gmail.com
% Universidade Federal de pwnas Gerais
% Started in:  07/2019
% Last update: 09/2020

%%
% First extract the data with: Extracting_LFPs_and_events.m
% ... then organize the data with the script: Pre_processing.m

%% Cut Time periods - pre behavior and during behavior events. Ignore "pos" period

% - cell          
% - 1st column    - > before onset
% - 2nd column    - > post onset
% - lines         - > behavioral events

% in each cell
% - rows          - > channels
% - columns       - > time

% Time vector
pw.behavior.time = (linspace(-parameters.behavior.Tpre,parameters.behavior.trialperiod+parameters.behavior.Tpos,size(data.data_behavior{1, 1},2)));

% Time index -->  behavior epoch begins / ends

pw.behavior.time_zero_idx = dsearchn(pw.behavior.time',0); % time zero index. Behavior Start.

for ii = 1:length(data.events.behavior.TS_LFPsec)

    temp  = data.events.behavior.TS_LFPsec(ii,2) - data.events.behavior.TS_LFPsec(ii,1); % time end index.  Behavior End.
    pw.behavior.time_end_idx(ii) = dsearchn(pw.behavior.time',temp');

end


% from...all behavior epochs over time
% Initialize variables
pw.behavior.data      = cell(6,2);

for jj = 1:parameters.behavior.NTrials
        
    pw.behavior.data{jj,1} = data.data_behavior{1,1}(:,1:pw.behavior.time_zero_idx-1,jj);              
    pw.behavior.data{jj,2} = data.data_behavior{1,1}(:,pw.behavior.time_zero_idx:pw.behavior.time_end_idx(jj),jj); 

end


clear( 'ii','jj','temp')

%% pwelch
%  [pxx,f] = pwelch(x,window,noverlap,f,fs)

% - cell          
% - 1st column    - > before onset
% - 2nd column    - > post onset
% - lines         - > behavioral events

% in each cell
% - rows          - > time 
% - columns       - > channels

% Time window
pw.behavior.timewin    = 1000; % in ms

% Convert time window to points
pw.behavior.timewinpnts  = round(pw.behavior.timewin/(1000/parameters.srate));

% nFFT
pw.behavior.nFFT = 2^nextpow2(pw.behavior.timewinpnts);

% Number of overlap samples
pw.behavior.overlap = 95;
pw.behavior.noverlap = floor(pw.behavior.overlap*0.01 * pw.behavior.timewinpnts);

for ii = 1:size(pw.behavior.data,2)
    for jj = 1:size(pw.behavior.data,1)
        for ll = 1:size(pw.behavior.data{1,1},1)
            
            if ii == 1 && jj == 1 && ll == 1 
                [pw.behavior.Pxx{jj,ii}(:,ll),pw.behavior.freq] = pwelch(pw.behavior.data{jj,ii}(ll,:),pw.behavior.nFFT,pw.behavior.nFFT*.5,pw.behavior.nFFT,parameters.srate);
            
            else                
                pw.behavior.Pxx{jj,ii}(:,ll) = pwelch(pw.behavior.data{jj,ii}(ll,:),pw.behavior.nFFT,pw.behavior.overlap,pw.behavior.nFFT,parameters.srate);
            
            end        
        end
    end
end


clear ('ii','jj','ll')

%% Extract the descriptive stats

% Define frequencies cutoff 

steps = diff(pw.behavior.freq); % frequency steps according to the fft time window

pw.behavior.freq2plot{1,1}   = 51.7:steps(1):55.7;    % 1 - modulator
pw.behavior.freq2plot{1,2}   = 1:steps(1):3;          % 2 - deltacutoff 
pw.behavior.freq2plot{1,3}   = 4:steps(1):10;         % 3 - thetacutoff1
pw.behavior.freq2plot{1,4}   = 9:steps(1):12;         % 4 - thetacutoff2
pw.behavior.freq2plot{1,5}   = 13:steps(1):15;        % 5 - alpha
pw.behavior.freq2plot{1,6}   = 16:steps(1):31;        % 6 - beta
pw.behavior.freq2plot{1,7}   = 30:steps(1):50;        % 7 - lowgamma
pw.behavior.freq2plot{1,8}   = 62:steps(1):100;       % 8 - highgamma
pw.behavior.freq2plot{1,9}   = 80:steps(1):140;       % 9 - extracutoff1
pw.behavior.freq2plot{1,10}  = 1:steps(1):100;        % 10 - extracutoff2
%pw.behavior.freq2plot{1,10} = 300:steps(1):3000;     % 11 - extracutoff3


% Cell columns --> frequencies cutoff 

% Cell Row 1 -> mean for each trial (only behavior period).
%   cell 2  
%   - 1st column    - > before onset
%   - 2nd column    - > post onset
%           in each cell
%           - rows          - > behavirol events
%           - columns       - > channels

% Cell Row 2 -> mean for each trial (only behavior period) normalize from baseline (pre behavior onset)
%   in each cell
%   - rows          - > behavirol events
%   - columns       - > channels

% Cell Row 3 -> total mean session (only behavior period)
%   cell 2  
%   - 1st column    - > before onset
%   - 2nd column    - > post onset
%           in each cell
%           - columns       - > channels

% Cell Row 4 -> total mean session (only behavior period) normalize from baseline (pre behavior onset)
%           in each cell
%           - columns       - > channels

% Initialize
pw.behavior.stats = cell(4, size(pw.behavior.freq2plot,2));

% loop over frequencies
for ii = 1:size(pw.behavior.freq2plot,2)
    for jj = 1:size(pw.behavior.Pxx,1)
    
    closestfreq = dsearchn(pw.freq,pw.behavior.freq2plot{1,ii}');

    % Mean
    pw.behavior.stats{1,ii}{1,1}(jj,:) = mean(pw.behavior.Pxx{jj,1}(closestfreq,:),1);
    pw.behavior.stats{1,ii}{1,2}(jj,:) = mean(pw.behavior.Pxx{jj,2}(closestfreq,:),1);

    pw.behavior.stats{2,ii}(jj,:) = mean(pw.behavior.Pxx{jj,2}(closestfreq,:),1)./mean(pw.behavior.Pxx{jj,1}(closestfreq,:),1);
    
    closestfreq = [];

    end
end

for ii = 1:size(pw.behavior.freq2plot,2)
       
    pw.behavior.stats{3,ii}{1,1} = mean(pw.behavior.stats{1, ii}{1, 1},1);
    pw.behavior.stats{3,ii}{1,2} = mean(pw.behavior.stats{1, ii}{1, 2},1);
    
    pw.behavior.stats{4,ii} = mean(pw.behavior.stats{3, ii}{1, 2},1)./mean(pw.behavior.stats{3, ii}{1, 1},1);

end

clear('closestfreq','ii','steps')

%% last update 14/09/2020 - 02:35am
%  listening: pj harvey - england