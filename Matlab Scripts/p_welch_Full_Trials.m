%%  Welch power spectral density estimate 

% Flavio Mourao. Nucleo de Neurociencias NNC.
% email: mourao.fg@gmail.com
% Universidade Federal de pwnas Gerais
% Started in:  07/2019
% Last update: 09/2020

%%
% First extract the data with: Extracting_LFPs_and_events.m
% ... then organize the data with the script: Pre_processing.m

%% Cut Time periods - pre trial and trial events. Ignore "pos" period

% - cell          
% - 1st column    - > before onset
% - 2nd column    - > post onset
% - lines         - > trials

% in each cell
% - rows          - > channels
% - columns       - > time

% Time vector
pw.full_trial.time = (linspace(-parameters.Tpre,parameters.trialperiod+parameters.Tpos,size(data.data_trials{1, 1},2)));

% Time index -->  trial epoch begins / ends

pw.full_trial.time_zero_idx = dsearchn(pw.full_trial.time',0); % time zero index. Behavior Start.

for ii = 1:length(data.events.idx)

    temp  = data.events.idx_t(ii,2) - data.events.idx_t(ii,1); % time end index. trial End.
    pw.full_trial.time_end_idx(ii) = dsearchn(pw.full_trial.time',temp');

end


% from...all trial epochs
% Initialize variables
pw.full_trial.data      = cell(parameters.NTrials,2);

for jj = 1:parameters.NTrials
        
    pw.full_trial.data{jj,1} = data.data_trials{1,1}(:,1:pw.full_trial.time_zero_idx-1,jj);              
    pw.full_trial.data{jj,2} = data.data_trials{1,1}(:,pw.full_trial.time_zero_idx:pw.full_trial.time_end_idx(jj),jj); 

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
pw.full_trial.timewin    = 1000; % in ms

% Convert time window to points
pw.full_trial.timewinpnts  = round(pw.full_trial.timewin/(1000/parameters.srate));

% nFFT
pw.full_trial.nFFT = 2^nextpow2(pw.full_trial.timewinpnts);

% Number of overlap samples
pw.full_trial.overlap = 50;
pw.full_trial.noverlap = floor(pw.full_trial.overlap*0.01 * pw.full_trial.timewinpnts);

for ii = 1:size(pw.full_trial.data,2)
    for jj = 1:size(pw.full_trial.data,1)
        for ll = 1:size(pw.full_trial.data{1,1},1)
            
            if ii == 1 && jj == 1 && ll == 1 
                [pw.full_trial.Pxx{jj,ii}(:,ll),pw.full_trial.freq] = pwelch(pw.full_trial.data{jj,ii}(ll,:),pw.full_trial.nFFT,pw.full_trial.nFFT*.5,pw.full_trial.nFFT,parameters.srate);
            
            else                
                pw.full_trial.Pxx{jj,ii}(:,ll) = pwelch(pw.full_trial.data{jj,ii}(ll,:),pw.full_trial.nFFT,pw.full_trial.overlap,pw.full_trial.nFFT,parameters.srate);
            
            end        
        end
    end
end


clear ('ii','jj','ll')

%% Extract the descriptive stats

% Define frequencies cutoff 

steps = diff(pw.full_trial.freq); % frequency steps according to the fft time window

pw.full_trial.freq2plot{1,1}   = 51.7:steps(1):55.7;    % 1 - modulator
pw.full_trial.freq2plot{1,2}   = 1:steps(1):3;          % 2 - deltacutoff 
pw.full_trial.freq2plot{1,3}   = 4:steps(1):10;         % 3 - thetacutoff1
pw.full_trial.freq2plot{1,4}   = 9:steps(1):12;         % 4 - thetacutoff2
pw.full_trial.freq2plot{1,5}   = 13:steps(1):15;        % 5 - alpha
pw.full_trial.freq2plot{1,6}   = 16:steps(1):31;        % 6 - beta
pw.full_trial.freq2plot{1,7}   = 30:steps(1):50;        % 7 - lowgamma
pw.full_trial.freq2plot{1,8}   = 62:steps(1):100;       % 8 - highgamma
pw.full_trial.freq2plot{1,9}   = 80:steps(1):140;       % 9 - extracutoff1
pw.full_trial.freq2plot{1,10}  = 1:steps(1):100;        % 10 - extracutoff2
%pw.full_trial.freq2plot{1,10} = 300:steps(1):3000;     % 11 - extracutoff3


% Cell columns --> frequencies cutoff 

% Cell Row 1 -> mean for each trial.
%   cell 2  
%   - 1st column    - > before onset
%   - 2nd column    - > post onset
%           in each cell
%           - rows          - > trial events
%           - columns       - > channels

% Cell Row 2 -> mean for each trial normalize from baseline (pre trial onset)
%   in each cell
%   - rows          - > trial events
%   - columns       - > channels

% Cell Row 3 -> total mean session
%   cell 2  
%   - 1st column    - > before onset
%   - 2nd column    - > post onset
%           in each cell
%           - columns       - > channels

% Cell Row 4 -> total mean session normalized from baseline (pre trial onset)
%           in each cell
%           - columns       - > channels

% Cell Row 5 ->  mean over trials full signal
%           - 1st column    - > before onset
%           - 2nd column    - > post onset
%           in each cell
%           - columns       - > channels
%           - rows          - > frequencies

% Cell Row 6 ->  mean over trials full signal decibel normaliztion
%           in each cell
%           - columns       - > channels
%           - rows          - > frequencies

% Initialize
pw.full_trial.stats = cell(5, size(pw.full_trial.freq2plot,2));

% loop over frequencies
for ii = 1:size(pw.full_trial.freq2plot,2)
    for jj = 1:size(pw.full_trial.Pxx,1)
    
    closestfreq = dsearchn(pw.full_trial.freq,pw.full_trial.freq2plot{1,ii}');

    % Mean
    pw.full_trial.stats{1,ii}{1,1}(jj,:) = mean(pw.full_trial.Pxx{jj,1}(closestfreq,:),1);
    pw.full_trial.stats{1,ii}{1,2}(jj,:) = mean(pw.full_trial.Pxx{jj,2}(closestfreq,:),1);

    pw.full_trial.stats{2,ii}(jj,:) = mean(pw.full_trial.Pxx{jj,2}(closestfreq,:),1)./mean(pw.full_trial.Pxx{jj,1}(closestfreq,:),1);
    
    closestfreq = [];

    end
end

for ii = 1:size(pw.full_trial.freq2plot,2)
       
    pw.full_trial.stats{3,ii}{1,1} = mean(pw.full_trial.stats{1, ii}{1, 1},1);
    pw.full_trial.stats{3,ii}{1,2} = mean(pw.full_trial.stats{1, ii}{1, 2},1);
    
    pw.full_trial.stats{4,ii} = mean(pw.full_trial.stats{3, ii}{1, 2},1)./mean(pw.full_trial.stats{3, ii}{1, 1},1);

end

    pw.full_trial.stats{5,1}{1,1} =  mean(cat(3,pw.full_trial.Pxx{:,1}), 3);
    pw.full_trial.stats{5,1}{1,2} =  mean(cat(3,pw.full_trial.Pxx{:,2}), 3);
    
    pw.full_trial.stats{6,1}      =  10*log10(pw.full_trial.stats{5,1}{1,2}./pw.full_trial.stats{5,1}{1,1});

    pw.full_trial.stats{7,1}      =  log(pw.full_trial.stats{5,1}{1,2});

clear('closestfreq','ii','steps')

%% Plot to check full session. Channels per substrate 
%   trial period

% Choose channel
ch = [2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17];

%Define frequencies to plot in each subplot
steps = diff(pw.full_trial.freq); % according to the fft time window

freq2plot = 1:steps(1):100;
closestfreq = dsearchn(pw.full_trial.freq,freq2plot');

figure

for ii = 1:length(ch)
    subplot(4,4,ii)
    %suptitle({'Amplitude Spectrum via short-window FFT';['(window = ' num2str(short_fft.timewin./1000) 's' ' - ' 'overlap = ' num2str(short_fft.overlap) '%)']}) 
    set(gcf,'color','white')

    plot(pw.full_trial.freq(closestfreq),pw.full_trial.stats{5, 1}{1, 2}(closestfreq,ch(ii)),'k', 'linew', 2);
    xlabel('Frequencies (Hz)','FontSize',14), ylabel('Power (mV^2/Hz)','FontSize',14)
    xlim([1 100])
    
    if ii <= 4
       ylim([0 1500])
       
    elseif ii >=5 && ii <=8
       ylim([0 2500]) 
       
    elseif ii >=9 && ii <=12
       ylim([0 2000]) 
    
    else 
       ylim([0 500])
       
    end
end

%% Plot to check full session. Channels per substrate 
%   decibel normalization

% Choose channel
ch = [2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17];

%Define frequencies to plot in each subplot
steps = diff(pw.full_trial.freq); % according to the fft time window

freq2plot = 1:steps(1):100;
closestfreq = dsearchn(pw.full_trial.freq,freq2plot');

figure

for ii = 1:length(ch)
    subplot(4,4,ii)
    %suptitle({'Amplitude Spectrum via short-window FFT';['(window = ' num2str(short_fft.timewin./1000) 's' ' - ' 'overlap = ' num2str(short_fft.overlap) '%)']}) 
    set(gcf,'color','white')

    plot(pw.full_trial.freq(closestfreq),pw.full_trial.stats{6, 1}(closestfreq,ch(ii)),'k', 'linew', 2);
    xlabel('Frequencies (Hz)','FontSize',14), ylabel('Power (mV^2/Hz)','FontSize',14)
    xlim([1 100])
    ylim([0   8])
    if ii <= 4
       ylim([-6 6])
       
    elseif ii >=5 && ii <=8
       ylim([-6 6]) 
       
    elseif ii >=9 && ii <=12
       ylim([-6 6]) 
    
    else 
       ylim([-6 20])
       
    end
end

%% last update 14/09/2020 - 02:35am
%  listening: pj harvey - england