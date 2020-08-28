
%%  Organize and plot data from spectrogram 
%   Considering the behavior events 

% Flavio Mourao. Nucleo de Neurociencias NNC.
% email: mourao.fg@gmail.com
% Universidade Federal de Minas Gerais
% Started in:  04/2020
% Last update: 04/2020

%%
% First extract the data with: Extracting_LFPs_and_events.m
% ... then organize the data with the script: Pre_processing.m
% and then run: sFFT_spectrogram.m

%% Define indexes from spectrogram - Considering the entire trial period

% rows -> trials
% columns 1 -> time before sound epoch 
% columns 2 -> sound start
% columns 3 -> sound stop
% columns 4 -> time after sound epoch

parameters.Tpre        = 30; % pre event in seconds
parameters.Tpos        = 10; % pos event in seconds

% Sound epochs
short_fft.time_idx_t      = (data.events.idx_t(:));
short_fft.time_idx(:,2:3) = reshape(dsearchn(short_fft.time',short_fft.time_idx_t),5,2);
short_fft.time_idx_t      = reshape(short_fft.time_idx_t,5,2);  % just to keep the same format m x n

% Pre sound
short_fft.time_idx(:,1) = dsearchn(short_fft.time',(short_fft.time_idx_t(:,1) - parameters.Tpre));

% Pos sound
short_fft.time_idx(:,4) = dsearchn(short_fft.time',(short_fft.time_idx_t(:,2) + parameters.Tpos));

%% Define indexes from spectrogram - Considering the behavior periods

% rows -> trials
% columns 1 -> time before behavior epoch 
% columns 2 -> behavior start
% columns 3 -> behavior stop
% columns 4 -> time after behavior epoch

% behavior epochs
short_fft.behavior.time_idx_t = (data.events.behavior.TS_LFPsec(:));
dsearchn(short_fft.time',short_fft.behavior.time_idx_t);

short_fft.behavior.time_idx(:,2:3) = reshape(dsearchn(short_fft.time',short_fft.behavior.time_idx_t),length(data.events.behavior.TS_LFPsec),2);
short_fft.behavior.time_idx_t = reshape(short_fft.behavior.time_idx_t,length(data.events.behavior.TS_LFPsec),2);  % just to keep the same format m x n

% Pre behavior
short_fft.behavior.time_idx(:,1) = dsearchn(short_fft.time',(short_fft.behavior.time_idx_t(:,1) - parameters.behavior.Tpre));

% Pos behavior
short_fft.behavior.time_idx(:,4) = dsearchn(short_fft.time',(short_fft.behavior.time_idx_t(:,2) + parameters.behavior.Tpos));

% Trial period in seconds. Considering the time of the biggest event
%short_fft.behavior.trialperiod = dsearchn(short_fft.time',parameters.behavior.trialperiod);


%% Organizing trials data during behavior period from spectrogram 

% Concatenate trial epochs (pre behavior, behavior, pos behavior) in fourth dimensions
% lines: frequencies / columns: time / third dimension: channels / fourth dimension: trials

short_fft.behavior.data_trials = complex(zeros(length(short_fft.freq),max(short_fft.behavior.time_idx(:,4)-short_fft.behavior.time_idx(:,1)) + 1,parameters.nch+1,parameters.behavior.NTrials));

for jj = 1:parameters.behavior.NTrials
    temp = short_fft.data(:,short_fft.behavior.time_idx(jj,1):short_fft.behavior.time_idx(jj,4),:);
    short_fft.behavior.data_trials(:,1:size(temp,2),:,jj) = temp;
    temp = [];
end

short_fft.behavior.data_trials(short_fft.behavior.data_trials==0) = nan;

% Time vector to plot
short_fft.behavior.time_trials = (linspace(-parameters.behavior.Tpre,parameters.behavior.trialperiod + parameters.behavior.Tpos,size(short_fft.behavior.data_trials,2)));



% Concatenate trial epochs (pre behavior) in fourth dimensions
% lines: frequencies / columns: time / third dimension: channels / fourth dimension: trials

short_fft.behavior.data_trials_pre = complex(zeros(length(short_fft.freq),max(short_fft.behavior.time_idx(:,2)-short_fft.behavior.time_idx(:,1)) - 1, parameters.nch+1,parameters.behavior.NTrials));

for jj = 1:parameters.behavior.NTrials
    temp = short_fft.data(:,short_fft.behavior.time_idx(jj,1):short_fft.behavior.time_idx(jj,2)-1,:);
    short_fft.behavior.data_trials_pre(:,1:size(temp,2),:,jj) = temp;
    temp = [];
end

short_fft.behavior.data_trials_pre(short_fft.behavior.data_trials_pre == 0) = nan;



% Concatenate trial epochs (behavior period) in fourth dimensions
% lines: frequencies / columns: time / third dimension: channels / fourth dimension: trials

short_fft.behavior.data_trials_behavior = complex(zeros(length(short_fft.freq),max(short_fft.behavior.time_idx(:,3)-short_fft.behavior.time_idx(:,2)) + 1,parameters.nch+1,parameters.behavior.NTrials));

for jj = 1:parameters.behavior.NTrials
    temp = short_fft.data(:,short_fft.behavior.time_idx(jj,2):short_fft.behavior.time_idx(jj,3),:);
    short_fft.behavior.data_trials_behavior (:,1:size(temp,2),:,jj) = temp;
    temp = [];
end

short_fft.behavior.data_trials_behavior(short_fft.behavior.data_trials_behavior==0) = nan;


% Time vector to plot
short_fft.behavior.time_trials_behavior = (linspace(0,parameters.behavior.trialperiod,size(short_fft.behavior.data_trials_behavior,2)));


clear ('temp','ii','jj')

%% Plot to check - Pre and pos behavior period

% Choose channel
ch = 16;

%Define frequencies to plot in each subplot
steps = diff(short_fft.freq); % according to the fft time window

freq2plot = 40:steps(1):70;
closestfreq = dsearchn(short_fft.freq,freq2plot');

%Define time to plot in each subplot
time2plot1 = [ short_fft.behavior.time_idx(:,1),short_fft.behavior.time_idx(:,2) ]; % pre sound
time2plot2 = [ short_fft.behavior.time_idx(:,2),short_fft.behavior.time_idx(:,3) ]; % sound period

% Spectrogram

% - dashed white lines : sound periods
% - red lines          : behavior start
% - black lines        : behavior stop

figure

subplot(3,parameters.behavior.NTrials,[1 parameters.behavior.NTrials])
suptitle({'Amplitude Spectrum via short-window FFT';['(window = ' num2str(short_fft.timewin./1000) 's' ' - ' 'overlap = ' num2str(short_fft.overlap) '%)'];[];['Channel - ' num2str(ch)]}) 
set(gcf,'color','white')

contourf(short_fft.time,short_fft.freq(closestfreq),abs(short_fft.data(closestfreq,:,ch)),80,'linecolor','none');
xlabel('Time (s)','FontSize',14), ylabel('Frequency (Hz)','FontSize',14)
xlim([short_fft.time(1) 600])
colorbar
%caxis([0 6.5*10^5])
%colorbar('Location','eastoutside','YTick',[]);


% Plot idx`s on spectrogram

z  = max(max(abs(short_fft.data(closestfreq,:,ch))));
zp = repmat (z, 1, length(short_fft.freq(closestfreq)));
    
for ii = 1:length(short_fft.behavior.time_idx_t)
    hold on   
    plot3(short_fft.behavior.time_idx_t(ii,1) + zeros(1,length(zp)),freq2plot, zp,'Color','[0.6350, 0.0780, 0.1840]','linew',2)
    plot3(short_fft.behavior.time_idx_t(ii,2) + zeros(1,length(zp)),freq2plot, zp,'k','linew',2)
end

for ii = 1:length(short_fft.time_idx_t)
    hold on
    plot3(short_fft.time_idx_t(ii,1) + zeros(1,length(zp)),freq2plot, zp,'w--','linew',2)
    plot3(short_fft.time_idx_t(ii,2) + zeros(1,length(zp)),freq2plot, zp,'w--','linew',2)

end


% Plot specs for each epoch

for ii = 1:parameters.behavior.NTrials
    subplot(3,parameters.behavior.NTrials,ii+parameters.behavior.NTrials)
    plot(short_fft.freq, mean(abs(short_fft.data(:,time2plot1(ii,1):time2plot1(ii,2)-1,ch)),2),'k')
    xlim([floor(freq2plot(1)) ceil(freq2plot(end))])
    %ylim([0 15*10^4])
    box off
    title(['event: ',num2str(ii)]) 
end 

for ii = 1:parameters.behavior.NTrials
    subplot(3,parameters.behavior.NTrials,ii+parameters.behavior.NTrials*2)
    plot(short_fft.freq, mean(abs(short_fft.data(:,time2plot2(ii,1):time2plot2(ii,2),ch)),2),'k')
    xlim([floor(freq2plot(1)) ceil(freq2plot(end))])
    %ylim([0 15*10^4])
    xlabel('Frequencies (Hz)')
    ylabel('Energy')
    box off
    title(['Event: ',num2str(ii)]) 
end 

clear ('ch','steps','freq2plot','closestfreq','time2plot1','time2plot2','z','zp', 'ii')

%% Plot to check - change from baseline

% Choose channel
ch = 16;

%Define frequencies to plot in each subplot
steps = diff(short_fft.freq); % according to the fft time window

freq2plot = 40:steps(1):70;
closestfreq = dsearchn(short_fft.freq,freq2plot');

%Define time to plot in each subplot
time2plot1 = [ short_fft.behavior.time_idx(:,1),short_fft.behavior.time_idx(:,2) ]; % pre behavior
time2plot2 = [ short_fft.behavior.time_idx(:,2),short_fft.behavior.time_idx(:,3) ]; % behavior period


% Spectrogram
% - dashed white lines : sound periods
% - red lines          : behavior start
% - black lines        : behavior stop

figure

subplot(2,parameters.behavior.NTrials,[1 parameters.behavior.NTrials])
suptitle({'Amplitude Spectrum via short-window FFT';['(window = ' num2str(short_fft.timewin./1000) 's' ' - ' 'overlap = ' num2str(short_fft.overlap) '%)'];[];['Channel -> ' num2str(ch)]}) 
set(gcf,'color','white')

contourf(short_fft.time,short_fft.freq(closestfreq),abs(short_fft.data(closestfreq,:,ch)),80,'linecolor','none');
xlabel('Time (s)','FontSize',14), ylabel('Frequency (Hz)','FontSize',14)
xlim([short_fft.time(1) 600])
colorbar
%caxis([0 .8*10^5])
%colorbar('Location','eastoutside','YTick',[]);


% Plot idx`s on spectrogram

z  = max(max(abs(short_fft.data(closestfreq,:,ch))));
zp = repmat (z, 1, length(short_fft.freq(closestfreq)));
    
for ii = 1:length(short_fft.behavior.time_idx_t)
    hold on   
    plot3(short_fft.behavior.time_idx_t(ii,1) + zeros(1,length(zp)),freq2plot, zp,'Color','[0.6350, 0.0780, 0.1840]','linew',2)
    plot3(short_fft.behavior.time_idx_t(ii,2) + zeros(1,length(zp)),freq2plot, zp,'k','linew',2)
end

for ii = 1:length(short_fft.time_idx_t)
    hold on
    plot3(short_fft.time_idx_t(ii,1) + zeros(1,length(zp)),freq2plot, zp,'w--','linew',2)
    plot3(short_fft.time_idx_t(ii,2) + zeros(1,length(zp)),freq2plot, zp,'w--','linew',2)

end


% Plot specs for each epoch

for ii = 1:parameters.behavior.NTrials
    subplot(2,parameters.behavior.NTrials,ii+parameters.behavior.NTrials)
    plot(short_fft.freq, mean(abs(short_fft.data(:,time2plot2(ii,1):time2plot2(ii,2),ch))./ mean(abs(short_fft.data(:,time2plot1(ii,1):time2plot1(ii,2)-1,ch)),2),2),'k');
    xlim([floor(freq2plot(1)) ceil(freq2plot(end))])
    ylim([0 8])
    xlabel('Frequencies (Hz)')
    ylabel('Change from baseline')
    box off
    title(['Trial: ',num2str(ii)]) 
end 

clear ('ch','steps','freq2plot','closestfreq','time2plot1','time2plot2', 'z', 'zp', 'ii')

%% Plot to check - decibel normalization

% Choose channel
ch = 11;

%Define frequencies to plot in each subplot
steps = diff(short_fft.freq); % according to the fft time window

freq2plot = 3:steps(1):9;
closestfreq = dsearchn(short_fft.freq,freq2plot');

%Define time to plot in each subplot
time2plot1 = [ short_fft.behavior.time_idx(:,1),short_fft.behavior.time_idx(:,2) ]; % pre behavior
time2plot2 = [ short_fft.behavior.time_idx(:,2),short_fft.behavior.time_idx(:,3) ]; % behavior period


% Spectrogram
% - dashed white lines : sound periods
% - red lines          : behavior start
% - black lines        : behavior stop

figure

subplot(2,parameters.behavior.NTrials,[1 parameters.behavior.NTrials])
suptitle({'Amplitude Spectrum via short-window FFT';['(window = ' num2str(short_fft.timewin./1000) 's' ' - ' 'overlap = ' num2str(short_fft.overlap) '%)'];[];['Channel -> ' num2str(ch)]}) 
set(gcf,'color','white')

contourf(short_fft.time,short_fft.freq(closestfreq),abs(short_fft.data(closestfreq,:,ch)),80,'linecolor','none');
xlabel('Time (s)','FontSize',14), ylabel('Frequency (Hz)','FontSize',14)
xlim([short_fft.time(1) 600])
colorbar
%caxis([0 6.5*10^5])
%colorbar('Location','eastoutside','YTick',[]);


% Plot idx`s on spectrogram

z  = max(max(abs(short_fft.data(closestfreq,:,ch))));
zp = repmat (z, 1, length(short_fft.freq(closestfreq)));
    
for ii = 1:length(short_fft.behavior.time_idx_t)
    hold on   
    plot3(short_fft.behavior.time_idx_t(ii,1) + zeros(1,length(zp)),freq2plot, zp,'Color','[0.6350, 0.0780, 0.1840]','linew',2)
    plot3(short_fft.behavior.time_idx_t(ii,2) + zeros(1,length(zp)),freq2plot, zp,'k','linew',2)
end

for ii = 1:length(short_fft.time_idx_t)
    hold on
    plot3(short_fft.time_idx_t(ii,1) + zeros(1,length(zp)),freq2plot, zp,'w--','linew',2)
    plot3(short_fft.time_idx_t(ii,2) + zeros(1,length(zp)),freq2plot, zp,'w--','linew',2)

end


% Plot specs for each epoch

for ii = 1:parameters.behavior.NTrials
    subplot(2,parameters.behavior.NTrials,ii+parameters.behavior.NTrials)
    plot(short_fft.freq, mean(10*log10(abs(short_fft.data(:,time2plot2(ii,1):time2plot2(ii,2),ch))./ mean(abs(short_fft.data(:,time2plot1(ii,1):time2plot1(ii,2)-1,ch)),2)),2),'k');
    xlim([floor(freq2plot(1)) ceil(freq2plot(end))])
    ylim([-5 15])
    xlabel('Frequencies (Hz)')
    ylabel('Log Power (DB)')
    box off
    title(['Trial: ',num2str(ii)]) 
end 


clear ('ch','steps','freq2plot','closestfreq','time2plot1','time2plot2', 'z', 'zp', 'ii')

%% Plot all trials and mean trials to check

% Choose channel
ch = 16;

% Time to plot
t = [ -5 10];

%Define frequencies to plot
steps = diff(short_fft.freq); % frequency steps according to the fft time window

freq2plot = 50:steps(1):60;
closestfreq = dsearchn(short_fft.freq,freq2plot');

figure
suptitle({'Amplitude Spectrum via short-window FFT';['(window = ' num2str(short_fft.timewin./1000) 's' ' - ' 'overlap = ' num2str(short_fft.overlap) '%) - Channel -> ' num2str(ch)];[];[]}) 
set(gcf,'color','white')

z  = max(max(abs(short_fft.data(closestfreq,:,ch))));
zp = repmat (z, 1, length(short_fft.freq(closestfreq)));

for ii=1:parameters.behavior.NTrials
    subplot(2,round((parameters.behavior.NTrials+1)/2),ii)
    
    contourf(short_fft.behavior.time_trials,short_fft.freq(closestfreq),(abs(short_fft.behavior.data_trials(closestfreq,:,ch, ii))),80,'linecolor','none');
    
    xlabel('Time (s)','FontSize',14), ylabel('Frequency (Hz)','FontSize',14)
    title(['Trial ', num2str(ii)]);
    colormap jet
    colorbar('Location','eastoutside','YTick',[]);
    %caxis([0 1.2*10^5])

    xlim(t)
    
    hold on    
    plot3(0 + zeros(1,length(zp)),freq2plot, zp,'w--','linew',2)

    
end 

subplot (2,round((parameters.behavior.NTrials+1)/2),parameters.behavior.NTrials+1);
contourf(short_fft.behavior.time_trials,short_fft.freq(closestfreq),nanmean(abs(short_fft.behavior.data_trials(closestfreq,:,ch,:)),4),80,'linecolor','none');
xlabel('Time (s)','FontSize',14), ylabel('Frequency (Hz)','FontSize',14)
title('Mean Trials');
colormap jet
colorbar('Location','eastoutside');
%caxis([0 1.0*10^5])
xlim(t)

hold on 
plot3(0 + zeros(1,length(zp)),freq2plot, zp,'w--','linew',2)


clear ('steps','freq2plot','closestfreq','ii','s','jj','ch','t');

%% last update 19/04/2020 - 22:56
%  listening: Caspian - circles in circles
