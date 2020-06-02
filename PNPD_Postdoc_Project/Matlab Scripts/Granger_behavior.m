%% Granger prediction
%  - Performs analysis considering only the behavior events

% The code relies on the following package:
% --> BSMART: A Matlab/C Toolbox for Analyzing Brain Circuits
% https://brain-smart.org/ 

% Code based on: 
% Analyzing Neural Time Series Data: Theory and Practice
% by Mike X Cohen


% Flavio Mourao. Nucleo de Neurociencias NNC.
% email: mourao.fg@gmail.com
% Universidade Federal de Minas Gerais
% Started in:  05/2020
% Last update: 06/2020

%%
% First extract the data with: Extracting_LFPs_and_events.m
% ... then organize the data with the script: Pre_processing.m

%%

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

%% Run each session sequentially

%% Cut Time periods - pre behavior and during behavior events. Ignore "pos" period

% Choose channels
Granger.ch1 = 16;
Granger.ch2 = 12;
Granger.data = data.data_behavior{1, 1}([Granger.ch1 Granger.ch2],:,:);


% Time index -->  behavior epoch begins / ends

Granger.time_zero_idx = dsearchn(data.time_behavior',0); % time zero index. Behavior Start.

for ii = 1:length(data.events.behavior.TS_LFPsec)

    temp  = data.events.behavior.TS_LFPsec(ii,2) - data.events.behavior.TS_LFPsec(ii,1); % time end index.  Behavior End.
    Granger.time_end_idx(ii) = dsearchn(data.time_behavior',temp');

end


% from...all behavior epochs over time
% Initialize variables

Granger.data_before = nan(2,max(Granger.time_zero_idx)-1,parameters.behavior.NTrials);
Granger.data_during = nan(2,max(Granger.time_end_idx - Granger.time_zero_idx)+1,parameters.behavior.NTrials);

for jj = 1:parameters.behavior.NTrials
        
        Granger.data_before(:,1:(Granger.time_zero_idx-1),jj) = Granger.data(:,1:Granger.time_zero_idx-1,jj);              
        Granger.data_during(:,1:(Granger.time_end_idx(jj) - Granger.time_zero_idx)+1,jj) = Granger.data(:,Granger.time_zero_idx:Granger.time_end_idx(jj),jj); 

end

Granger.data = [Granger.data_before Granger.data_during];

clear( 'ii','jj','temp')

%% Parameters

% Time to analyse
Granger.Tpre = 5;
Granger.Tduring = 5;

% Granger prediction parameters
Granger.timewin = 1000; % in ms

% Overlap
Granger.overlap = .9;
Granger.overlap_t = (Granger.timewin - (Granger.overlap * Granger.timewin))/1000;

% Convert parameters to indices
Granger.timewin_points = round(Granger.timewin/(1000/parameters.srate));

% Temporal down-sample results (but not data!)
Granger.times2save = -Granger.Tpre:Granger.overlap_t:Granger.Tduring; % in seconds

% Convert requested times to indices. Considering all trials with the same time window
Granger.times2saveidx = dsearchn(data.time_behavior',Granger.times2save');


%% Test Bayes info criteria(BIC) for optimal model order at each time point

% initialize
Granger.bic = zeros(length(Granger.times2save),40); % Bayes info criteria (hard-coded to order=15)


for timei=1:length(Granger.times2save)
        
    % data from all trials in this time window
    Granger.tempdata = Granger.data(:,Granger.times2saveidx(timei)-floor(Granger.timewin_points/2):Granger.times2saveidx(timei)+floor(Granger.timewin_points/2)-mod(Granger.timewin_points+1,2),:);

    %  zscore all data 
    for triali=1:parameters.behavior.NTrials
        
        Granger.tempdata(1,:,triali) = zscore(detrend(Granger.tempdata(1,:,triali)));
        Granger.tempdata(2,:,triali) = zscore(detrend(Granger.tempdata(2,:,triali)));
   
    end
    
        Granger.tempdata = reshape(Granger.tempdata,2,Granger.timewin_points*parameters.behavior.NTrials);
    
    for bici=1:size(Granger.bic,2)
        
        % run model
        [Axy,E] = armorf(Granger.tempdata,parameters.behavior.NTrials,Granger.timewin_points,bici);
        
        % compute Bayes Information Criteria
        Granger.bic(timei,bici) = log(det(E)) + (log(length(Granger.tempdata))*bici*2^2)/length(Granger.tempdata);

    end
end

% Plot

figure
set(gcf,'color','white')

subplot(121)
plot((1:size(Granger.bic,2))*(1000/parameters.srate),mean(Granger.bic,1),'--.')
xlabel('Order (converted to ms)')
ylabel('Mean BIC over all time points')

[Granger.bestbicVal,Granger.bestbicIdx]=min(mean(Granger.bic,1));
hold on
plot(Granger.bestbicIdx*(1000/parameters.srate),Granger.bestbicVal,'mo','markersize',15)

title([ 'Optimal order is ' num2str(Granger.bestbicIdx) ' (' num2str(Granger.bestbicIdx*(1000/parameters.srate)) ' ms)' ])

subplot(122)
[junk,Granger.bic_per_timepoint] = min(Granger.bic,[],2);
plot(Granger.times2save,Granger.bic_per_timepoint*(1000/parameters.srate),'--.')
xlabel('Time (ms)')
ylabel('Optimal order (converted to ms)')

title('Optimal order (in ms) at each time point')

clear ('junk','Axy','bici','E','timei', 'triali')

%% Model for each trial

% Define model order in time points
 
Granger.order  =  Granger.bestbicIdx;

% initialize
[Granger.x2y,Granger.y2x] = deal(zeros(1,length(Granger.times2save))); % the function deal assigns inputs to all outputs
    
for timei=1:length(Granger.times2save)
    
    % data from all trials in this time window
    Granger.tempdata = Granger.data(:,Granger.times2saveidx(timei)-floor(Granger.timewin_points/2):Granger.times2saveidx(timei)+floor(Granger.timewin_points/2)-mod(Granger.timewin_points+1,2),:);

    % detrend and zscore all data
    for triali=1:parameters.behavior.NTrials
        
        Granger.tempdata(1,:,triali) = zscore(detrend(Granger.tempdata(1,:,triali)));
        Granger.tempdata(2,:,triali) = zscore(detrend(Granger.tempdata(2,:,triali)));
   
        % At this point with real data, you should check for stationarity
        % and possibly discard or mark data epochs that are extreme stationary violations.
   
    end
    
    Granger.tempdata = reshape(Granger.tempdata,2,Granger.timewin_points*parameters.behavior.NTrials);
        
    % fit AR models (model estimation from bsmart toolbox)
    
    % i.e = [Ax,Ex] = armorf(data,trials,timewin_points,order_points);
    % Ax = polynomial coefficients A corresponding to the AR model estimate of matrix X using Morf's method
    % Ex = prediction error E (the covariance matrix of the white noise of the AR model).
    
    [Ax,Ex] = armorf(Granger.tempdata(1,:),parameters.behavior.NTrials,Granger.timewin_points,Granger.order);
    [Ay,Ey] = armorf(Granger.tempdata(2,:),parameters.behavior.NTrials,Granger.timewin_points,Granger.order);
    [Axy,E] = armorf(Granger.tempdata     ,parameters.behavior.NTrials,Granger.timewin_points,Granger.order);
    
    % causal estimate reorganized in another variable as well
    Granger.y2x(timei)=log(Ex/E(1,1));
    Granger.x2y(timei)=log(Ey/E(2,2));
    
   
end

clear ('Ax','Axy','Ay','E','Ex','Ey','timei', 'triali')

% Plot

figure
set(gcf,'color','white')

plot(Granger.times2save,Granger.x2y,'Color','[0.5, 0.5, 0.5]','linew',2)
hold
plot(Granger.times2save,Granger.y2x,'k','linew',2)

title([ '\fontsize{14} Window length: ' num2str(Granger.timewin) ' ms, ' 'overlap: ' num2str(Granger.overlap*100) '%, ' 'order: ' num2str(Granger.order*(1000/parameters.srate)) ' ms' ])
xlabel('\fontsize{14} Time (s)')
ylabel('\fontsize{14} Granger prediction estimate')
xlim ([-5 5])

legend('\fontsize{14}PFC -> AMY','\fontsize{14}AMY -> PFC','box','off')

box off

%% last update 01/06/2020 - 2:38am
%  listening: Snowman - WYS
