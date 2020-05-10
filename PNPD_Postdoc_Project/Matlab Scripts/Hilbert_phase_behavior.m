
%% Phase analyses based on Hilbert Transform
%  Considering the behavior events

% The code relies on the following package:
% --> circular-statistics-toolbox
% https://github.com/circstat/circstat-matlab 


% Flavio Mourao. Nucleo de Neurociencias NNC.
% email: mourao.fg@gmail.com
% Universidade Federal de Minas Gerais
% Started in:  04/2019
% Last update: 04/2020

%% Run each session sequentially

%%
% Hilbert Transform

% Loop over channels and make hilbert transform

% Initializing with NaN
hilb.behavior.coefficients = cell(size(data.data_behavior));
hilb.behavior.coefficients(:,:) = {nan(size(data.data_behavior{1,1}))};

for ii = 1:size(data.data,2)
    for jj = 1:size(data.data_behavior{1, ii},3)
        for ll = 1:size(data.data_behavior{1, ii},1)
        temp = data.data_behavior{1, ii}(ll,:,jj);
        temp(isnan(temp(:,:)))=[]; % it was necessary to remove the nan for the function to work
        hilb.behavior.coefficients{1, ii}(ll,1:length(temp),jj) = hilbert(temp);
        end
    end 
end 

% Extract phase

hilb.behavior.phase = cell(size(hilb.behavior.coefficients));
hilb.behavior.energy = cell(size(hilb.behavior.coefficients));

for ii = 1:size(hilb.behavior.coefficients,2)
    hilb.behavior.phase{1,ii}  = angle(hilb.behavior.coefficients{1,ii});
    hilb.behavior.energy{1,ii} = abs(hilb.behavior.coefficients{1,ii});
end


clear('ii','jj','ll','temp')

%%  - CHANNELS MAP - 
% *just in case to check

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

%% Delta phase from Euler representation of angles

% Difference between channels
% all possible combinations between all channels:

% hilb.behavior.delta_phase_all = cell matriz 17x17.
%                               each cell -> rows    -> trials
%                                            columns -> times

% Create a channels grid for combinations
% recording channels (rows 2->17) and CS modulating signal (row 1)
[X,Y] = meshgrid((1:parameters.nch+1),(1:parameters.nch+1));

% Choose filter column to analyse
ff= 3;

% Initilize variable
hilb.behavior.delta_phase = cell(parameters.nch+1);
hilb.behavior.delta_phase(:,:) = {complex(nan(size(hilb.behavior.phase{1,1},3),size(hilb.behavior.phase{1,1},2)))};

% Delta phase for each trial
for jj = 1:size(X,1)
    for ii = 1:size(Y,1)
        hilb.behavior.delta_phase{jj,ii} = squeeze(exp(1i*(hilb.behavior.phase{1,ff}(jj,:,:) - hilb.behavior.phase{1,ff}(ii,:,:))))';
    end
end 

clear('ff','X','Y','ch','jj','ii')

%% Extracts relative phase and length of circular variance (PLV) 
%  Measure of phase synchronization

% Time window
hilb.behavior.time_window     =  .5; % sec.
hilb.behavior.time_window_idx = round(hilb.behavior.time_window*parameters.srate);

% Overlap
hilb.behavior.timeoverlap    = .9; % percentage
overlap = round((hilb.behavior.time_window_idx)-(hilb.behavior.timeoverlap*hilb.behavior.time_window_idx));

% Time epochs
time2save_idx = (1:overlap:length(hilb.behavior.delta_phase{1,1})-hilb.behavior.time_window_idx);

hilb.behavior.phase_win              = cell(size(hilb.behavior.delta_phase)); % all events
hilb.behavior.PLV_win                = cell(size(hilb.behavior.delta_phase)); % all events
hilb.behavior.stats.phase_win_mean   = cell(size(hilb.behavior.delta_phase)); % mean events
hilb.behavior.stats.PLV_win_mean  	 = cell(size(hilb.behavior.delta_phase)); % mean events

hilb.behavior.stats.PLV_win_median   = cell(size(hilb.behavior.delta_phase)); % median events

for ii = 1:size(hilb.behavior.delta_phase,2)
    for jj = 1:size(hilb.behavior.delta_phase,1)
        for ll = 1:length(time2save_idx)
            
        temp1   = nanmean(hilb.behavior.delta_phase{ii,jj}(:,time2save_idx(ll):(time2save_idx(ll) + hilb.behavior.time_window_idx - 1)),2);
        
        hilb.behavior.phase_win{ii,jj}(:,ll) = angle(temp1); % time epoch for all trials over time
        hilb.behavior.PLV_win{ii,jj}(:,ll)   = abs(temp1);   % time epoch for all trials over time
        
        hilb.behavior.stats.phase_win_mean{ii,jj}(:,ll)  = angle(nanmean(temp1,1)); % one value for each time epoch. Average relative phase
        hilb.behavior.stats.PLV_win_mean{ii,jj}(:,ll)    = nanmean(abs(temp1),1);   % one value for each time epoch. Average PLV
        
        hilb.behavior.stats.PLV_win_median{ii,jj}(:,ll)  = nanmedian(abs(temp1),1); % one value for each time epoch. Median PLV
        
        end
    end
end 

% Time vector
hilb.behavior.time = (linspace(-parameters.behavior.Tpre,parameters.behavior.trialperiod+parameters.behavior.Tpos,size(hilb.behavior.phase_win{1, 1},2)));


clear ('overlap','time2save_idx','ii','jj','ll','temp1')

%% Plot to check average angles in the sliding window

figure

% choose par channels to compare
ch = [2 11];

% Time to plot
t = [-2 5];

suptitle({['\Delta Phase over trials' ' / ' 'Time window = ' num2str(hilb.behavior.time_window) 's ' '- ' 'Overlap = ' num2str(hilb.behavior.timeoverlap*100) '%'];[]}) 
set(gcf,'color','white')

sb = 1; %subplot number

for ii = 1:size(ch,1)
    for jj = 1:parameters.behavior.NTrials        
        subplot(size(ch,1),parameters.behavior.NTrials,sb)
        
        plot(hilb.behavior.time, rad2deg(hilb.behavior.phase_win{ch(ii,1), ch(ii,2)}(jj,:)),'k','linew',1)
        
        hold all
        
        plot([0 0],[-180 180],'r--')
        plot([data.events.behavior.TS_LFPsec(jj,2) - data.events.behavior.TS_LFPsec(jj,1) data.events.behavior.TS_LFPsec(jj,2) - data.events.behavior.TS_LFPsec(jj,1)],[-180 180],'r--')

        ylim([-180 180])       
        %xlim([ -parameters.behavior.Tpre (data.events.behavior.TS_LFPsec(jj,2) - data.events.behavior.TS_LFPsec(jj,1)) + parameters.behavior.Tpos])      
        xlim(t)
        
        title (['Event ' num2str(jj)])        
        
        xlabel('Time (s)')      
        ylabel('\Delta Phase (^{o})')
        
        sb = sb + 1;
    end
end     

figure
       
suptitle({['\Delta Phase over trials' ' / ' 'Time window = ' num2str(hilb.behavior.time_window) 's ' '- ' 'Overlap = ' num2str(hilb.behavior.timeoverlap*100) '%'];[]}) 
set(gcf,'color','white')

for ii = 1:size(ch,1)
    subplot(1,size(ch,1),ii)
    
    plot(hilb.behavior.time,(rad2deg(hilb.behavior.stats.phase_win_mean{ch(ii,1), ch(ii,2)})),'k','linew',1)
    
    hold
    plot([0 0],[-180 180],'r--')

    ylim([-180 180])
    xlim(t)
        
    %title (titles{ii})    
    title (['Channels: ' num2str(ch(ii,1)) ' <-> ' num2str(ch(ii,2))])
    
    xlabel('Time (s)')
    ylabel('\Delta Phase (^{o})')
end

clear ('ch','sb','ii','jj','t')

%% Plot to check PLV in the sliding window

figure

% choose par channels to compare
ch = [2 11];

% Time to plot
t = [ -3 10];

suptitle({[' PLV over trials' ' / ' 'Time window = ' num2str(hilb.behavior.time_window) 's ' '- ' 'Overlap = ' num2str(hilb.behavior.timeoverlap*100) '%'];[]}) 
set(gcf,'color','white')

sb = 1; %subplot counter

for ii = 1:size(ch,1)
    for jj = 1:parameters.behavior.NTrials        
        subplot(size(ch,1),parameters.behavior.NTrials,sb)
        
        plot(hilb.behavior.time, hilb.behavior.PLV_win{ch(ii,1), ch(ii,2)}(jj,:),'k','linew',1)
       
        hold all
        plot([0 0],[0 1],'r--')
        plot([data.events.behavior.TS_LFPsec(jj,2) - data.events.behavior.TS_LFPsec(jj,1) data.events.behavior.TS_LFPsec(jj,2) - data.events.behavior.TS_LFPsec(jj,1)],[0 1],'r--')
        
        ylim([0 1])      
        %xlim([ -parameters.behavior.Tpre (data.events.behavior.TS_LFPsec(jj,2) - data.events.behavior.TS_LFPsec(jj,1)) + parameters.behavior.Tpos])
        xlim(t)
        
        title (['Event ' num2str(jj)])        
        
        xlabel('Time (s)')
        ylabel('Phase Synchronization')
        
        sb = sb + 1;
    end
end     

figure
       
suptitle({[' PLV over trials' ' / ' 'Time window = ' num2str(hilb.behavior.time_window) 's ' '- ' 'Overlap = ' num2str(hilb.behavior.timeoverlap*100) '%'];[]}) 
set(gcf,'color','white')

for ii = 1:size(ch,1)
    subplot(1,size(ch,1),ii)
    
    plot(hilb.behavior.time,hilb.behavior.stats.PLV_win_mean{ch(ii,1), ch(ii,2)},'k','linew',1)    
    
    hold all
    plot([0 0],[0 1],'r--')

    ylim([0 1])
    xlim(t)
    
    title (['Channels: ' num2str(ch(ii,1)) ' <-> ' num2str(ch(ii,2))])
    
    xlabel('Time (s)')
    ylabel('Phase Synchronization')
end

clear ('ch','sb','ii','jj','t')

%% Plot to check PLV in a color map

% Choose the channel that will be compared with all the others
ch = 11;
ch_compare = cell2mat(hilb.behavior.stats.PLV_win_mean(:,ch));

% Choose channels to plot
chs = [2 16];

% Time to plot
t = [ -3 10];

figure
set(gcf,'color','white')

subplot 121
contourf(hilb.behavior.time,1:length(chs),ch_compare(chs,:),80,'linecolor','none');
xlabel('Time (s)','FontSize',14), ylabel('channels','FontSize',14)
set(gca,'yticklabel',[])
xlim(t)
colorbar
caxis([.3 .8])

subplot 122
imagesc(hilb.behavior.time,chs,ch_compare(chs,:))
xlabel('Time (s)','FontSize',14), ylabel('channels','FontSize',14)
set(gca,'yticklabel',[])
xlim(t)
colorbar
caxis([.3 .8])

clear ('ch','chs','ch_compare','t')

%% Cut Time periods - pre behavior and behavior events. Ignore "pos" period

% Time index -->  behavior epoch begins / ends

hilb.behavior.time_zero_idx = dsearchn(hilb.behavior.time',0'); % time zero index. Behavior Start.

for ii = 1:length(data.events.behavior.TS_LFPsec)

    temp  = data.events.behavior.TS_LFPsec(ii,2) - data.events.behavior.TS_LFPsec(ii,1); % time end index.  Behavior End.
    hilb.behavior.time_end_idx(ii) = dsearchn(hilb.behavior.time',temp');

end


% from...all behavior epochs over time
% Initialize variables
hilb.behavior.phase_win_before      = cell(length(hilb.behavior.delta_phase),length(hilb.behavior.delta_phase));
hilb.behavior.phase_win_before(:,:) = {nan(parameters.behavior.NTrials,max(hilb.behavior.time_zero_idx)-1)};

hilb.behavior.phase_win_during      = cell(length(hilb.behavior.delta_phase),length(hilb.behavior.delta_phase));
hilb.behavior.phase_win_during(:,:) = {nan(parameters.behavior.NTrials,max(hilb.behavior.time_end_idx - hilb.behavior.time_zero_idx)+1)};

hilb.behavior.PLV_win_before        = cell(length(hilb.behavior.delta_phase),length(hilb.behavior.delta_phase));
hilb.behavior.PLV_win_before(:,:)   = {nan(parameters.behavior.NTrials,max(hilb.behavior.time_zero_idx)-1)};

hilb.behavior.PLV_win_during        = cell(length(hilb.behavior.delta_phase),length(hilb.behavior.delta_phase));
hilb.behavior.PLV_win_during(:,:)   = {nan(parameters.behavior.NTrials,max(hilb.behavior.time_end_idx - hilb.behavior.time_zero_idx)+1)};

for ii = 1:length(hilb.behavior.delta_phase)*length(hilb.behavior.delta_phase)
    for jj = 1:parameters.behavior.NTrials
        
        hilb.behavior.phase_win_before{ii}(jj,1:(hilb.behavior.time_zero_idx-1)) = hilb.behavior.phase_win{ii}(jj,1:hilb.behavior.time_zero_idx-1);              
        hilb.behavior.phase_win_during{ii}(jj,1:(hilb.behavior.time_end_idx(jj) - hilb.behavior.time_zero_idx)+1)  = hilb.behavior.phase_win{ii}(jj,hilb.behavior.time_zero_idx:hilb.behavior.time_end_idx(jj)); 

        hilb.behavior.PLV_win_before{ii}(jj,1:(hilb.behavior.time_zero_idx-1))    = hilb.behavior.PLV_win{ii}(jj,1:hilb.behavior.time_zero_idx-1);                
        hilb.behavior.PLV_win_during{ii}(jj,1:(hilb.behavior.time_end_idx(jj) - hilb.behavior.time_zero_idx)+1)    = hilb.behavior.PLV_win{ii}(jj,hilb.behavior.time_zero_idx:hilb.behavior.time_end_idx(jj));
    
    end
end

clear('ii', 'jj', 'temp')

%% Average behavior epochs and total sessions

% define time window to analyse

% - IMPORTANT -
% Whereas the behavioral events are asymmetric in the time domain, 
% consider the maximum  time window to analysis according to the minimum behavior window. 
% In other words, smaller behavior window = 5s, this will be the maximum
% time to analysis during the behavior period. This avoids asymmetries in the average time

% choose time window
t = 3;
time_zero_idx = (dsearchn(hilb.behavior.time',t')) - hilb.behavior.time_zero_idx-1;


% Relative Phase and PLV - mean behavior events over time
hilb.behavior.stats.phase_win_mean_before = cell(length(hilb.behavior.delta_phase),length(hilb.behavior.delta_phase));
hilb.behavior.stats.phase_win_mean_during = cell(length(hilb.behavior.delta_phase),length(hilb.behavior.delta_phase));

hilb.behavior.stats.PLV_win_mean_before   = cell(length(hilb.behavior.delta_phase),length(hilb.behavior.delta_phase));
hilb.behavior.stats.PLV_win_mean_during   = cell(length(hilb.behavior.delta_phase),length(hilb.behavior.delta_phase));

for ii = 1:length(hilb.behavior.delta_phase)*length(hilb.behavior.delta_phase)
        
    hilb.behavior.stats.phase_win_mean_before{ii}  = circ_mean(hilb.behavior.phase_win_before{ii},[],1);                
    hilb.behavior.stats.phase_win_mean_during{ii}  = circ_mean(hilb.behavior.phase_win_during{ii}(:,1:time_zero_idx),[],1);

    hilb.behavior.stats.PLV_win_mean_before{ii}    = mean(hilb.behavior.PLV_win_before{ii},1);                   
    hilb.behavior.stats.PLV_win_mean_during{ii}    = mean(hilb.behavior.PLV_win_during{ii}(:,1:time_zero_idx),1);

end


% Relative Phase and PLV - one mean value for each behavior event
hilb.behavior.stats.phase_mean_before = cell(length(hilb.behavior.delta_phase),length(hilb.behavior.delta_phase));
hilb.behavior.stats.phase_mean_during = cell(length(hilb.behavior.delta_phase),length(hilb.behavior.delta_phase));

hilb.behavior.stats.PLV_mean_before   = cell(length(hilb.behavior.delta_phase),length(hilb.behavior.delta_phase));
hilb.behavior.stats.PLV_mean_during   = cell(length(hilb.behavior.delta_phase),length(hilb.behavior.delta_phase));

for ii = 1:length(hilb.behavior.delta_phase)*length(hilb.behavior.delta_phase)
        
    hilb.behavior.stats.phase_mean_before{ii}  = circ_mean(hilb.behavior.phase_win_before{ii},[],2);                
    hilb.behavior.stats.phase_mean_during{ii}  = circ_mean(hilb.behavior.phase_win_during{ii}(:,1:time_zero_idx),[],2);

    hilb.behavior.stats.PLV_mean_before{ii}    = nanmean(hilb.behavior.PLV_win_before{ii},2);                   
    hilb.behavior.stats.PLV_mean_during{ii}    = nanmean(hilb.behavior.PLV_win_during{ii}(:,1:time_zero_idx),2);
   
end


% Relative Phase and PLV - Only one mean value for the whole session
hilb.behavior.stats.phase_Total_mean_before = zeros(length(hilb.behavior.delta_phase),length(hilb.behavior.delta_phase));
hilb.behavior.stats.phase_Total_mean_during = zeros(length(hilb.behavior.delta_phase),length(hilb.behavior.delta_phase));

hilb.behavior.stats.PLV_Total_mean_before   = zeros(length(hilb.behavior.delta_phase),length(hilb.behavior.delta_phase));
hilb.behavior.stats.PLV_Total_mean_during   = zeros(length(hilb.behavior.delta_phase),length(hilb.behavior.delta_phase));

for ii = 1:length(hilb.behavior.delta_phase)*length(hilb.behavior.delta_phase)
    
    hilb.behavior.stats.phase_Total_mean_before(ii)    = circ_mean(hilb.behavior.stats.phase_mean_before{ii},[],1); % one value for each session. Average events during pre behavior period
    hilb.behavior.stats.phase_Total_mean_during(ii)    = circ_mean(hilb.behavior.stats.phase_mean_during{ii},[],1); % one value for each session. Average events during behavior period
    
    hilb.behavior.stats.PLV_Total_mean_before(ii)      = mean(hilb.behavior.stats.PLV_mean_before{ii},1); % one value for each session. Average events during pre behavior period
    hilb.behavior.stats.PLV_Total_mean_during(ii)      = mean(hilb.behavior.stats.PLV_mean_during{ii},1); % one value for each session. Average events during behavior period

end

clear('ii')

%% Plot polar plots - Phase Coherence values - Each trial

figure

% choose par channels to compare
ch = [3 12];

suptitle({'\Delta Phase average over trials. Before and during behavior event';['Time window = ' num2str(hilb.behavior.time_window) 's ' '- ' 'Overlap = ' num2str(hilb.behavior.timeoverlap*100) '%'];[]}) 
set(gcf,'color','white')

sb = 1; %subplot counter

for ii = 1:size(ch,1)
    for jj = 1:parameters.behavior.NTrials 
        subplot(2,parameters.behavior.NTrials,sb)

        polarplot([zeros(size(hilb.behavior.phase_win_before{ch(ii,1), ch(ii,2)}(jj,:))), hilb.behavior.phase_win_before{ch(ii,1), ch(ii,2)}(jj,:)]',repmat([0 1],1,length(hilb.behavior.phase_win_before{ch(ii,1), ch(ii,2)}(jj,:)))','k');
        hold all
        polarplot([0,hilb.behavior.stats.phase_mean_before{ch(ii,1), ch(ii,2)}(jj)]',[0 1]','Color','[0.6350, 0.0780, 0.1840]','linew',2);
        
        sb = sb + 1;
    end
end  


for ii = 1:size(ch,1)
    for jj = 1:parameters.behavior.NTrials 
        subplot(2,parameters.behavior.NTrials,sb)

        polarplot([zeros(size(hilb.behavior.phase_win_during{ch(ii,1), ch(ii,2)}(jj,1:time_zero_idx))), hilb.behavior.phase_win_during{ch(ii,1), ch(ii,2)}(jj,1:time_zero_idx)]',repmat([0 1],1,length(hilb.behavior.phase_win_during{ch(ii,1), ch(ii,2)}(jj,1:time_zero_idx)))','k');
        hold all
        polarplot([0,hilb.behavior.stats.phase_mean_during{ch(ii,1), ch(ii,2)}(jj)]',[0 1]','Color','[0.6350, 0.0780, 0.1840]','linew',2);
        
        sb = sb + 1;
    end
end

clear ('ii','jj','sb','ch')

%% Plot polar plots - Phase Coherence values - Total events mean over time

figure

% choose par channels to compare
ch = [2 16 ; 3 16 ; 6 16 ; 12 16];

suptitle({'Total \Delta Phase average over time.';['Time window = ' num2str(hilb.behavior.time_window) 's ' '- ' 'Overlap = ' num2str(hilb.behavior.timeoverlap*100) '%'];[]}) 
set(gcf,'color','white')

for ii = 1:length(ch)

    subplot(2,length(ch),ii)

    polarplot([zeros(size(hilb.behavior.stats.phase_win_mean_before{ch(ii,1), ch(ii,2)})), hilb.behavior.stats.phase_win_mean_before{ch(ii,1), ch(ii,2)}]',repmat([0 1],1,length(hilb.behavior.stats.phase_win_mean_before{ch(ii,1), ch(ii,2)}))','k');
    hold all
    polarplot([0,hilb.behavior.stats.phase_Total_mean_before(ch(ii,1), ch(ii,2))]',[0 1]','Color','[0.6350, 0.0780, 0.1840]','linew',2);
    
end  

for ii = 1:length(ch)

    subplot(2,length(ch),ii+length(ch))

    polarplot([zeros(size(hilb.behavior.stats.phase_win_mean_during{ch(ii,1), ch(ii,2)})), hilb.behavior.stats.phase_win_mean_during{ch(ii,1), ch(ii,2)}]',repmat([0 1],1,length(hilb.behavior.stats.phase_win_mean_during{ch(ii,1), ch(ii,2)}))','k');
    hold all
    polarplot([0,hilb.behavior.stats.phase_Total_mean_during(ch(ii,1), ch(ii,2))]',[0 1]','Color','[0.6350, 0.0780, 0.1840]','linew',2);
      
end

clear ('ii','jj','sb','ch')

%% Plot all channels - Phase Coherence values in a color map - Each event

% Choose channels to plot
chs = [3 6 12 16];

figure
set(gcf,'color','white')

temp1 = zeros(length(hilb.behavior.delta_phase));
temp2 = zeros(length(hilb.behavior.delta_phase));

for ii = 1:parameters.behavior.NTrials 
    for jj = 1:length(hilb.behavior.delta_phase)*length(hilb.behavior.delta_phase)  
        
        temp1(jj) = hilb.behavior.stats.PLV_mean_before{jj}(ii);
        temp2(jj) = hilb.behavior.stats.PLV_mean_during{jj}(ii);

    end 
    
        subplot (2,parameters.behavior.NTrials,ii)
        imagesc(chs,chs,temp1(chs,chs))
        xlabel('channels','FontSize',14), ylabel('channels','FontSize',14)
        colorbar
        caxis([0 1])
        
        subplot (2,parameters.behavior.NTrials,ii+parameters.behavior.NTrials)
        imagesc(chs,chs,temp2(chs,chs))
        xlabel('channels','FontSize',14), ylabel('channels','FontSize',14)
        colorbar
        caxis([0 1])
        
end

clear ('temp1','temp2','chs')

%% Plot all channels - Phase Coherence values in a color map - Total Session

% Choose channels to plot
chs = [3 6 12 16];

figure
set(gcf,'color','white')

subplot 121
imagesc(1:length(chs),1:length(chs),hilb.behavior.stats.PLV_Total_mean_before(chs,chs))
xlabel('channels','FontSize',14), ylabel('channels','FontSize',14)
colorbar
box off
caxis([0 1])
set (gca,'visible','off')

subplot 122
imagesc(1:length(chs),1:length(chs),hilb.behavior.stats.PLV_Total_mean_during(chs,chs))
xlabel('channels','FontSize',14), ylabel('channels','FontSize',14)
colorbar
box off
caxis([0 1])
set (gca,'visible','off')

clear ('chs')

%%
% save('')

%% last update 25/04/2020 - 2:19am
%  listening: Aurora - life on mars

