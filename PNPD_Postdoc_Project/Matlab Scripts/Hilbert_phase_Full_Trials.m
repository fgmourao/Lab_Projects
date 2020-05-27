
%% Phase analyses based on Hilbert Transform
%  - Performs analysis considering the trial periods - sound on/off

%  Considering the entire trial period

% The code relies on the following package:
% --> circular-statistics-toolbox
%     https://github.com/circstat/circstat-matlab 


% Flavio Mourao. Nucleo de Neurociencias NNC.
% email: mourao.fg@gmail.com
% Universidade Federal de Minas Gerais
% Started in:  02/2019
% Last update: 04/2020

%%
% First extract the data with: Extracting_LFPs_and_events.m
% ... then organize the data with the script: Pre_processing.m

%% Run each session sequentially

%% Hilbert Transform

% Loop over channels and make hilbert transform

% Initializing with NaN
hilb.coefficients_trials = cell(size(data.data_trials));
hilb.coefficients_trials(:,:) = {nan(size(data.data_trials{1,1}))};

for ii = 1:size(data.data_trials,2)
    for jj = 1:size(data.data_trials{1, ii},3)
        for ll = 1:size(data.data_trials{1, ii},1)
        temp = data.data_trials{1, ii}(ll,:,jj);
        temp(isnan(temp(:,:)))=[]; % it was necessary to remove the NaN for the function to work
        hilb.coefficients_trials{1, ii}(ll,1:length(temp),jj) = hilbert(temp);
        end
    end 
end 

% Extract phase

hilb.phase_trials  = cell(size(hilb.coefficients_trials));
hilb.energy_trials = cell(size(hilb.coefficients_trials));

for ii = 1:size(hilb.coefficients_trials,2)
    hilb.phase_trials{1,ii}   = angle(hilb.coefficients_trials{1,ii});
    hilb.energy_trials{1,ii}  = abs(hilb.coefficients_trials{1,ii});
end


clear('ii','jj','ll','temp')

%% - CHANNELS MAP - 
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

% hilb.delta_phase_trials = cell matriz 17x17.
%                               each cell -> rows    -> trials
%                                            columns -> times

% Create a channels grid for combinations
% recording channels (rows 2->17) and CS modulating signal (row 1)
[X,Y] = meshgrid((1:parameters.nch+1),(1:parameters.nch+1));

% Choose frequency bands index according to the Pre_processing.m

% 2  - modulator / 3 - delta / 4 - lowtheta / 5 - hightheta
% 6  - alpha     / 7 - beta  / 8 - lowgamma / 9 - highgamma
% 10 - extra

ff= 2;

% Initilize variable
hilb.delta_phase_trials = cell(parameters.nch+1);
hilb.delta_phase_trials(:,:) = {complex(nan(size(hilb.phase_trials{1,1},3),size(hilb.phase_trials{1,1},2)))};

% Delta phase for each trial
for jj = 1:size(X,1)
    for ii = 1:size(Y,1)
        hilb.delta_phase_trials{jj,ii} = squeeze(exp(1i*(hilb.phase_trials{1,ff}(jj,:,:) - hilb.phase_trials{1,ff}(ii,:,:))))';
    end
end 

clear('ff','X','Y','ch','jj','ii')

%% Extracts relative phase and length of circular variance (PLV) 
%  Measure of phase synchronization

% Time window
hilb.time_window     =  2; % sec.
hilb.time_window_idx = round(hilb.time_window*parameters.srate);

% Overlap
hilb.timeoverlap    = .8; % percentage
overlap = round((hilb.time_window_idx)-(hilb.timeoverlap*hilb.time_window_idx));

% Time epochs
time2save_idx = (1:overlap:length(hilb.delta_phase_trials{1,1})-hilb.time_window_idx);

hilb.phase_win_trials              = cell(size(hilb.delta_phase_trials)); % all trials
hilb.PLV_win_trials                = cell(size(hilb.delta_phase_trials)); % all trials
hilb.stats.phase_win_mean_trials   = cell(size(hilb.delta_phase_trials)); % mean trials
hilb.stats.PLV_win_mean_trials     = cell(size(hilb.delta_phase_trials)); % mean trials

hilb.stats.PLV_win_median_trials   = cell(size(hilb.delta_phase_trials)); % median trials

for ii = 1:size(hilb.delta_phase_trials,2)
    for jj = 1:size(hilb.delta_phase_trials,1)
        for ll = 1:length(time2save_idx)
            
        temp1   = mean(hilb.delta_phase_trials{ii,jj}(:,time2save_idx(ll):(time2save_idx(ll) + hilb.time_window_idx -1)),2);
        
        hilb.phase_win_trials{ii,jj}(:,ll) = angle(temp1); % time epoch for all trials over time
        hilb.PLV_win_trials{ii,jj}(:,ll)   = abs(temp1);   % time epoch for all trials over time
        
        hilb.stats.phase_win_mean_trials{ii,jj}(:,ll) = angle(mean(temp1,1)); % one value for each time epoch. Average relative phase
        hilb.stats.PLV_win_mean_trials{ii,jj}(:,ll)   = mean(abs(temp1),1);   % one value for each time epoch. Average PLV
        
        hilb.stats.PLV_win_median_trials{ii,jj}(:,ll) = median(abs(temp1),1); % one value for each time epoch. Median PLV

        end
    end
end 

% Time vector
hilb.time_trials = (linspace(-parameters.Tpre,parameters.trialperiod+parameters.Tpos,size(hilb.phase_win_trials{1, 1},2)));

clear ('overlap','time2save_idx','ii','jj','ll','temp1')

%% Plot to check average angles in the sliding window

figure

% choose par channels to compare
ch = [3 12 ; 16 11];

suptitle({'\Delta Phase over trials';['Time window = ' num2str(hilb.time_window) 's ' '- ' 'Overlap = ' num2str(hilb.timeoverlap*100) '%'];[];[]}) 
set(gcf,'color','white')

sb = 1; %subplot number

for ii = 1:size(ch,1)
    for jj = 1:parameters.NTrials        
        subplot(size(ch,1),parameters.NTrials,sb)
        
        plot(hilb.time_trials, rad2deg(hilb.phase_win_trials{ch(ii,1), ch(ii,2)}(jj,:)),'k','linew',1)
        
        hold all
        
        plot([0 0],[-180 180],'r--')
        plot([30 30],[-180 180],'r--')

        set(gca,'xlim',[hilb.time_trials(1) hilb.time_trials(end)],'ylim',[-180 180])
        
        title (['Trial ' num2str(jj)])  
        
        xlabel('Time (s)')
        ylabel('\Delta Phase (^{o})')
        
        sb = sb + 1;
    end
end     

figure
       
 suptitle({'\Delta Phase average between trials';['Time window = ' num2str(hilb.time_window) 's ' '- ' 'Overlap = ' num2str(hilb.timeoverlap*100) '%'];[]}) 
 set(gcf,'color','white')

for ii = 1:size(ch,1)
    subplot(1,size(ch,1),ii)
    
    plot(hilb.time_trials,(rad2deg(hilb.stats.phase_win_mean_trials{ch(ii,1), ch(ii,2)})),'k','linew',1)
    
    hold
    plot([0 0],[-180 180],'r--')
    plot([30 30],[-180 180],'r--')
    
    set(gca,'xlim',[hilb.time_trials(1) hilb.time_trials(end)],'ylim',[-180 180])
    
    title (['Channels: ' num2str(ch(ii,1)) ' <-> ' num2str(ch(ii,2))])
    
    xlabel('Time (s)')
    ylabel('\Delta Phase (^{o})')
end

clear ('ch','sb','ii','jj')

%% Plot to check PLV in the sliding window

figure

% choose par channels to compare
ch = [3 12 ; 16 1];

suptitle({[' PLV over trials'];['Time window = ' num2str(hilb.time_window) 's ' '- ' 'Overlap = ' num2str(hilb.timeoverlap*100) '%'];[]}) 
set(gcf,'color','white')

sb = 1; %subplot counter

for ii = 1:size(ch,1)
    for jj = 1:parameters.NTrials        
        subplot(size(ch,1),parameters.NTrials,sb)
        
        plot(hilb.time_trials, hilb.PLV_win_trials{ch(ii,1), ch(ii,2)}(jj,:),'k','linew',1)
        
        hold all
        plot([0 0],[0 1],'r--')
        plot([30 30],[0 1],'r--')
        
        set(gca,'xlim',[hilb.time_trials(1) hilb.time_trials(end)],'ylim',[0 1])
        
        title (['Trial ' num2str(ii)])
        
        xlabel('Time (s)')
        ylabel('Phase Synchronization')
        
        sb = sb + 1;
    end
end     

figure
       
suptitle({'PLV average between trials';['Time window = ' num2str(hilb.time_window) 's ' '- ' 'Overlap = ' num2str(hilb.timeoverlap*100) '%'];[]}) 
set(gcf,'color','white')

for ii = 1:size(ch,1)
    subplot(1,size(ch,1),ii)
    
    plot(hilb.time_trials,hilb.stats.PLV_win_mean_trials{ch(ii,1), ch(ii,2)},'k','linew',1)    
    
    hold all
    plot([0 0],[0 1],'r--')
    plot([30 30],[0 1],'r--')
    
    set(gca,'xlim',[hilb.time_trials(1) hilb.time_trials(end)],'ylim',[0 1])
    
    title (['Channels: ' num2str(ch(ii,1)) ' <-> ' num2str(ch(ii,2))])
    
    xlabel('Time (s)')
    ylabel('Phase Synchronization')
end

clear ('ch','sb','ii','jj')

%% Plot to check PLV in a color map

% Choose the channel that will be compared with all the others
ch = 16;
ch_compare = cell2mat(hilb.stats.PLV_win_mean_trials(:,ch));

% Choose channels to plot
chs = [2 5 12 1];

figure
set(gcf,'color','white')

subplot 121
contourf(hilb.time_trials,1:length(chs),ch_compare(chs,:),80,'linecolor','none');
xlabel('Time (s)','FontSize',14), ylabel('channels','FontSize',14)
set(gca,'yticklabel',[])
xlim([-30 40])
colorbar
caxis([0 0.5])

subplot 122
imagesc(hilb.time_trials,chs,ch_compare(chs,:))
xlabel('Time (s)','FontSize',14), ylabel('channels','FontSize',14)
set(gca,'yticklabel',[])
xlim([-30 40])
colorbar
caxis([0 0.5])

clear ('ch','chs','ch_compare')

%% Cut Time vectors - pre trials and trials periods. Ignore "pos" period

% Time index --> 0s: sound sound begins / 30s: sound ends
hilb.time_zero_idx = dsearchn(hilb.time_trials',0'); % time zero index. Sound Start.
hilb.time_end_idx  = dsearchn(hilb.time_trials',30'); % time  30s index. Sound End.


% from...all trials over time
% Initialize variables
hilb.phase_win_trials_Presound      = cell(length(hilb.delta_phase_trials),length(hilb.delta_phase_trials));
hilb.phase_win_trials_Presound(:,:) = {nan(parameters.NTrials,max(hilb.time_zero_idx)-1)};

hilb.phase_win_trials_Sound         = cell(length(hilb.delta_phase_trials),length(hilb.delta_phase_trials));
hilb.phase_win_trials_Sound(:,:)    = {nan(parameters.NTrials,max(hilb.time_end_idx - hilb.time_zero_idx)+1)};

hilb.PLV_win_trials_Presound        = cell(length(hilb.delta_phase_trials),length(hilb.delta_phase_trials));
hilb.PLV_win_trials_Presound(:,:)   = {nan(parameters.NTrials,max(hilb.time_zero_idx)-1)};

hilb.PLV_win_trials_Sound           = cell(length(hilb.delta_phase_trials),length(hilb.delta_phase_trials));
hilb.PLV_win_trials_Sound(:,:)      = {nan(parameters.NTrials,max(hilb.time_end_idx - hilb.time_zero_idx)+1)};

for ii = 1:length(hilb.delta_phase_trials)*length(hilb.delta_phase_trials)
    hilb.phase_win_trials_Presound{ii}   = hilb.phase_win_trials{ii}(:,1:hilb.time_zero_idx-1);               
    hilb.phase_win_trials_Sound{ii}      = hilb.phase_win_trials{ii}(:,hilb.time_zero_idx:hilb.time_end_idx);

    hilb.PLV_win_trials_Presound{ii}    = hilb.PLV_win_trials{ii}(:,1:hilb.time_zero_idx-1);                
    hilb.PLV_win_trials_Sound{ii}       = hilb.PLV_win_trials{ii}(:,hilb.time_zero_idx:hilb.time_end_idx);
end

%% Average behavior epochs and total sessions

% Relative Phase and PLV - mean behavior events over time
hilb.stats.phase_win_mean_trials_Presound = cell(length(hilb.delta_phase_trials),length(hilb.delta_phase_trials));
hilb.stats.phase_win_mean_trials_Sound    = cell(length(hilb.delta_phase_trials),length(hilb.delta_phase_trials));

hilb.stats.PLV_win_mean_trials_Presound   = cell(length(hilb.delta_phase_trials),length(hilb.delta_phase_trials));
hilb.stats.PLV_win_mean_trials_Sound      = cell(length(hilb.delta_phase_trials),length(hilb.delta_phase_trials));

for ii = 1:length(hilb.delta_phase_trials)*length(hilb.delta_phase_trials)
    
    hilb.stats.phase_win_mean_trials_Presound{ii}   = hilb.stats.phase_win_mean_trials{ii}(:,1:hilb.time_zero_idx-1);               
    hilb.stats.phase_win_mean_trials_Sound{ii}      = hilb.stats.phase_win_mean_trials{ii}(:,hilb.time_zero_idx:hilb.time_end_idx);

    hilb.stats.PLV_win_mean_trials_Presound{ii}     = hilb.stats.PLV_win_mean_trials{ii}(:,1:hilb.time_zero_idx-1);                  
    hilb.stats.PLV_win_mean_trials_Sound{ii}        = hilb.stats.PLV_win_mean_trials{ii}(:,hilb.time_zero_idx:hilb.time_end_idx);

end


% Relative Phase and PLV - one mean value for each behavior event
hilb.stats.phase_mean_trials_Presound = cell(length(hilb.delta_phase_trials),length(hilb.delta_phase_trials));
hilb.stats.phase_mean_trials_Sound    = cell(length(hilb.delta_phase_trials),length(hilb.delta_phase_trials));

hilb.stats.PLV_mean_trials_Presound   = cell(length(hilb.delta_phase_trials),length(hilb.delta_phase_trials));
hilb.stats.PLV_mean_trials_Sound      = cell(length(hilb.delta_phase_trials),length(hilb.delta_phase_trials));

for ii = 1:length(hilb.delta_phase_trials)*length(hilb.delta_phase_trials)
    
    hilb.stats.phase_mean_trials_Presound{ii}   = circ_mean(hilb.phase_win_trials_Presound{ii},[],2); 
    hilb.stats.phase_mean_trials_Sound{ii}      = circ_mean(hilb.phase_win_trials_Sound{ii},[],2);   

    hilb.stats.PLV_mean_trials_Presound{ii}    = mean(hilb.PLV_win_trials_Presound{ii},2);  
    hilb.stats.PLV_mean_trials_Sound{ii}       = mean(hilb.PLV_win_trials_Sound{ii},2);     

end


% Relative Phase and PLV - Only one mean value for the whole session
hilb.stats.phase_Total_mean_Presound = zeros(length(hilb.delta_phase_trials),length(hilb.delta_phase_trials));
hilb.stats.phase_Total_mean_Sound    = zeros(length(hilb.delta_phase_trials),length(hilb.delta_phase_trials));

hilb.stats.PLV_Total_mean_Presound   = zeros(length(hilb.delta_phase_trials),length(hilb.delta_phase_trials));
hilb.stats.PLV_Total_mean_Sound      = zeros(length(hilb.delta_phase_trials),length(hilb.delta_phase_trials));

for ii = 1:length(hilb.delta_phase_trials)*length(hilb.delta_phase_trials)
    
    hilb.stats.phase_Total_mean_Presound(ii) = circ_mean(hilb.stats.phase_mean_trials_Presound{ii},[],1); % one value for each session. Average trials during pre sound period
    hilb.stats.phase_Total_mean_Sound(ii)    = circ_mean(hilb.stats.phase_mean_trials_Sound{ii},[],1);    % one value for each session. Average trials during pre sound period
    
    hilb.stats.PLV_Total_mean_Presound(ii)   = mean(hilb.stats.PLV_mean_trials_Presound{ii},1); % one value for each session. Average trials during pre sound period
    hilb.stats.PLV_Total_mean_Sound(ii)      = mean(hilb.stats.PLV_mean_trials_Sound{ii},1);    % one value for each session. Average trials during sound period

end

clear('ii')

%% Plot polar plots - Phase Coherence values - Each trial

figure

% choose par channels to compare
ch = [3 12];

suptitle({'\Delta Phase average over trials. Pre Sound period';['Time window = ' num2str(hilb.time_window) 's ' '- ' 'Overlap = ' num2str(hilb.timeoverlap*100) '%'];[]}) 
set(gcf,'color','white')

sb = 1; %subplot counter

for ii = 1:size(ch,1)
    for jj = 1:parameters.NTrials 
        subplot(2,parameters.NTrials,sb)

        polarplot([zeros(size(hilb.phase_win_trials_Presound{ch(ii,1), ch(ii,2)}(jj,:))), hilb.phase_win_trials_Presound{ch(ii,1), ch(ii,2)}(jj,:)]',repmat([0 1],1,length(hilb.phase_win_trials_Presound{ch(ii,1), ch(ii,2)}(jj,:)))','k');
        hold all
        polarplot([0,hilb.stats.phase_mean_trials_Presound{ch(ii,1), ch(ii,2)}(jj)]',[0 1]','Color','[0.6350, 0.0780, 0.1840]','linew',2);
        
        sb = sb + 1;
    end
end  


for ii = 1:size(ch,1)
    for jj = 1:parameters.NTrials 
        subplot(2,parameters.NTrials,sb)

        polarplot([zeros(size(hilb.phase_win_trials_Sound{ch(ii,1), ch(ii,2)}(jj,:))), hilb.phase_win_trials_Sound{ch(ii,1), ch(ii,2)}(jj,:)]',repmat([0 1],1,length(hilb.phase_win_trials_Sound{ch(ii,1), ch(ii,2)}(jj,:)))','k');
        hold all
        polarplot([0,hilb.stats.phase_mean_trials_Sound{ch(ii,1), ch(ii,2)}(jj)]',[0 1]','Color','[0.6350, 0.0780, 0.1840]','linew',2);
        
        sb = sb + 1;
    end
end

clear ('ii','jj','sb','ch')

%% Plot polar plots - Phase Coherence values - Total over time

figure

% choose par channels to compare
ch = [1 16 ; 3 16 ; 6 16 ; 12 16];

suptitle({'Total \Delta Phase average over time.';['Time window = ' num2str(hilb.time_window) 's ' '- ' 'Overlap = ' num2str(hilb.timeoverlap*100) '%'];[]}) 
set(gcf,'color','white')

for ii = 1:length(ch)

    subplot(2,length(ch),ii)

    polarplot([zeros(size(hilb.stats.phase_win_mean_trials_Presound{ch(ii,1), ch(ii,2)})), hilb.stats.phase_win_mean_trials_Presound{ch(ii,1), ch(ii,2)}]',repmat([0 1],1,length(hilb.stats.phase_win_mean_trials_Presound{ch(ii,1), ch(ii,2)}))','k');
    hold all
    polarplot([0,hilb.stats.phase_Total_mean_Presound(ch(ii,1), ch(ii,2))]',[0 1]','Color','[0.6350, 0.0780, 0.1840]','linew',2);
    
end  

for ii = 1:length(ch)

    subplot(2,length(ch),ii+length(ch))

    polarplot([zeros(size(hilb.stats.phase_win_mean_trials_Sound{ch(ii,1), ch(ii,2)})), hilb.stats.phase_win_mean_trials_Sound{ch(ii,1), ch(ii,2)}]',repmat([0 1],1,length(hilb.stats.phase_win_mean_trials_Sound{ch(ii,1), ch(ii,2)}))','k');
    hold all
    polarplot([0,hilb.stats.phase_Total_mean_Sound(ch(ii,1), ch(ii,2))]',[0 1]','Color','[0.6350, 0.0780, 0.1840]','linew',2);
        
end

clear ('ii','jj','sb','ch')

%% Plot all channels - Phase Coherence values in a color map - Each trial

% Choose channels to plot
chs = [3 6 12 16];

figure
set(gcf,'color','white')

temp1 = zeros(length(hilb.delta_phase_trials));
temp2 = zeros(length(hilb.delta_phase_trials));

for ii = 1:parameters.NTrials 
    for jj = 1:length(hilb.delta_phase_trials)*length(hilb.delta_phase_trials)  
        
        temp1(jj) = hilb.stats.PLV_mean_trials_Presound{jj}(ii);
        temp2(jj) = hilb.stats.PLV_mean_trials_Sound{jj}(ii);

    end 
    
        subplot (2,parameters.NTrials,ii)
        imagesc(chs,chs,temp1(chs,chs))
        xlabel('channels','FontSize',14), ylabel('channels','FontSize',14)
        colorbar
        caxis([0 1])
        
        subplot (2,parameters.NTrials,ii+parameters.NTrials)
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
imagesc(1:length(chs),1:length(chs),hilb.stats.PLV_Total_mean_Presound(chs,chs))
xlabel('channels','FontSize',14), ylabel('channels','FontSize',14)
colorbar
box off
caxis([0 1])
set (gca,'visible','off')

subplot 122
imagesc(1:length(chs),1:length(chs),hilb.stats.PLV_Total_mean_Sound(chs,chs))
xlabel('channels','FontSize',14), ylabel('channels','FontSize',14)
colorbar
box off
caxis([0 1])
set (gca,'visible','off')

clear ('chs')

%%
% save('')

%% last update 26/04/2020 - 01:19
%  listening: Early day miners - Placer found

