
%% Phase-amplitude Cross-frequency coupling measure
%  Considering the behavior events

% - Performs analysis with raw and surrogate values

% The code relies on the following functions:
% --> ModIndex.m    - by Adriano Tort, Instituto do Cerebro - Universidade Federal do Rio Grande do Norte
% --> shuffle_esc.m - by Rafal Bogacz, Angela Onslow, May 2010

% See Tort et al, 2010 -> 10.1152/jn.00106.2010 
%     Tort et al, 2008 -> 10.1073/pnas.0810524105


% Flavio Mourao. Nucleo de Neurociencias NNC.
% email: mourao.fg@gmail.com
% Universidade Federal de Minas Gerais
% Started in:  05/2020
% Last update: 05/2020

%%
% First extract the data with: Extracting_LFPs_and_events.m
% ... then organize the data with the script: Pre_processing.m

%% Run each session sequentially

%% Hilbert Transform

% in each cell -> frequencies cutoff
% - rows          - > channels
% - columns       - > time
% - 3nd dimension - > behavioral events

% Choose frequency bands index according to the Pre_processing.m

% 2  - modulator / 3 - delta / 4 - lowtheta / 5 - hightheta
% 6  - alpha     / 7 - beta  / 8 - lowgamma / 9 - highgamma
% 10 - extra

% 1st index - > phase
% 2nd index - > amplitude

f_idx = [4 2];

% Loop over channels and make hilbert transform
% Initializing with NaN
MI.behavior.hilb_coefficients = cell(size(f_idx));
MI.behavior.hilb_coefficients(:,:) = {nan(size(data.data_behavior{1,1}))};

for ii = 1:length(f_idx)
    for jj = 1:size(data.data_behavior{1, f_idx(ii)},3)
        for ll = 1:size(data.data_behavior{1, f_idx(ii)},1)
        
            temp = data.data_behavior{1, f_idx(ii)}(ll,:,jj);
            temp(isnan(temp(:,:)))=[]; % it was necessary to remove the nan for the function to work
            MI.behavior.hilb_coefficients{1, ii}(ll,1:length(temp),jj) = hilbert(temp);
    
        end 
    end 
end

% Extract amplitude and phase

% MI.behavior.data - colunm 1: phase / colunm 2: amplitude

MI.behavior.data{1,1} = angle(MI.behavior.hilb_coefficients{1,1});
MI.behavior.data{1,2} = abs(MI.behavior.hilb_coefficients{1,2});

% Time vector
MI.behavior.time = (linspace(-parameters.behavior.Tpre,parameters.behavior.trialperiod+parameters.behavior.Tpos,size(MI.behavior.data{1, 1},2)));

clear('ii','jj','ll','temp','f_idx')

%% Cut Time periods - pre behavior and during behavior events. Ignore "pos" period

% Time index -->  behavior epoch begins / ends

MI.behavior.time_zero_idx = dsearchn(MI.behavior.time',0); % time zero index. Behavior Start.

for ii = 1:length(data.events.behavior.TS_LFPsec)

    temp  = data.events.behavior.TS_LFPsec(ii,2) - data.events.behavior.TS_LFPsec(ii,1); % time end index.  Behavior End.
    MI.behavior.time_end_idx(ii) = dsearchn(MI.behavior.time',temp');

end


% from...all behavior epochs over time
% Initialize variables
MI.behavior.data_before      = cell(1,length(MI.behavior.data));
MI.behavior.data_before(:,:) = {nan(parameters.nch+1,max(MI.behavior.time_zero_idx)-1,parameters.behavior.NTrials)};

MI.behavior.data_during      = cell(1,length(MI.behavior.data));
MI.behavior.data_during(:,:) = {nan(parameters.nch+1,max(MI.behavior.time_end_idx - MI.behavior.time_zero_idx)+1,parameters.behavior.NTrials)};

for ii = 1:length(MI.behavior.data)
    for jj = 1:parameters.behavior.NTrials
        
        MI.behavior.data_before{1,ii}(:,1:(MI.behavior.time_zero_idx-1),jj) = MI.behavior.data{1,ii}(:,1:MI.behavior.time_zero_idx-1,jj);              
        MI.behavior.data_during{1,ii}(:,1:(MI.behavior.time_end_idx(jj) - MI.behavior.time_zero_idx)+1,jj) = MI.behavior.data{1,ii}(:,MI.behavior.time_zero_idx:MI.behavior.time_end_idx(jj),jj); 

    end
end


clear( 'ii','jj','temp')

%% - CHANNELS MAP - 
% * just in case to check

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

%% Phase-amplitude cross-frequency coupling measure

% MI_values 
% - rows          - > events
% - 1st column    - > before onset
% - 2nd column    - > post onset
% - 3nd dimension - > rearrangements

% MeanAmp 
% - rows          - > events
% - coluns        - > amplitude distribution over phase bins
% - 3nd dimension - > 1st: before onset
%                     2nd: post onset
% - 3nd dimension - > rearrangements

% Choose channels
MI.behavior.params.ch_phase = 11;
MI.behavior.params.ch_amp   = 16;

% Define time to analyse
MI.behavior.params.time2analise = 5;
time2samples = MI.behavior.params.time2analise * parameters.srate; % in samples

% Define number number of phase bins 
MI.behavior.params.nbin = 20; 

% variable not centered in the phase bin (rad)
MI.behavior.params.position=zeros(1,MI.behavior.params.nbin); 

MI.behavior.params.winsize = 2*pi/MI.behavior.params.nbin;

for jj = 1:MI.behavior.params.nbin
    MI.behavior.params.position(jj) = -pi+(jj-1)*MI.behavior.params.winsize;
end

% Loop over behavior events
MI.behavior.MI_value = nan(parameters.behavior.NTrials,2);
MI.behavior.MeanAmp  = nan(parameters.behavior.NTrials,MI.behavior.params.nbin,2);

for ii = 1:parameters.behavior.NTrials

[MI.behavior.MI_value(ii,1),MI.behavior.MeanAmp(ii,:,1)] = ModIndex(MI.behavior.data_before{1,1}(MI.behavior.params.ch_phase,:,ii), MI.behavior.data_before{1,2}(MI.behavior.params.ch_amp,:,ii), MI.behavior.params.position);
[MI.behavior.MI_value(ii,2),MI.behavior.MeanAmp(ii,:,2)] = ModIndex(MI.behavior.data_during{1,1}(MI.behavior.params.ch_phase,1:time2samples,ii), MI.behavior.data_during{1,2}(MI.behavior.params.ch_amp,1:time2samples,ii), MI.behavior.params.position);

end

% Stats

% mean events
MI.behavior.stats.MI_mean = mean(MI.behavior.MI_value);

clear('ii','jj','time2samples')

%% Plot MI values

figure
set(gcf,'color','w');

suptitle({' ';'\fontsize{14}Phase-amplitude cross-frequency coupling during behavior events'; 
    ['Phase channel: ' num2str(MI.behavior.params.ch_phase),' / Amplitude channel: ' num2str(MI.behavior.params.ch_amp)];' ' });

subplot 121
b1 = bar(MI.behavior.MI_value);
b1(1).FaceColor = 'w'; b1(2).FaceColor = [0.8, 0.8, 0.8];

for ii = 1:parameters.behavior.NTrials
    labels(ii)={['Event ' num2str(ii)]};
end

set(gca,'xticklabel',labels(1:1:end));
ylabel('\fontsize{14} Modulation Index');

legend({'Before';'After'},'Box','off');
axis square
box off

subplot 122
b2 = bar(MI.behavior.stats.MI_mean);
axis square

title('Mean Events')
b2.FaceColor = 'flat';
b2.CData(1,:) = 'w' ; b2.CData(2,:) = [0.8, 0.8, 0.8];

labels={'Before';'After'};
set(gca,'xticklabel',labels(1:1:end));
ylabel('\fontsize{14} Modulation Index');

box off

clear('b1','b2','labels','ii')

%% Plot phase x amplitude distributions

figure
set(gcf,'color','w');

xvalue1 = rad2deg(MI.behavior.params.position) + 180;
xvalue2 = xvalue1 + 360;

count = 1;
for ii = 1:2:parameters.behavior.NTrials*2
    
    subplot(parameters.behavior.NTrials,2,ii)
    b3 = bar([xvalue1 xvalue2],[MI.behavior.MeanAmp(count,:,1)./sum(MI.behavior.MeanAmp(count,:,1)) MI.behavior.MeanAmp(count,:,1)./sum(MI.behavior.MeanAmp(count,:,1))]);
    b3.FaceColor = 'w';
    
    if count == 1
       title({'pre-onset';'   '})    
    end
    
    ylabel('\fontsize{11} amplitude');
    
    box off
    count = count + 1;
    %ylim([0 60])
end

    xlabel('\fontsize{11} phase');

count = 1;
for ii = 2:2:parameters.behavior.NTrials*2
    
    subplot(parameters.behavior.NTrials,2,ii)
    b4 = bar([xvalue1 xvalue2],[MI.behavior.MeanAmp(count,:,2)./sum(MI.behavior.MeanAmp(count,:,2)) MI.behavior.MeanAmp(count,:,2)./sum(MI.behavior.MeanAmp(count,:,2))]);
    b4.FaceColor = [0.8, 0.8, 0.8];
    
    if count == 1
       title({'post-onset';'  '})    
    end

    box off
    count = count + 1;
    %ylim([0 60])
end

xlabel('\fontsize{11} phase');
    
clear('xvalue1','xvalue2','ii','b3','b4', 'count')

%% Surrogate phase vectors to compare

% MI_shuffle_values 
% - rows          - > events
% - 1st column    - > before onset
% - 2nd column    - > post onset
% - 3nd dimension - > rearrangements

% MeanAmp_shuffle 
% - rows          - > events
% - coluns        - > amplitude distribution over phase bins
% - 3nd dimension - > 1st: before onset
%                     2nd: post onset
% - 3nd dimension - > rearrangements

numshf  = 200; % number of shuffled segments
nsurrog = 200; % number of rearrangements

time2samples = MI.behavior.params.time2analise * parameters.srate; % in samples

% Loop over behavior events 
MI.behavior.MI_shuffle_values = nan(parameters.behavior.NTrials,2,nsurrog);
MI.behavior.MeanAmp_shuffle   = nan(parameters.behavior.NTrials,MI.behavior.params.nbin,2,nsurrog);


before_shuffle = [];
during_shuffle = [];

for jj = 1:nsurrog
    
    for ii = 1:parameters.behavior.NTrials
    
        before_shuffle(1,:,ii) = shuffle_esc(MI.behavior.data_before{1,1}(MI.behavior.params.ch_phase,:,ii),parameters.srate,numshf);
        during_shuffle(1,:,ii) = shuffle_esc(MI.behavior.data_during{1,1}(MI.behavior.params.ch_phase,1:time2samples,ii),parameters.srate,numshf);
   
    end
    
    for ii = 1:parameters.behavior.NTrials
    
        [MI.behavior.MI_shuffle_values(ii,1,jj),MI.behavior.MeanAmp_shuffle(ii,:,1,jj)] = ModIndex(before_shuffle(1,:,ii), MI.behavior.data_before{1,2}(MI.behavior.params.ch_amp,:,ii), MI.behavior.params.position);
        [MI.behavior.MI_shuffle_values(ii,2,jj),MI.behavior.MeanAmp_shuffle(ii,:,2,jj)] = ModIndex(during_shuffle(1,1:time2samples,ii), MI.behavior.data_during{1,2}(MI.behavior.params.ch_amp,1:time2samples,ii), MI.behavior.params.position);
    
    end
    
    before_shuffle = [];
    during_shuffle = [];
    
end

% z surrogated values
MI.behavior.stats.z_MI_shuffle_values = zscore(MI.behavior.MI_shuffle_values);

% z real(observed) values
MI.behavior.stats.z_MI_value = (MI.behavior.MI_value - (mean(MI.behavior.MI_shuffle_values,3)))./std(MI.behavior.MI_shuffle_values,[],3);

% p values
MI.behavior.stats.p_MI_value = 0.5 * erfc(MI.behavior.stats.z_MI_value ./ sqrt(2));


%clear('jj','ii','nsurrog','numshf','time2samples','before_shuffle','during_shuffle')

%% Plot statistical distribuitions

figure
set(gcf,'color','w');

MI.behavior.params.nbins = 30;

count = 1;
for ii = 1:2:parameters.behavior.NTrials*2
    
    subplot(parameters.behavior.NTrials,2,ii)
    
    h = histogram(MI.behavior.stats.z_MI_shuffle_values(count,1,:),MI.behavior.params.nbins);
    h.FaceColor = 'w';
    hold on
    plot([MI.behavior.stats.z_MI_value(count,1) MI.behavior.stats.z_MI_value(count,1)],[0 30],'Color', '[0.6350, 0.0780, 0.1840]','linewidth',4);
    %ylim ([0 30])
    
    % just to make arrow
    %[figx,figy] = dsxy2figxy([MI.behavior.stats.z_MI_value(ii,1) MI.behavior.stats.z_MI_value(ii,1)],[10 0]); % Transform point or MI.behavior.params.position from data space 
                                                                                            % coordinates into normalized figure coordinates . 
                                                                                            % >>> https://uk.mathworks.com/matlabcentral/fileexchange/30347-sigplot?focused=5178148&tab=function
    %annotation('textarrow',figx,figy,'Color', '[0.6350, 0.0780, 0.1840]','linewidth',4);
    
    if count == 1
       title({'pre-onset';'  ';['MI = ' num2str(MI.behavior.MI_value(count,1)*10^3),'  ','\rm\fontsize{12}p_v_a_l_u_e = ' num2str(MI.behavior.stats.p_MI_value(count,2))]})    
    
    else
       title(['MI = ' num2str(MI.behavior.MI_value(count,1)*10^3),'  ','\rm\fontsize{12}p_v_a_l_u_e = ' num2str(MI.behavior.stats.p_MI_value(count,1)) ])
    
    end
    
    ylabel('\fontsize{14}Frequency')
    %title(['MI = ' num2str(MI.behavior.MI_value(count,1)),'  ','\rm\fontsize{12}p_v_a_l_u_e = ' num2str(MI.behavior.stats.p_MI_value(count,1)) ])
    %legend({'\fontsize{12}Random';'\fontsize{12}Observed'},'box','off');

    box off

    count = count + 1;

end

xlabel('\fontsize{12}z values')

count = 1;
for ii = 2:2:parameters.behavior.NTrials*2
    
    subplot(parameters.behavior.NTrials,2,ii)
    
    h = histogram(MI.behavior.stats.z_MI_shuffle_values(count,2,:),MI.behavior.params.nbins);
    h.FaceColor = [0.8, 0.8, 0.8];
    hold on
    plot([MI.behavior.stats.z_MI_value(count,2) MI.behavior.stats.z_MI_value(count,2)],[0 30],'Color', '[0.6350, 0.0780, 0.1840]','linewidth',4);
    %ylim ([0 30])

    % just to make arrow
    %[figx,figy] = dsxy2figxy([MI.behavior.stats.z_MI_value(ii,1) MI.behavior.stats.z_MI_value(ii,1)],[10 0]); % Transform point or MI.behavior.params.position from data space 
                                                                                            % coordinates into normalized figure coordinates . 
                                                                                            % >>> https://uk.mathworks.com/matlabcentral/fileexchange/30347-sigplot?focused=5178148&tab=function
    % annotation('textarrow',figx,figy,'Color', '[0.6350, 0.0780, 0.1840]','linewidth',4);
    
    if count == 1
       title({'post-onset';'  ';['MI = ' num2str(MI.behavior.MI_value(count,2)*10^3),'  ','\rm\fontsize{12}p_v_a_l_u_e = ' num2str(MI.behavior.stats.p_MI_value(count,2))]})    
    else
       title(['MI = ' num2str(MI.behavior.MI_value(count,2)*10^3),'  ','\rm\fontsize{12}p_v_a_l_u_e = ' num2str(MI.behavior.stats.p_MI_value(count,2)) ])
    end
    
    %title(['MI = ' num2str(MI.behavior.MI_value(count,2)),'  ','\rm\fontsize{12}p_v_a_l_u_e = ' num2str(MI.behavior.stats.p_MI_value(count,2)) ])
    %legend({'\fontsize{12}Random';'\fontsize{12}Observed'},'box','off');

    box off

    count = count + 1;

end

xlabel('\fontsize{12}z values')
legend({'\fontsize{12} Distribuition';'\fontsize{12}Observed'},'box','off')%,'location','bestoutside')

clear ('ii','jj','count','h','figx','figy');

%% last update 14/05/2020 - 11:37am
%  listening: The Smiths - Stretch Out and Wait
