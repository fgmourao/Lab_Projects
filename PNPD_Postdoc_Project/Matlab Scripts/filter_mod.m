function [data_filtered] = filter_mod(data,f_cut,srate)
% Filter using the matlab filt filt function.

% Parameters define by hand.

%  Works very well to filter the modulated envelope
%  Building the filter by hand avoided deformations to 53.71 modulated frequency. 
%  Furthermore, it guarantees smooth edges in the CS modulating transitions.

% - Inputs:
%   data
%   f_cut : bandscutoff [start end]
%   srate : sampling frequency

% - Outputs:
%   data_filtered
%   parameters : used parameters

% by Flavio Mourao. Nucleo de Neurociencias - NNC.
% email: mourao.fg@gmail.com
% Universidade Federal de Minas Gerais. 
% Started in:  07/2019
% Last update: 05/2020


% Specify nyquist frequency
nyquistS = srate/2;

% filter frequency band
filter.filtbound(1,:) = f_cut; % Hz

% transition width
filter.trans_width(1,:) = 0.20; % fraction of 1, thus 20%

% filter order
filter.filt_order(1,:) = round(3*(srate/filter.filtbound(1,1)));

% frequency vector (as fraction of nyquist)
filter.ffrequencies(1,:)  = [ 0 (1-filter.trans_width(1,:))*filter.filtbound(1,1) filter.filtbound(1,:) (1+filter.trans_width(1,:))*filter.filtbound(1,2) nyquistS ]/nyquistS;

% shape of filter (must be the same number of elements as frequency vector
filter.idealresponse = [ 0 0 1 1 0 0 ];

% get filter weights
filter.filterweights = cell(size(filter.filtbound,1),1);
filter.filterweights{1} = firls(filter.filt_order(1,:),filter.ffrequencies(1,:),filter.idealresponse);


% Apply filter to the data
data_filtered = filtfilt(filter.filterweights{1},1,double(data));


%% Plot filter parameters for visual inspection

% figure, clf
% set(gcf,'color','white')
% 
% subplot(1,2,1)
% plot(filter.ffrequencies(1,:)*nyquistS,filter.idealresponse,'k--o','markerface','[0.6350, 0.0780, 0.1840]')
% set(gca,'ylim',[-.1 1.1],'xlim',[-2 nyquistS+2])
% xlabel('Frequencies (Hz)'), ylabel('Response amplitude')
% legend({['Bandcuts',' ',num2str(filter.filtbound(1,1)),'-',num2str(filter.filtbound(1,2))]},'Location','northeast')
% legend('boxoff')
% 
% subplot(1,2,2)
% plot((0:filter.filt_order(1,:))*(1000/srate),filter.filterweights{1},'Color','[0.6350, 0.0780, 0.1840]','linew',2)
% xlabel('Time (ms)'), ylabel('Amplitude')
% legend({['filter Weights',' ',num2str(filter.filtbound(1,1)),'-',num2str(filter.filtbound(1,2))]},'Location','northeast')
% legend('boxoff')

end
