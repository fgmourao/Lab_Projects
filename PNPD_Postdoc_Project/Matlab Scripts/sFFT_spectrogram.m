
%% Short-time FFT by matlab built function spectrogram

% Flavio Mourao. Nucleo de Neurociencias NNC.
% email: mourao.fg@gmail.com
% Universidade Federal de Minas Gerais
% Started in:  07/2019
% Last update: 04/2020

%%
% First extract the data with: Extracting_LFPs_and_events.m
% ... then organize the data with the script: Pre_processing.m

%%

% Time window
short_fft.timewin    = 5200; % in ms

% Convert time window to points
short_fft.timewinpnts  = round(short_fft.timewin/(1000/parameters.srate));

% Number of overlap samples
short_fft.overlap = 95;
short_fft.noverlap = floor(short_fft.overlap*0.01*short_fft.timewinpnts);

% nFFT
short_fft.nFFT = 2^nextpow2(short_fft.timewinpnts);

% Spectrogram
% lines: frequencies / columns: time / third dimension: channels

for ii = 1:size(data.data{1,1},1)
    if ii ==1
       [short_fft.data(:,:,ii),short_fft.freq,short_fft.time] = spectrogram(data.data{1,1}(ii,:),short_fft.timewinpnts,short_fft.noverlap,short_fft.nFFT,parameters.srate);
    else
        short_fft.data(:,:,ii) = spectrogram(data.data{1,1}(ii,:),short_fft.timewinpnts,short_fft.noverlap,short_fft.nFFT,parameters.srate);
    end
end


clear ('ii','jj')


%% last update 07/04/2020 - 20:43
%  listening: Alice in Chains - Nutshell

