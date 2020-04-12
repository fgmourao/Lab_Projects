
%% Plot electrophysiological parameters on the displacement map

% Flavio Mourao. Nucleo de Neurociencias NNC.
% email: mourao.fg@gmail.com
% Universidade Federal de Minas Gerais
% Started in:  04/2020
% Last update: 04/2020

% How to do:

% 1) Identify the beginning and the end of trials through the Track Video GUIDE
% 2) Extract the X and Y coordinates acquired in the Bonsai
% 3) Choose a channel of interest and run the analyses.
% 4) Normalize vectors. Usually the video and the record start separately. First the record then the video
% 5) Perform signal interpolation in order to normalize the number of points
% 6 To plot:
%   - scatter3 -> Simplest option
%   - cline    -> Third party function. Similar to plot3 but attaches color map
%                 https://la.mathworks.com/matlabcentral/fileexchange/14677-cline

%% Beginning and the end of trials

% Put all in one struct variable
behavior.TS_LFPindex = TS_LFPindex;
behavior.TS_LFPsec   = TS_LFPsec;
behavior.TSframesc   = TSframes;
behavior.TSseconds   = TSseconds;

% Correction factor to normalize time vector
behavior.time_factor_correction = mean(mean(behavior.TS_LFPsec-behavior.TSseconds));

clear('TS_LFPindex','TS_LFPsec','TSframes','TSseconds')

%% X and Y coordinates

% Put all in one struct variable
behavior.data = data;
clear('data', 'Header')

%% sFFT from spectrogram 

% Frequency steps according to the fft time window
steps       = diff(short_fft.freq); 

% Choose a frequency band to plot
freq2plot   = 1:steps:5;    
closestfreq = dsearchn(short_fft.freq,freq2plot');

% Perform the average of frequency band
behavior.sfft     = mean(abs(short_fft.data(closestfreq,:)),1);

% Time vector
behavior.timesfft = short_fft.time;

clear('steps','freq2plot','closestfreq')

%%  Normalize vectors.

% Remove the time spoch which is not matched
to_remove = dsearchn(behavior.timesfft',behavior.time_factor_correction );

behavior.sfft = behavior.sfft(to_remove:end);
behavior.timesfft = behavior.timesfft(to_remove:end);

clear ('to_remove')

%%  Main Time Vector

% Record time / number of frames
behavior.timev = linspace(1,short_fft.time(end),length(behavior.data));

%% signal interpolation

behavior.sfft_int = interp1((1:numel(behavior.sfft)), behavior.sfft, linspace(1, numel(behavior.sfft), numel(behavior.timev)), 'linear')';

%% Plot

figure
set(gcf,'color','white')

% First Option
% scatter3(behavior.data(:,1),behavior.data(:,2),behavior.sfft_int(:),30,behavior.sfft_int(:),'filled')

% Second Option
% https://la.mathworks.com/matlabcentral/fileexchange/14677-cline

c = cline(behavior.data(:,1),behavior.data(:,2),behavior.sfft_int(:),behavior.sfft_int(:));
c.LineWidth = 6;
axis off
colorbar

clear('c')

%% last update 11/04/2020 - 21:20am
%  listening: Milton Nascimento - Fé cega / faca amolada

