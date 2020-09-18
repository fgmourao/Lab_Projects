
%% Phase-amplitude cross-frequency coupling measure

% Quantifying the amount of amplitude modulation by means of anormalized entropy index 
% Tort et al, 2008 -> 10.1073/pnas.0810524105
% Tort et al, 2010 -> 10.1152/jn.00106.2010

% Inputs:
% - Phase = phase time series
% - Amp = amplitude time series
% - position = phase bins (left boundary)

% Outputs:
% - MI = modulation index
% - MeanAmp = amplitude distribution over phase bins (non-normalized)
%             to normalize, do MeanAmp = MeanAmp/sum(MeanAmp)

% Adapted by:
% Flavio Mourao. Nucleo de Neurociencias NNC.
% email: mourao.fg@gmail.com
% Universidade Federal de Minas Gerais
% 05/2020

% [MI,MeanAmp]=ModIndex(Phase, Amp, position)

function [MI,MeanAmp]=ModIndex(Phase, Amp, position)

nbin=length(position);  
winsize = 2*pi/nbin;
 
% now we compute the mean amplitude in each phase:

MeanAmp=zeros(1,nbin); 

for j=1:nbin   
I = find(Phase <  position(j)+winsize & Phase >=  position(j));
MeanAmp(j)= nanmean(Amp(I)); 
end
 
% the center of each bin (for plotting purposes) is position+winsize/2

%  Modulation index
MI=(log(nbin)-(-sum((MeanAmp/sum(MeanAmp)).*log((MeanAmp/sum(MeanAmp))))))/log(nbin);

end
