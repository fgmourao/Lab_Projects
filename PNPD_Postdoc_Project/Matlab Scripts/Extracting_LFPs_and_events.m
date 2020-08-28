
function [data, parameters] = Extracting_LFPs_and_events()

% Extracting LFPs and Events from Intan/Open Ephys

% - extract, organize and save data from Intan/Open Ephys:  *.continuous and  *.events
% - The code relies on the following functions : load_open_ephys_data.m (https://github.com/open-ephys/analysis-tools)

% - Option: down sampling data

% - Outputs:

%   "data" 
%   -> data.raw       -> raw data. Original sample rate
%                        Columns: Channels x  Rows:Time
%   -> data.timev_raw -> time vector. Original sample rate

%   -> data.data      -> Cell -> First cell column: signal decimated.
%                        Each cell: Rows: Channels x Columns: Time
%   -> data.timev     -> time vector. Signal decimated

%   -> events         -> External TTls. Events that are detected in a continuous data stream
%                        _ Supports up to 8 inputs. For digital inputs - labels: 0 - 7. 
%                        _ ts -> All timestamps
%                        _ ts_sort .... -> sorted according to the labels
%                                          each cell column corresponds to a recorded event type
%                                          within the cell -> lines -> timestamps(seconds)


%   "parameters"      ->  Record informations and parameters


% by Flavio Mourao. Nucleo de Neurociencias - NNC.
% email: mourao.fg@gmail.com
% Universidade Federal de Minas Gerais. 
% Started in:  11/2019
% Last update: 04/2020

%%
% Load files (*.continuous -> LFP and *.events -> Events)
[FilesLoaded,parameters.Path] = uigetfile({'*.continuous; *.events'},'MultiSelect', 'on'); % Define file type *.*

% Define and save number of channels loaded
parameters.nch = sum(contains(FilesLoaded, 'continuous'));

% Define a struct with files informations from dir organization'
% BEWARE ! This organization changes according to the operating system.
% parameters.FilesLoaded = repmat(struct('name',[],'folder',[],'date',[],'bytes',[],'isdir',[],'datenum',[]), 1, length(FilesLoaded));

% Filename can be organize as a single char or a group char in a cell depending on the number os files selected
if ischar(FilesLoaded)
   parameters.FilesLoaded = dir(fullfile(parameters.Path, FilesLoaded)); % condition for a single file selected       
else    
   for ii = 1:length(FilesLoaded) % loop over multiple files selected
       parameters.FilesLoaded(ii) = dir(fullfile(parameters.Path, char(FilesLoaded(ii))));
   end 
end  

% Optional - Uncomment the line below for sort data. Channels based on a specific file properties. 
% data.Channels = nestedSortStruct(parameters.FilesLoaded,'name',1); % Perform a nested sort of a struct array based on multiple fields. 
                                                                     % >>> https://uk.mathworks.com/matlabcentral/fileexchange/28573-nested-sort-of-structure-arrays?focused=5166120&tab=function

%% Choose factor to LFP down sampling
% - Manually - 
% parameters.downsampling = 6; 


if parameters.nch > 0
    
   % - Request from user -
   prompt        = {'Decimation Factor:'};
   dlgtitle      = 'Please enter';
   dims          = [1 30];
   default_input = {'6'};

   input = inputdlg(prompt,dlgtitle,dims,default_input); %gui

   parameters.downsampling = str2double(input{1,1});

   clear ('prompt','dlgtitle','dims','default_input','input')
   
end

%% Loop to extract data
% Required function: load_open_ephys_data.m
% https://github.com/open-ephys/analysis-tools

for jj = 1:length(parameters.FilesLoaded)
    baseFileName = parameters.FilesLoaded(jj).name;
    fullFileName = fullfile(parameters.Path,baseFileName);
    
    %Identify the file extension
    [~, ~, fExt] = fileparts(baseFileName);
    
    
    switch lower(fExt)
                
        % Case for load channels
        
        case '.continuous'
    
        % Identify the channel number and print the name on the Command Window:
        % channels   1 to 16 and/or 17 to 32

        channel = split(baseFileName,{'100_CH','.continuous'});
        fprintf(1, '\nExtracting LFP from Channel %s\n', channel{2, 1}); 
        
        
        if      jj == 1 && parameters.downsampling == 1
            
                % Load datafiles (*.continuous), timestamps e record info.
                % Raw data - Rows: Time  x Columns: Channels
                % Remove linear trend (Matlab build function 'detrend'/ Method: 'constant' - subtract the mean from the data)
                [data_temp, data.timev_raw, info] = load_open_ephys_data(fullFileName);
                data.raw = detrend(data_temp, 'constant');  % Raw data           
                parameters.header = info.header;            % Data File Header
        
        elseif  jj == 1 && parameters.downsampling > 1
            
                % Load datafiles (*.continuous), timestamps e record info.
                % Raw data - Rows: Time  x Columns: Channels
                % Remove linear trend (Matlab build function 'detrend'/ Method: 'constant' - subtract the mean from the data)
                [data_temp, data.timev_raw, info] = load_open_ephys_data(fullFileName);
                data.raw = detrend(data_temp, 'constant');  % Raw data           
                parameters.header = info.header;   
            
                % Downsampling with Matlab decimate function
                % data - Rows: Channels x Columns: Time
                data(1,1).data{1,1}  = zeros(parameters.nch, ceil(length(data.raw)/parameters.downsampling));
                data.data{1,1}(jj,:) = decimate(data.raw,parameters.downsampling); 
                                           
                % Organize parameters according to the downsampling information
                parameters.srate  = info.header.sampleRate./parameters.downsampling;  % Sampling frequency after downsamplig(Hz)
                parameters.header.downsampling = parameters.downsampling; 
                
                % Normalizing time vector from down sampling data
                data.timev  = (data.timev_raw(1:parameters.downsampling:end)) - min(data.timev_raw);  % Time stamp (sec)
          
        elseif  jj > 1 && parameters.downsampling == 1
            
                % Load datafiles (*.continuous).
                % Rows: Time x Columns: Channels
                % Remove linear trend (function 'detrend')               
                data.raw(:,jj)  = detrend(load_open_ephys_data(fullFileName),'constant'); % Raw data
        
        elseif  jj > 1 && parameters.downsampling > 1
                
                % Load datafiles (*.continuous).
                % Rows: Time x Columns: Channels
                % Remove linear trend (function 'detrend')               
                data.raw(:,jj)  = detrend(load_open_ephys_data(fullFileName),'constant'); % Raw data
            
                % Downsampling with Matlab decimate function
                % data - Rows: Channels x Columns: Time               
                data.data{1,1}(jj,:) = decimate(data.raw(:,jj),parameters.downsampling);            % parameters.downsampling with Matlab decimate function
      
        end

        
        % Case for load events
         
        case '.events'
                
        % Identify TTL events file and print the name on the Command Window:   
        fprintf(1, '\nExtracting %s\n', 'all_channels.events'); 
        
        % Load datafiles (*.continuous), timestamps e record info.
        [data.events.labels, data.events.ts, parameters.events.info] = load_open_ephys_data(fullFileName);
        
        % Sort Events
        % Trigger/events labels according to the digital inputs
        labels = 0:7;

        % data.events.ts_sort -> each cell column corresponds to a recorded event type
        %                        within the cell -> lines -> timestamps(seconds)

        for ii= 1:length(labels)
            data.events.ts_sort{1,ii} = data.events.ts(data.events.labels(:,1) == labels(ii));
        end
        
    end
end  


    
fprintf('\n Done. \n');

end

%% last update 08/04/2020 - 21:19
%  listening: Elliott Smith - Angeles
