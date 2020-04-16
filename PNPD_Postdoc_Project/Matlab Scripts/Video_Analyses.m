
%%  Behavioral Analysis. 

% Files need:
% - *.csv files imported from Bonsai software
% - frame trials index identified through the Video Guide

% Common *.csv files organization for full arena analyse
%        column 1 -> x coordinates (pixels)
%        column 2 -> y coordinates (pixels)
%        column 3 -> arena width   (pixels)
%        column 4 -> arena height  (pixels)

% - Analyse one or more *.csv files

% - Analyse the whole arena or divided into quadrant according to the bonsai output

% Outputs:

% Variables : Cell -> Analyze.Arena_full 
%             each cell corresponds to onde *.csv file

%   - x and y coordinates in pixels              -> vec_y/x
%   - x and y coordinates in cm                  -> vec_y/x_cm
%   - x / y derivative in cm                     -> d_vec_y/x_cm
%   - Displacement in cm                         -> displacement 
%   - Accumulated distance over time in cm       -> Accumulate_distance
%   - Total distance covered in cm               -> Total_distance
%   - Time vector                                -> Time_vector

%   - Total distance covered in cm in each trial (sound on)  -> Total_distance_Trials
%   - Total distance covered in all trials (sound on)        -> Total_distance_ALL_Trials
%     usually : column 1 pre-conditioning
%               column 2 test
%               Row 1 = in cm
%               Row 2 = in % of total distance
%   - Total distance covered out of trials (sound off)       -> Total_distance_ALL_Trials
%     usually : column 1 pre-conditioning
%               column 2 test
%               Row 1 = in cm
%               Row 2 = in % of total distance


% PLots
% - Tracking map (pivels or cm)
% - Displacement (cm)
% - Accumulated distance over time
% - Trials period  with different colors


% by Flavio Mourao. Nucleo de Neurociencias - NNC.
% email: mourao.fg@gmail.com
% Universidade Federal de Minas Gerais. 
% Started in:  05/2019
% Last update: 04/2020

%%
% Load datafiles (*.csv)
[FileName,PathName] = uigetfile({'*.csv'},'MultiSelect', 'on'); % Define file type *.*

% Filename can be organize as a single char or a group char in a cell depending on the number os files selected

if ischar(FileName)
   Header.FilePattern = dir(fullfile(PathName, FileName)); % condition for a single file selected 
   fprintf(1, 'Reading %s\n', FileName);
   data = readtable([PathName '/' FileName],'Delimiter',',');
   data = table2array(data);
%    val_name = FileName(1:regexp(FileName,'.csv')-1);
%    assignin('base',val_name, data);
   
else   
    
   data = cell(size(FileName));

   for ii = 1:length(FileName) % loop over multiple files selected
       Header.FilePattern(ii) = dir(fullfile(PathName, char(FileName(ii))));
       fprintf(1, 'Reading %s\n', FileName{ii});
       output = readtable([PathName '/'  FileName{ii}],'Delimiter',',');
       data{ii} = table2array(output);
%        val_name = FileName{ii};
%        val_name = val_name(1:regexp(val_name,'.csv')-1);
%        assignin('base',val_name, data);
   end 
end 

clear ('FileName', 'ii', 'val_name','PathName','output');

% Sort data based on file properties
% BehavFiles = nestedSortStruct(Header.FilePattern,'name',1); % Perform a nested sort of a struct array based on multiple fields.

%% Beginning and the end of trials

% load('')

% usually
% column 1 pre-conditioning
% column 2 test

% Put all in one struct variable. Change columns manually
Analyze.TSframes{1,2}   = TSframes;
Analyze.TSseconds{1,2}   = TSseconds;

clear ('TSframes', 'TSseconds','TS_LFPindex','TS_LFPsec')

%% Analysis

% Set the video parameters

% Frame Rate
Header.Num_frames = 30; 

% Pixels area (check the bonsai values defined)
Header.Video_width  = 461; % data{1, 1}(1,3);
Header.Video_height = 454; % data{1, 1}(1,4);

% Arena Dimensions
Header.Arena_width  = 30; % in cm
Header.Arena_height = 30; % in cm

% Define conversion Analyze.factor (cm from pixels)
Analyze.factor_w = Header.Arena_width/Header.Video_width;   % width
Analyze.factor_h = Header.Arena_height/Header.Video_height; % height


%% Analyse considering the entire arena area

% Extracting parameters

for ii = 1:length(data)

Analyze.Arena_Full.vec_y{ii} = data{ii}(:,2); % y axis - in pixels
Analyze.Arena_Full.vec_x{ii} = data{ii}(:,1); % x axis - in pixels

Analyze.Arena_Full.vec_y_cm{ii} = data{ii}(:,2).* Analyze.factor_h; % y axis - convert y values to cm from pixels
Analyze.Arena_Full.vec_x_cm{ii} = data{ii}(:,1).* Analyze.factor_w; % x axis - convert x values to cm from pixels

Analyze.Arena_Full.d_vec_y_cm{ii} = [0 ; diff(data{ii}(:,2))].* Analyze.factor_h; % convert y Displacement to cm from pixels
Analyze.Arena_Full.d_vec_x_cm{ii} = [0 ; diff(data{ii}(:,1))].* Analyze.factor_w; % convert x Displacement to cm from pixels

Analyze.Arena_Full.Displacement{ii}        = sqrt(Analyze.Arena_Full.d_vec_x_cm{ii} .^2 + Analyze.Arena_Full.d_vec_y_cm{ii} .^2); % Displacement in cm
Analyze.Arena_Full.Accumulate_distance{ii} = cumsum(Analyze.Arena_Full.Displacement{ii});                                         % Accumulate distance in cm
Analyze.Arena_Full.Total_distance(ii)      = Analyze.Arena_Full.Accumulate_distance{ii}(end);                                     % Total distance in cm
Analyze.Arena_Full.Time_vector{ii}         = linspace(0,length(data{ii})/Header.Num_frames,length(data{ii}));                     % Time vector in sec.

    for jj = 1:size(Analyze.TSframes{ii},1)
        Analyze.Arena_Full.Total_distance_Trials{ii}(jj,1) = Analyze.Arena_Full.Accumulate_distance{ii}(Analyze.TSframes{ii}(jj,2)) - Analyze.Arena_Full.Accumulate_distance{ii}(Analyze.TSframes{ii}(jj,1)); % Total distance in cm per trial (sound on)
    end 
    
Analyze.Arena_Full.Total_distance_ALL_Trials(1,ii) = sum(Analyze.Arena_Full.Total_distance_Trials{ii});                                        % Total distance in cm in all trials - first row in cm (sound on)
Analyze.Arena_Full.Total_distance_ALL_Trials(2,ii) = Analyze.Arena_Full.Total_distance_ALL_Trials(1,ii)/Analyze.Arena_Full.Total_distance(ii); % Total distance in cm in all trials - second row in % (sound on)

Analyze.Arena_Full.Total_distance_OUT_Trials(1,ii) = Analyze.Arena_Full.Total_distance(ii) - Analyze.Arena_Full.Total_distance_ALL_Trials(1,ii); % Total distance in cm out of trials - first row in cm (sound off)
Analyze.Arena_Full.Total_distance_OUT_Trials(2,ii) = Analyze.Arena_Full.Total_distance_OUT_Trials(1,ii)/Analyze.Arena_Full.Total_distance(ii);   % Total distance in cm out of trials - second row in % (sound off)

end



clear ('ii', 'jj')

%% Plots

figure;

set(gcf,'color','w');

titles = {'Pre CS';'Test'};

for ii = 1:length(data)
    subplot(3,2,ii)
    hold on
    
    plot(Analyze.Arena_Full.vec_x_cm{ii},Analyze.Arena_Full.vec_y_cm{ii},'Color',[0.6, 0.6, 0.6],'linewidth',2)
    title(titles{ii},'FontSize',14);
    
    for jj = 1:size(Analyze.TSframes{ii},1)
        plot(Analyze.Arena_Full.vec_x_cm{ii}(Analyze.TSframes{ii}(jj,1):Analyze.TSframes{ii}(jj,2)),Analyze.Arena_Full.vec_y_cm{ii}(Analyze.TSframes{ii}(jj,1):Analyze.TSframes{ii}(jj,2)),'Color','[0.6350, 0.0780, 0.1840]','linewidth',2)
                
    end
    
     xlim ([0 Header.Arena_width]);  % in cm
     ylim ([0 Header.Arena_height]); % in cm

    axis off
    
end

for ii = 1:length(data)
subplot(3,2,ii+2)
hold on
plot(Analyze.Arena_Full.Time_vector{ii},Analyze.Arena_Full.Displacement{ii}','Color',[0.6, 0.6, 0.6],'linewidth',2)

    for jj = 1:size(Analyze.TSframes{ii},1)
        plot(Analyze.Arena_Full.Time_vector{ii}(Analyze.TSframes{ii}(jj,1):Analyze.TSframes{ii}(jj,2)),Analyze.Arena_Full.Displacement{ii}(Analyze.TSframes{ii}(jj,1):Analyze.TSframes{ii}(jj,2)),'Color','[0.6350, 0.0780, 0.1840]','linewidth',2)
        
%         plot([Analyze.TSseconds(jj,1) Analyze.TSseconds(jj,1)],[0 max(Analyze.Arena_Full.Displacement{1})/2],'--','Color','[0.6350, 0.0780, 0.1840]','linewidth',1)
%         plot([Analyze.TSseconds(jj,2) Analyze.TSseconds(jj,2)],[0 max(Analyze.Arena_Full.Displacement{1})/2],'--','Color','[0.6350, 0.0780, 0.1840]','linewidth',1)               
    end

xlabel('Time (s)')
ylabel('Displacement (cm)')

xlim ([0 Analyze.Arena_Full.Time_vector{ii}(end)]); 
ylim ([0 max(Analyze.Arena_Full.Displacement{1})]);

box off

end

for ii = 1:length(data)
subplot(3,2,ii+4)
hold on
plot(Analyze.Arena_Full.Time_vector{ii},Analyze.Arena_Full.Accumulate_distance{ii},'Color',[0.6, 0.6, 0.6],'linewidth',3)
    
    for jj = 1:size(Analyze.TSframes{ii},1)
        plot(Analyze.Arena_Full.Time_vector{ii}(Analyze.TSframes{ii}(jj,1):Analyze.TSframes{ii}(jj,2)),Analyze.Arena_Full.Accumulate_distance{ii}(Analyze.TSframes{ii}(jj,1):Analyze.TSframes{ii}(jj,2)),'Color','[0.6350, 0.0780, 0.1840]','linewidth',4)

%         plot([Analyze.TSseconds(jj,1) Analyze.TSseconds(jj,1)],[0 max(Analyze.Arena_Full.Accumulate_distance{1})],'--','Color','[0.6350, 0.0780, 0.1840]','linewidth',1)
%         plot([Analyze.TSseconds(jj,2) Analyze.TSseconds(jj,2)],[0 max(Analyze.Arena_Full.Accumulate_distance{1})],'--','Color','[0.6350, 0.0780, 0.1840]','linewidth',1)               
    end

xlim ([0 Analyze.Arena_Full.Time_vector{ii}(end)]); 
ylim ([0 max(Analyze.Arena_Full.Accumulate_distance{1})]);

xlabel('Time (s)')
ylabel('Accumulated distance (cm)')

legend({'sound off', 'sound on'},'FontSize',12, 'location', 'southeast')
legend('boxoff')

box off

end

clear ('titles','ii','jj')

%% Analysis by quadrants

% THIS SESSION NEEDS TO BE UPDATE

% Number of quadrants
% n_q = 2;
% 
% % Analysis of each quadrant
% sq = (M1(:,4:5)); 
% z = zeros(1,n_q);
% idx = false(size(sq));
% 
% % Pre allocate for data analysis for each quadrant
% Arena_Quadrants = cell(5,4);
% 
% for ii=1:n_q
%     z(1,ii) = ii;
%     idx = ismember(sq,z,'rows');
%     
%     Arena_Quadrants{1,ii} = M1(idx,1);                                                                                    % Frames
%     Arena_Quadrants{2,ii} = [Arena_Full.vec_x(idx,1) Arena_Full.vec_y(idx,1)];                                            % x y coordinates in cm
%     Arena_Quadrants{3,ii} = [Arena_Full.d_vec_x(idx,1) Arena_Full.d_vec_y(idx,1)];                                        % x y Displacement in cm
%     Arena_Quadrants{4,ii} = cumsum(Arena_Full.Displacement(idx,1));                                                       % Accumulate distance in cm
%     Arena_Quadrants{5,ii} = linspace(0,length(Arena_Quadrants{1,ii})/Header.Num_frames,length(Arena_Quadrants{1,ii}))';   % Time vector
%     Arena_Quadrants{6,ii} = Arena_Quadrants{4,ii}./ Arena_Quadrants{5,ii} ;                                               % Velocity 
%     Arena_Quadrants{7,ii} = round(length(Arena_Quadrants{1,ii})./Header.Num_frames);                                      % Total time
%     Arena_Quadrants{8,ii} = Arena_Quadrants{4,ii}(end);                                                                   % Total distance traveled 
%     Arena_Quadrants{9,ii} = sum(diff(Arena_Quadrants{1,ii}) > 1);                                                         % Crossings
%     
%     z(1,ii) = 0;
% end
% 
% clear ('Analyze.factor_h', 'Analyze.factor_w', 'idx', 'ii', 'n_q', 'sq', 'z')
% 
% %
% 
% Track_Fig = figure('position', [0, 0, 600, 400],'resize', 'on');
% set(gcf,'color','w');
% hold
% plot(Arena_Quadrants{2,1}(:,1),Arena_Quadrants{2,1}(:,2),'Color',[0.6, 0.6, 0.6],'linewidth',2)
% plot(Arena_Quadrants{2,2}(:,1),Arena_Quadrants{2,2}(:,2),'Color',[0.6, 0.6, 0.6],'linewidth',2)
% plot(Arena_Quadrants{2,3}(:,1),Arena_Quadrants{2,3}(:,2),'Color',[0.6, 0.6, 0.6],'linewidth',2)
% plot(Arena_Quadrants{2,4}(:,1),Arena_Quadrants{2,4}(:,2),'Color',[0.6, 0.6, 0.6],'linewidth',2)
% 
% subplot 221
% plot(Arena_Quadrants{5,1},Arena_Quadrants{6,1},'Color',[0.6, 0.6, 0.6],'linewidth',2)
% subplot 222
% plot(Arena_Quadrants{5,2},Arena_Quadrants{6,2},'Color',[0.6, 0.6, 0.6],'linewidth',2)
% subplot 223
% plot(Arena_Quadrants{5,3},Arena_Quadrants{6,3},'Color',[0.6, 0.6, 0.6],'linewidth',2)
% subplot 224
% plot(Arena_Quadrants{5,4},Arena_Quadrants{6,4},'Color',[0.6, 0.6, 0.6],'linewidth',2)
% 
%% last update 16/04/2020 - 00:59
% %  listening: Early day miners - East Berlin at night
