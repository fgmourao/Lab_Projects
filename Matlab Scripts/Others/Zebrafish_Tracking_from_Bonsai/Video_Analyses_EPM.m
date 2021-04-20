
%% Analysis of behavioral changes in zebrafish (Danio rerio) larvae. 
%  Elevated plus maze - One larvae

% X and Y coordinates (pixels) related to displacement and quadrants entries (value 1) are extracted from the bonsai software
%   - https://bonsai-rx.org/
%   - https://doi.org/10.3389/fninf.2015.00007

% Files needed:
% - *.csv files from Bonsai software
% - *.png files from Bonsai software

% - Common *.csv files organization. 
%        first  columns -> x coordinates (pixels)
%        second columns -> y coordinates (pixels)
%           3th column  -> upper arm
%           4th column  -> lower arm
%           5th column  -> right arm
%           6th column  -> left  arm
%           7th column  -> center
%                lines  -> values 

% - Analyse one or more *.csv/*.png files

% - Along the trajectory some values may fail so,
%    the empty spaces are filled through a autoregressive modeling. 
%   'Fillgaps' built function


% Outputs:

% Variables : 

% data       -> raw values
%               each cell correspond to a different experiment

% FileName   -> files downloaded
% Header     -> files information and parameters to analyse
% last_frame -> image of the last frame to take the resolution parameters and plot 
%               each cell correspond to a different image (other experiment)

% Analyze

% Arena_full -> in each variable each cell corresponds to one *.csv file download

%            ANALYSED DATA                          VARIABLE
%   - x and y coordinates in pixels              -> vec_y/x
%   - x and y coordinates in cm                  -> vec_y/x_cm
%   - x / y first derivative in cm               -> d_vec_y/x_cm
%   - Displacement in cm                         -> Displacement in cm above threshold
%   - Accumulated distance over time in cm       -> Accumulate_distance
%   - Total distance covered in cm               -> Total_distance
%   - Time vector                                -> Time_vector
%   - Velocity for each movement over time       -> Velocity
%   - Mean Velocity for each experiment          -> Mean_Velocity
%   - Total number of movements                  -> Movements
%   - Total time in movement                     -> Time_Movement
%   - Total time in resting                      -> Time_Resting

%   --------

% quadrants -> in each variable each cell corresponds to one *.csv file download

%            ANALYSED DATA                              VARIABLE
%   - number of entries in each quadrant (raw value)     -> Q

%     Obs. When the animal is at the intersection its presence may appear in more than one quadrant
%     ...then to normalize: 

%   1) the number of entries between the quadrants       -> colDif
%      was sequentially subtracted
%      Q1 - Q2 - Q3 - Q4 - Q5

%   2) only full entries were considered                 -> idx(index)
%      == 1

%   3) "time" spent "over time"                          -> TS

%   4) Total time spent in seconds (pixels/frame rate)   -> Total_TS
%      in each quadrant

%   5) Total number of entries                           -> Total_Entries_Q
%      in each quadrant


% Plots
% - Fillgaps -> Original values (black) and estimated values (red)
% - Tracking map
% - "time" spent "over time"
% - Number of entries in each quadrant

% by Flavio Mourao. Nucleo de Neurociencias - NNC.
% email: mourao.fg@gmail.com
% Universidade Federal de Minas Gerais. 
% Started in:  01/2021
% Last update: 02/2021

%% Load and organize data.

% Load datafiles (*.csv)
[FileName,PathName] = uigetfile({'*.*'},'MultiSelect', 'on'); % Define file type *.*

% Filename can be organize as a single char or a group char in a cell depending on the number os files selected

% data       = each cell corresponds to data from 1 experiment
% last_frame = each cell corresponds to the last frame from specific filename

data = cell(size(FileName));
last_frame = cell(size(FileName));

for ii = 1:length(FileName)
    
    if length(FileName) == 1 % condition for a single file selected 
       Header.FilePattern = dir(fullfile(PathName, char(FileName)));
       [~, ~, fExt] = fileparts(FileName);
        
       switch lower(fExt)

       case '.csv' 
       output = readtable([PathName '/' FileName],'Delimiter',',');
       data{1,ii} = table2array(output);
       Header.Filename_csv(ii) = FileName{ii};
       fprintf(1, 'Reading %s\n', FileName{ii});

       case '.png' 
       last_frame{ii} = imread([PathName '/' FileName]);
       Header.Filename_png(ii) = FileName{ii};
       fprintf(1, 'Reading %s\n', FileName{ii});
       
       end

    else      % condition for multiple files selected
        Header.FilePattern = dir(fullfile(PathName, char(FileName{ii})));
        [~, ~, fExt] = fileparts(FileName{ii});

        switch lower(fExt)

        case '.csv' 
        output = readtable([PathName '/' FileName{ii}],'Delimiter',',');
        data{1,ii} = table2array(output);
        Header.Filename_csv{ii} = FileName{ii};
        fprintf(1, 'Reading %s\n', FileName{ii});

        case '.png' 
        last_frame{ii} = rgb2gray(imread([PathName '/' FileName{ii}])); % open image and convert to grey scale
        Header.Filename_png{ii} = FileName{ii};
        fprintf(1, 'Reading %s\n', FileName{ii});
        
        end
   end
end 

% Remove empty cells
data = data(~cellfun('isempty',data));
last_frame = last_frame(~cellfun('isempty',last_frame));
Header.Filename_csv = Header.Filename_csv(~cellfun('isempty',Header.Filename_csv));
Header.Filename_png = Header.Filename_png(~cellfun('isempty',Header.Filename_png));

%    val_name = FileName(1:regexp(FileName,'.csv')-1);
%    assignin('base',val_name, data);

clear ('fExt', 'ii', 'val_name','PathName','output','temp','FileName');

%% Running data

for ii = 1:length(data)
    
    %  just calling ...
    file = regexprep(Header.Filename_csv{ii} ,'.csv','...');
    fprintf(1, 'Analyzing %s\n', file);
    
%% Fill empty spaces using autoregressive modeling. 
    
    % Fillgaps
    temp_x_y = cell(size(data));    
    temp_x_y{ii} = fillgaps(data{ii}(:,1:2));


    % Plot to check each weel. One experiment per figure
    
    figure
       
    plot(temp_x_y{ii}(:,1),temp_x_y{ii}(:,2),'r-','linew',2)
    hold
    plot(data{ii}(:,1),data{ii}(:,2),'k-','linew',2)    
    title({'Signal Reconstruction. Autoregressive modeling. "Fillgaps" matlab function';[]})

    legend('Estimated values','location','southoutside')

    % Replace original data with estimated values
    data{ii}(:,1:2) = temp_x_y{ii};
    
    % Save Figure
    newStr = regexprep(Header.Filename_png{ii} ,'.png','_');
    
    name = strcat(Header.FilePattern.folder,'\',newStr,'Estimated_values');
    saveas(gcf,name,'png')
    
    close all
    clear('file','jj','temp_x_y','newStr','name')


%% Analysis

    % Set the video parameters
    % - Request from user -
    prompt        = {'Frames per second:','Arena size -> Width (cm):','Arena size -> Height (cm):'};
    dlgtitle      = 'Define values to convertion';
    dims          = [1 40];
    default_input = {'30','6','6'};

    input = inputdlg(prompt,dlgtitle,dims,default_input); %gui

    % Frame Rate
    Header.Num_frames(ii)   = str2double(input{1, 1});

    % Pixels area (check the bonsai values defined)    
    Header.Video_width(ii)  = size(last_frame{ii},2);
    Header.Video_height(ii) = size(last_frame{ii},1);
    
    % Arena Dimensions
    Header.Arena_width(ii)  = str2num(input{2, 1}); % in cm
    Header.Arena_height(ii) = str2num(input{3, 1}); % in cm

    % Define conversion Analyze.factor (cm from pixels)
    Analyze.factor_w(ii) = Header.Arena_width(ii)/Header.Video_width(ii);   % width
    Analyze.factor_h(ii) = Header.Arena_height(ii)/Header.Video_height(ii); % height
    
    % NORMALIZATION :     
    % Size in cm of one fucking pixel
    % Header.Arena_width (in cm) ------- Header.Video_width (in pixels)
    %          x -------------------------------- 1 pixel
    
    % Motion detection threshold - considering the average size of a larvae equal to 3mm
    Analyze.Mov_threshold{ii} = 0.02; %(cm) - 7% seven percent of displacement considering the average size of the larvae
    
    clear ('default_input','dims','dlgtitle','input','prompt')

%% Analyse considering the entire arena

    % Extracting parameters
    
    Analyze.numFrames = data{ii}(:,1);

    Analyze.Arena_Full.vec_y{ii} = data{ii}(:,2); % y axis - in pixels
    Analyze.Arena_Full.vec_x{ii} = data{ii}(:,1); % x axis - in pixels

    Analyze.Arena_Full.vec_y_cm{ii} = data{ii}(:,2).* Analyze.factor_h(1,ii); % y axis - convert y values to cm from pixels
    Analyze.Arena_Full.vec_x_cm{ii} = data{ii}(:,1).* Analyze.factor_w(1,ii); % x axis - convert x values to cm from pixels

    Analyze.Arena_Full.d_vec_y_cm{ii} = [0 ; diff(data{ii}(:,2))].* Analyze.factor_h(1,ii); % convert y Displacement to cm from pixels
    Analyze.Arena_Full.d_vec_x_cm{ii} = [0 ; diff(data{ii}(:,1))].* Analyze.factor_w(1,ii); % convert x Displacement to cm from pixels

    Analyze.Arena_Full.Displacement_raw_values{ii} = sqrt(Analyze.Arena_Full.d_vec_x_cm{ii} .^2 + Analyze.Arena_Full.d_vec_y_cm{ii} .^2); % Displacement in cm
    Analyze.Arena_Full.Displacement{ii}            = Analyze.Arena_Full.Displacement_raw_values{ii};
    
    Analyze.Arena_Full.Displacement{ii}(Analyze.Arena_Full.Displacement{ii} < Analyze.Mov_threshold{ii}) = 0;                             % Displacement in cm above threshold

    Analyze.Arena_Full.Accumulate_distance{ii} = cumsum(Analyze.Arena_Full.Displacement{ii});                                             % Accumulated distance in cm
    Analyze.Arena_Full.Total_distance{ii}      = max(Analyze.Arena_Full.Accumulate_distance{ii});                                         % Total distance in cm
    Analyze.Arena_Full.Time_vector{ii}         = linspace(0,length(data{ii})/Header.Num_frames(1,ii),length(data{ii}));                   % Time vector in sec.

    Analyze.Arena_Full.Velocity{ii}            = Analyze.Arena_Full.Displacement{ii}./[0 diff(Analyze.Arena_Full.Time_vector{ii})]';      % Velocity in cm/s
    Analyze.Arena_Full.Mean_Velocity{ii}       = nanmean(Analyze.Arena_Full.Velocity{ii});                                                % Mean Velocity
    
    Analyze.Arena_Full.Movements{ii}           =  sum(Analyze.Arena_Full.Displacement{ii}>0);                                             % Total number of movements
    Analyze.Arena_Full.Time_Movement{ii}       =  sum(Analyze.Arena_Full.Displacement{ii}>0)*(1/Header.Num_frames(ii));                   % Total time in movement    
    Analyze.Arena_Full.Time_Resting{ii}        =  sum(Analyze.Arena_Full.Displacement{ii}==0)*(1/Header.Num_frames(ii));                  % Total time in resting
   
    % Kernel smoothing function to compute the probability density estimate (PDE) for each point. 
    Analyze.Arena_probability_density{ii}      = ksdensity([Analyze.Arena_Full.vec_x{1},Analyze.Arena_Full.vec_y{1}],...
                                                [Analyze.Arena_Full.vec_x{1},Analyze.Arena_Full.vec_y{1}]);

    % Extracting parameters and analize entries in each quadrant
 
    Analyze.quadrants.Q{ii}(:,1) = data{ii}(:,3); % upper arm   
    Analyze.quadrants.Q{ii}(:,2) = data{ii}(:,4); % lower arm   
    Analyze.quadrants.Q{ii}(:,3) = data{ii}(:,5); % right arm   
    Analyze.quadrants.Q{ii}(:,4) = data{ii}(:,6); % left  arm
    Analyze.quadrants.Q{ii}(:,5) = data{ii}(:,7); % center  
 

    % Normalized Entries

    Analyze.quadrants.colDif{ii} = abs(Analyze.quadrants.Q{ii}(:,1) - Analyze.quadrants.Q{ii}(:,2) - Analyze.quadrants.Q{ii}(:,3) - Analyze.quadrants.Q{ii}(:,4) - Analyze.quadrants.Q{ii}(:,5));
    Analyze.quadrants.idx{ii} = Analyze.quadrants.colDif{ii} == 1; % only full entries were considered
    
    Analyze.quadrants.TS{ii} = Analyze.quadrants.Q{ii}(Analyze.quadrants.idx{ii},:); % time spent over time
    Analyze.quadrants.Entries_Q{ii} = diff(Analyze.quadrants.TS{ii});                % entries over time

    Analyze.quadrants.Total_TS{ii} =  sum(Analyze.quadrants.TS{ii})./Header.Num_frames(1,ii);   % Total time spent (pixels/frame rate)
    Analyze.quadrants.Total_Entries_Q{ii} = sum(Analyze.quadrants.Entries_Q{ii}>0);             % Total entries



%% Plot one maze

    figure
    set(gcf,'color','w');
        
    subplot(3,3,2)
    plot(Analyze.quadrants.TS{ii}(:,1),'k','linew',2)
    ylim([0 1.5])
    title('upper arm')

    subplot(3,3,4)
    plot(Analyze.quadrants.TS{ii}(:,4),'k','linew',2)
    ylim([0 1.5])
    title('left  arm')

    subplot(3,3,5)
    plot(Analyze.quadrants.TS{ii}(:,5),'k','linew',2)
    ylim([0 1.5])
    title('center')

    subplot(3,3,6)
    plot(Analyze.quadrants.TS{ii}(:,3),'k','linew',2)
    ylim([0 1.5])
    title('right arm')

    subplot(3,3,8)
    plot(Analyze.quadrants.TS{ii}(:,2),'k','linew',2)
    ylim([0 1.5])
    title('lower arm')
    
    suptitle ('Time spent in each arm')
    
    % Save Figure
    newStr = regexprep(Header.Filename_png{ii} ,'.png','_');
    
    name = strcat(Header.FilePattern.folder,'\',newStr,'Time_spent_in_each_arm');
    saveas(gcf,name,'png')
    
    close all
    
    figure
    set(gcf,'color','w');
        
    subplot(3,3,2)
    plot(Analyze.quadrants.Entries_Q{ii}(:,1),'k','linew',2)
    ylim([0 1.5])
    title('upper arm')

    subplot(3,3,4)
    plot(Analyze.quadrants.Entries_Q{ii}(:,4),'k','linew',2)
    ylim([0 1.5])
    title('left  arm')

    subplot(3,3,5)
    plot(Analyze.quadrants.Entries_Q{ii}(:,5),'k','linew',2)
    ylim([0 1.5])
    title('center')

    subplot(3,3,6)
    plot(Analyze.quadrants.Entries_Q{ii}(:,3),'k','linew',2)
    ylim([0 1.5])
    title('right arm')

    subplot(3,3,8)
    plot(Analyze.quadrants.Entries_Q{ii}(:,2),'k','linew',2)
    ylim([0 1.5])
    title('lower arm')
    
    suptitle ('Number of entries')
    
    % Save Figure
    newStr = regexprep(Header.Filename_png{ii} ,'.png','_');
    
    name = strcat(Header.FilePattern.folder,'\',newStr,'Number_of_entries');
    saveas(gcf,name,'png')
    
    close all
    
    clear('newStr','name')
    
%% Plot track full arena

    figure;
    set(gcf,'color','w');

    imshow(last_frame{ii})
    hold on
    plot(Analyze.Arena_Full.vec_x{ii},Analyze.Arena_Full.vec_y{ii},'r-','linewidth',2)
    axis off
%     xlim ([0 Header.Arena_width]);  % in cm
%     ylim ([0 Header.Arena_height]); % in cm
    
    % Scale bar size
    sb = round(.5 / Analyze.factor_w(1,ii)); % 5 mm
    
    % Scale bar location
    x1 = size(last_frame{1,1},2) - 150;
    x2 = x1 + sb;
    y = size(last_frame{1,1},2) - 150;
    
    plot([x1 x2],[y y],'k-','linew',3)
    
    suptitle(['Displacement Tracking - Experiment ' num2str(ii)])
    
    % Save Figure
    newStr = regexprep(Header.Filename_png{ii} ,'.png','_');
    
    name = strcat(Header.FilePattern.folder,'\',newStr,'Single_tracking');
    saveas(gcf,name,'png')
    
    close all
    
    clear('sb','x1','x2','y','newStr','name')
    
%% Plot tracking as a function of probability density estimate
%  Plot where each point is colored by the spatial density of nearby points. 

    pic = repmat(last_frame{1},[1 1 3]);       % simulating RGB to plot colors 
    p = Analyze.Arena_probability_density{ii}; % doubling the variable just to shorten the name

    figure
    set(gcf,'color','w');

    imshow(pic)
    hold on
    % plot3(a(:,1),a(:,2),c,'linew',10)

    % plot3 only supports a single color.
    % You can use other graphics functions. The "traditional" one for this purpose 
    % is actually rather surprising. It's the patch function, which is designed for drawing filled polygons.
    % What's going on here is the following:
    % - The 4th arg says to use Z for color data
    % - The EdgeColor=interp says to interpolate the color data as the edge color
    % - The FaceColor=none says to not fill the polygon. Just draw the edges
    % - The nans say not to connect the last point to the first point.

    % Solution from here: https://www.mathworks.com/matlabcentral/answers/267468-change-colour-of-points-in-plot3-with-increasing-z-value

    patch([Analyze.Arena_Full.vec_x{ii} nan(size(Analyze.Arena_Full.vec_x{ii}))],...
        [Analyze.Arena_Full.vec_y{ii} nan(size(Analyze.Arena_Full.vec_y{ii}))],[p nan(size(p))],[p nan(size(p))],'EdgeColor','interp','FaceColor','none','LineWidth',2)

    cb = colorbar();
    cb.Label.String = 'Probability density estimate';
    colormap(jet(4096))

    % Scale bar size in pixels
    scale_size = 0.5; %(cm)
    sb = round(scale_size / Analyze.factor_w(1,ii)); % 5 mm
    
    %           1(one ficking pixel) - Analyze.factor_w (size in cm of 1 pixel)
    %                    x  ----------      0.5 (cm)
    
    % Scale bar location
    x1 = size(last_frame{1,1},2) - 150;
    x2 = x1 + sb;
    y = size(last_frame{1,1},2) - 150;
    
    plot([x1 x2],[y y],'k-','linew',3)
    
    suptitle(['Displacement Tracking - Experiment ' num2str(ii)])
    
    % Save Figure
    newStr = regexprep(Header.Filename_png{ii} ,'.png','_');
    
    name = strcat(Header.FilePattern.folder,'\',newStr,'PDE_tracking');
    saveas(gcf,name,'png')
    
    close all

    clear('p','pic','cb','sb','x1','x2','y','name','newStr','scale_size')

%% Plot heatmap as a function of probability density estimate

    pic = repmat(last_frame{1},[1 1 3]);       % simulating RGB to plot colors 
    p = Analyze.Arena_probability_density{ii}; % doubling the variable just to shorten the name

    % Perform interpolation to smooth the data. Built in Function:'scatteredInterpolant'
    % scatteredInterpolant uses a Delaunay triangulation of the scattered sample points to perform interpolation

    % Create the interpolant.
    F = scatteredInterpolant(Analyze.Arena_Full.vec_x{1},Analyze.Arena_Full.vec_y{1},p);

    % Optional to increase resolution 
    factor = 1; % change multiplication factor

    newNumberOfRows = factor*(size(last_frame{1,1},1));
    newNumberOfCols = factor*(size(last_frame{1,1},2));

    % Grid of query points based on Video Resolution
    [xq, yq] = meshgrid(linspace(1, size(last_frame{1,1},2), newNumberOfCols), linspace(1, size(last_frame{1,1},1), newNumberOfRows));

    % Evaluate the interpolant at query locations (xq,yq).
    % Use the 'nearest', 'linear', or 'natural' methods

    F.Method = 'natural';
    pq = F(xq,yq);

    % Remove negative values
    pq(pq<0) = 0;

    figure;
    set(gcf,'color','w');
    
    imshow(pic) % plot the original picture below just to ensure the correct orientation
    hold on

    ax = newplot;
    surf(xq,yq,pq,...
             'EdgeColor','none',...
             'FaceColor','interp');    
    colorbar('Location','eastoutside')%,'YTick'),[]);
    view(ax,2);
    grid(ax,'off');
    colormap jet
    caxis([0.4*10^-5 4*10^-5]) % Be careful here. Set a pattern for all figures
    cb = colorbar();
    cb.Label.String = 'Probability density estimate';

    % Drawing squares to simulate the maze. 
    annotation('rectangle',...
        [0.571049136786188 0.144215530903329 0.233731739707835 0.256735340729002],...
        'Color',[1 1 1],...
        'FaceColor',[0.941176470588235 0.941176470588235 0.941176470588235]);

    annotation('rectangle',...
        [0.114209827357237 0.14421553090333 0.233731739707835 0.256735340729002],...
        'Color',[1 1 1],...
        'FaceColor',[0.941176470588235 0.941176470588235 0.941176470588235]);

    annotation('rectangle',...
        [0.114209827357238 0.654516640253566 0.233731739707835 0.256735340729002],...
        'Color',[1 1 1],...
        'FaceColor',[0.941176470588235 0.941176470588235 0.941176470588235]);

    annotation('rectangle',...
        [0.571049136786188 0.654516640253567 0.233731739707835 0.256735340729002],...
        'Color',[1 1 1],...
        'FaceColor',[0.941176470588235 0.941176470588235 0.941176470588235]);

    % Scale bar location in matlab units
    pos = get(gca, 'Position');
    
    % Scale bar size as function of Matlab units
    unit = pos(3)/size(last_frame{1,1},2); % define size of one unit normalizing by the figure width 
                                           % default width of the plot within the figure in normalized coordinates 
                                           % (i.e. the figure / the window containing the plot has a width of 1 length units).  
    
                                           
%      NORMALIZATION :                                          
%      0.6910(default width of the plot in units ---  1158(figure size in pixels) ------ 6 (figure size in cm)
%                   x ------------------------------------------ 1 (one fucking pixel) - Analyze.factor_w (size in cm of 1 pixels)
%                
%                   x = 5.9672e-04 ( size of 1 unit)
                  

    % define scale bar size in matlab units
    scale_size = 0.5; %(cm)
    sb = (scale_size * unit) / Analyze.factor_w(1,ii); % 5 mm
 
%         unit - Analyze.factor_w
%           x  -      0.5
    
    x1 = pos(3) - sb;
    x2 = pos(3);
    y = 1 - pos(4);
    
    annotation('line',[x1 x2],...
    [y y],'Color',[0 0 0],'LineWidth',3);

    %plot([x1 x2],[y y],'k-','linew',3)
    
    % Save Figure
    newStr = regexprep(Header.Filename_png{ii} ,'.png','_');
    
    name = strcat(Header.FilePattern.folder,'\',newStr,'PDE_Interp_tracking');
    saveas(gcf,name,'png')
    
    close all

    clear('p','pic','cb','sb','x1','x2','y','F','ax','scale_size','newStr','name','newNumberOfCols','newNumberOfRows','pos','pq','unit','xq','yq','factor')


    
end


%% Save
    
% newStr = regexprep(Header.Filename_csv{ii} ,'.csv','');  % save 1 experiment 
newStr = EPM_maze;  
name = strcat(Header.FilePattern.folder,'\',newStr);
save(name,'Analyze','-v7.3')

clear('newStr','name')

%% Save Header
    
% newStr = regexprep(Header.Filename_csv{ii} ,'.csv','_Header'); % save 1 experiment  
newStr = EPM_maze_header;
name = strcat(Header.FilePattern.folder,'\',newStr);
save(name,'Header','-v7.3')
    
clear('newStr','name','ii')  

%% last update 10/02/2021 - 12:12
% listening: this will destroy you - they move on tracks...
