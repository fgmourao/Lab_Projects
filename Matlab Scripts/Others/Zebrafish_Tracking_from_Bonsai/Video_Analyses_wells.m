
%% Analysis of behavioral changes in zebrafish larvae. 
%  Experiments on x well plate - One larva per well

% X and Y coordinates related to displacement are extracted from the bonsai software
%   - https://bonsai-rx.org/
%   - https://doi.org/10.3389/fninf.2015.00007

% Files needed:
% - *.csv files from Bonsai software
% - *.png files from Bonsai software

% - Common *.csv files organization. 
%        odd  columns -> x coordinates (pixels)
%        even columns -> y coordinates (pixels)
%             lines   -> values 

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

% Arena_full -> to each variable each cell corresponds to one *.csv file download

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

% Plots
% - Fillgaps -> Original values (black) and estimated values (red)
% - Tracking map to each well
% - Displacement over time (cm)
% - Accumulated distance over time
% - Velocity to each displacement over time

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


%%

for ii = 1:length(data)

    %  just calling ...
    file = regexprep(Header.Filename_csv{ii} ,'.csv','...');
    fprintf(1, 'Analyzing %s\n', file);    

%% Analysis Settings
   
%     % Set the video parameters
%     % Request from user
    prompt        = {'Frames per second:','Plate Width (cm):',...
        'Plate Height (cm):','Well Diameter (cm)','Empty Wells (Enter space separated values)' };

    dlgtitle      = ['Experiment number ' num2str(ii)];
    dims          = [1 45];
    default_input = {'30','12.8','8.6','2.1','10 11 12'};

    input = inputdlg(prompt,dlgtitle,dims,default_input); %gui


    % Frame Rate
    Header.Num_frames(ii)   = str2double(input{1, 1});

    % Pixels area (check the bonsai values defined)
    Header.Video_width(ii)  = size(last_frame{ii},2);
    Header.Video_height(ii) = size(last_frame{ii},1);

    % Plate area (check Tissue culture plate, 12 wells datasheed)
    Header.Plate_width(ii)  = str2num(input{3, 1});
    Header.Plate_height(ii) = str2num(input{2, 1});

    % Arena Dimensions
    Header.Well_Diameter(ii)  = str2num(input{4, 1}); % in cm

    % Empty wells
    Header.Empty_Wells{ii}    = str2num(input{5, 1}); % in cm


%     % % Set the video parameters
%     % % - Manually -
%    
%     % Frame Rate
%     Header.Num_frames(ii)   = 30;
% 
%     % Pixels area (check the bonsai values defined)
%     Header.Video_width(ii)  = size(last_frame{ii},2);
%     Header.Video_height(ii) = size(last_frame{ii},1);
% 
%     % Plate area (check Tissue culture plate, 12 wells datasheed)
%     Header.Plate_width(ii)  = 12.8;
%     Header.Plate_height(ii) = 8.6;
% 
%     % Arena Dimensions
%     Header.Well_Diameter(ii)  = 2.1; % in cm
% 
%     % Empty wells
%     Header.Empty_Wells{ii}    = [9 10 11 12]; % in cm


    % Clear empty wells

        for jj = 1:length(Header.Empty_Wells{ii})
            data{1,ii}(:,Header.Empty_Wells{ii}(1,jj)*2 - 1 : Header.Empty_Wells{ii}(1,jj)*2)  = nan(size(length(data{1,ii}),2));
        end


    % Define conversion factor -> Analyze.factor (cm from pixels)
    Analyze.factor_w(ii) = Header.Plate_width(ii)/Header.Video_width(ii);   % width
    Analyze.factor_h(ii) = Header.Plate_height(ii)/Header.Video_height(ii);  % height
    
    % NORMALIZATION :
    % Size in cm of one fucking pixel
    % Header.Arena_width (in cm) ------- Header.Video_width (in pixels)
    %          x -------------------------------- 1 pixel
    
    % Motion detection threshold 
    Analyze.Mov_threshold{ii} = 0.03; %(cm)
    
    % Considering that the average size of a larvae equal to ~3 mm (REF), 
    % one single movement was defined when the larvae displacement was greater 
    % than ~ 3 pixels (0.2 mm / 7% total larva size) (REF = 10.1016/j.beproc.2010.12.003). 
   
    % other possible ref:                          
    % https://www.jove.com/t/54431/using-touch-evoked-response-locomotion-assays-to-assess-muscle
    % inactivity if the larva moves less than 6 mm per second (camera 100 fps)
            
    clear ('default_input','dims','dlgtitle','input','prompt','jj')
    
%% Fill empty spaces using autoregressive modeling. 
    
    % Fillgaps
    temp_x_y = cell(size(data));    
    temp_x_y{ii} = fillgaps(data{ii});


    % Plot to check each weel. One experiment per figure
    
    titles{ii} = 1:size(data{ii},2)/2;

    figure
    suptitle({'Data Reconstruction';[]})

    for jj = 2:2:size(data{ii},2)
        subplot(3,4,jj/2) 
        plot(temp_x_y{ii}(:,jj-1),temp_x_y{ii}(:,jj),'r-','linew',2)
        hold
        plot(data{ii}(:,jj-1),data{ii}(:,jj),'k-','linew',2)
        title(titles{ii}(1,jj/2),'FontSize',11);        
    end
    
    legend('Estimated values','location','southoutside')

    % Replace original data with estimated values
    data{ii} = temp_x_y{ii};

    % Save Figure
    newStr = regexprep(Header.Filename_png{ii} ,'.png','_');
    
    name = strcat(Header.FilePattern.folder,'\',newStr,'Estimated_values');
    saveas(gcf,name,'png')
    
    close all
    clear('file','jj','temp_x_y','newStr','name','titles')

%% Analyse considering the entire arena
 
    % Extracting parameters
 
    for jj = 2:2:size(data{ii},2)
 
        Analyze.Arena_Full.vec_y{ii}(:,jj/2) = data{ii}(:,jj);   % y axis - in pixels
        Analyze.Arena_Full.vec_x{ii}(:,jj/2) = data{ii}(:,jj-1); % x axis - in pixels

        Analyze.Arena_Full.vec_y_cm{ii}(:,jj/2) = data{ii}(:,jj).* Analyze.factor_h(1,ii);   % y axis - convert y values to cm from pixels
        Analyze.Arena_Full.vec_x_cm{ii}(:,jj/2) = data{ii}(:,jj-1).* Analyze.factor_w(1,ii); % x axis - convert x values to cm from pixels

        Analyze.Arena_Full.d_vec_y_cm{ii}(:,jj/2) = [0 ; diff(data{ii}(:,jj))].* Analyze.factor_h(1,ii); % convert y Displacement to cm from pixels
        Analyze.Arena_Full.d_vec_x_cm{ii}(:,jj/2) = [0 ; diff(data{ii}(:,jj-1))].* Analyze.factor_w(1,ii); % convert x Displacement to cm from pixels

    end
    
    Analyze.Arena_Full.Displacement_raw_values{ii} = sqrt(Analyze.Arena_Full.d_vec_x_cm{ii} .^2 + Analyze.Arena_Full.d_vec_y_cm{ii} .^2); % Displacement in cm
    Analyze.Arena_Full.Displacement{ii}            = Analyze.Arena_Full.Displacement_raw_values{ii};
    
    Analyze.Arena_Full.Displacement{ii}(Analyze.Arena_Full.Displacement{ii} < Analyze.Mov_threshold{ii}) = 0;                             % Displacement in cm above threshold

    Analyze.Arena_Full.Accumulate_distance{ii} = cumsum(Analyze.Arena_Full.Displacement{ii});                                             % Accumulated distance in cm
    Analyze.Arena_Full.Total_distance{ii}      = max(Analyze.Arena_Full.Accumulate_distance{ii});                                         % Total distance in cm
    Analyze.Arena_Full.Time_vector{ii}         = linspace(0,length(data{ii})/Header.Num_frames(ii),length(data{ii}));                     % Time vector in sec.

    Analyze.Arena_Full.Velocity{ii}            = Analyze.Arena_Full.Displacement{ii}./[0 diff(Analyze.Arena_Full.Time_vector{ii})]';      % Velocity in cm/s
    Analyze.Arena_Full.Mean_Velocity{ii}       = nanmean(Analyze.Arena_Full.Velocity{ii});                                                % Mean Velocity
    
    Analyze.Arena_Full.Movements{ii}           =  sum(Analyze.Arena_Full.Displacement{ii}>0);                                             % Total number of movements
    Analyze.Arena_Full.Time_Movement{ii}       =  sum(Analyze.Arena_Full.Displacement{ii}>0)*(1/Header.Num_frames(ii));                    % Total time movement    
    Analyze.Arena_Full.Time_Resting{ii}        =  sum(Analyze.Arena_Full.Displacement{ii}==0)*(1/Header.Num_frames(ii));                  % Total time in resting
    
    clear ('jj')
 
%% Plot track in each real well

    figure('units','normalized','outerposition',[0 0 1 1]);
    set(gcf,'color','w');
 
    imshow(last_frame{ii})
    title(['Experiment ' num2str(ii)],'FontSize',12);
    hold on
    
    for jj = 1:size(data{ii},2)/2
        %subplot(3,4,jj)
    
        plot(Analyze.Arena_Full.vec_x{ii}(:,jj),Analyze.Arena_Full.vec_y{ii}(:,jj),'r-','linewidth',1)

        axis off
    end  
 
    suptitle('Displacement Tracking in each well')
    
    % Save Figure
    newStr = regexprep(Header.Filename_png{ii} ,'.png','_');
    
    name = strcat(Header.FilePattern.folder,'\',newStr,'track');
    saveas(gcf,name,'png')
    
    close all
    
    clear('jj','newStr','name')
    
%% Plot track in each real well as a function of velocity
    
    pic = repmat(last_frame{1},[1 1 3]);       % simulating RGB to plot colors 
    v = Analyze.Arena_Full.Velocity{ii};       % doubling the variable just to shorten the name

    figure('units','normalized','outerposition',[0 0 1 1]);
    set(gcf,'color','w');
 
    imshow(pic)
    title(['Experiment ' num2str(ii)],'FontSize',12);
    hold on
    
    for jj = 1:size(data{ii},2)/2
        %subplot(3,4,jj)
    
        plot(Analyze.Arena_Full.vec_x{ii}(:,jj),Analyze.Arena_Full.vec_y{ii}(:,jj),'r-','linewidth',1)
        
        patch([Analyze.Arena_Full.vec_x{ii}(:,jj) nan(size(Analyze.Arena_Full.vec_x{ii}(:,jj)))],...
            [Analyze.Arena_Full.vec_y{ii}(:,jj) nan(size(Analyze.Arena_Full.vec_y{ii}(:,jj)))],[v(:,jj) nan(size(v(:,jj)))],[v(:,jj) nan(size(v(:,jj)))],'EdgeColor','interp','FaceColor','none','LineWidth',2)

    cb = colorbar();
    cb.Label.String = 'Velocity (cm/s)';
    colormap(jet(4096))
    caxis([0 6])
        axis off
    end  
 
    suptitle('Displacement tracking as a function of velocity')
    
    % Save Figure
    newStr = regexprep(Header.Filename_png{ii} ,'.png','_');
    
    name = strcat(Header.FilePattern.folder,'\',newStr,'track_velocity');
    saveas(gcf,name,'png')
    
    close all
    
    clear('v','pic','jj','cb','newStr','name')
    
%% Plot Displacement
 
    figure('units','normalized','outerposition',[0 0 1 1]);
    set(gcf,'color','w');
    
    titles{ii} = 1:size(data{ii},2)/2;

    for jj = 1:size(data{ii},2)/2
        subplot(3,4,jj)
        hold on
        plot(Analyze.Arena_Full.Time_vector{ii},Analyze.Arena_Full.Displacement{ii}(:,jj)','Color',[0.6, 0.6, 0.6],'linewidth',1)
        title(['Well ' num2str(titles{ii}(1,jj))],'FontSize',11);
 
        xlabel('Time (s)')
        ylabel('Displacement (cm)')
        
        xlim ([0 Analyze.Arena_Full.Time_vector{ii}(end)]); 
        %ylim ([0 max(max(Analyze.Arena_Full.Displacement{ii}))]);
        
        box off
        
       % plot([0 Analyze.Arena_Full.Time_vector{ii}(end)],[0.0200 0.0200]) % threshold

    end

    % Save Figure
    newStr = regexprep(Header.Filename_png{ii} ,'.png','_');
    
    name = strcat(Header.FilePattern.folder,'\',newStr,'displacement');
    saveas(gcf,name,'png')
    
    close all
    
    clear('titles','jj','newStr','name')    
 
%% Plot Displacement as a function of Probablity Density Distribuition
 
    figure('units','normalized','outerposition',[0 0 1 1]);
    set(gcf,'color','w');
    
    titles{ii} = 1:size(data{ii},2)/2;

    for jj = 1:size(data{ii},2)/2
        subplot(3,4,jj)
        hold on
        
        d = Analyze.Arena_Full.Displacement{ii}(:,jj);% doubling the variable just to shorten the name

        h = histogram(d(d>0),50,'Normalization','pdf');
        h.FaceColor = [1 1 1];
        h.EdgeColor = 'k';
        
        y = min(d(d>0)):0.001:max(d(d>0));
        mu = mean(d(d>0));
        sigma = std(d(d>0),0,1);
        f = exp(-(y-mu).^2./(2*sigma^2))./(sigma*sqrt(2*pi));
        plot(y,f,'LineWidth',1.5)

        title(['Well ' num2str(titles{ii}(1,jj))],'FontSize',11);
 
        xlabel('Displacement (cm)')
        ylabel('Probability density function')
        
        
        box off
        
       % plot([0 Analyze.Arena_Full.Time_vector{ii}(end)],[0.0200 0.0200]) % threshold

    end

 
    % Save Figure
    newStr = regexprep(Header.Filename_png{ii} ,'.png','_');
    
    name = strcat(Header.FilePattern.folder,'\',newStr,'displacement_PDE');
    saveas(gcf,name,'png')
    
    close all
    
    clear('titles','d','jj','newStr','name','h','mu','f','sigma','y') 


%% Plot Accumulated Distance
 
    figure('units','normalized','outerposition',[0 0 1 1]);
    set(gcf,'color','w');
 
    titles{ii} = 1:size(data{ii},2)/2;

    for jj = 1:size(data{ii},2)/2
        subplot(3,4,jj)

        plot(Analyze.Arena_Full.Time_vector{ii},Analyze.Arena_Full.Accumulate_distance{ii}(:,jj)','Color',[0.6, 0.6, 0.6],'linewidth',1)
        title(['Well ' num2str(titles{ii}(1,jj))],'FontSize',11);

        xlim ([0 Analyze.Arena_Full.Time_vector{ii}(end)]); 
        %ylim ([0 max(max(Analyze.Arena_Full.Accumulate_distance{ii}))]);

        xlabel('Time (s)')
        ylabel('Accumulated distance (cm)')

        box off
        
    end

 
    clear ('jj','titles')

    % Save Figure
    newStr = regexprep(Header.Filename_png{ii} ,'.png','_');
    
    name = strcat(Header.FilePattern.folder,'\',newStr,'Accumulated_Distance');
    saveas(gcf,name,'png')
    
    close all
    
    clear('titles','jj','newStr','name')
    
%% Plot Velocity
 
    figure('units','normalized','outerposition',[0 0 1 1]);
    set(gcf,'color','w');
    
    titles{ii} = 1:size(data{ii},2)/2;
    
    for jj = 1:size(data{ii},2)/2
        subplot(3,4,jj)
 
        plot(Analyze.Arena_Full.Time_vector{ii},Analyze.Arena_Full.Velocity{ii}(:,jj)','Color',[0.6, 0.6, 0.6],'linewidth',1)
        title(['Well ' num2str(titles{ii}(1,jj))],'FontSize',11);
 
        xlabel('Time (s)')
        ylabel('Velocity (cm/s)')
 
        xlim ([0 Analyze.Arena_Full.Time_vector{ii}(end)]); 
        %ylim ([0 max(max(Analyze.Arena_Full.Velocity{ii}))]);
 
        box off

    end

    % Save Figure
    newStr = regexprep(Header.Filename_png{ii} ,'.png','_');
    
    name = strcat(Header.FilePattern.folder,'\',newStr,'Velocity');
    saveas(gcf,name,'png')
    
    close all
    
    clear('titles','jj','newStr','name')

end  

%% Save
    
% newStr = regexprep(Header.Filename_csv{ii} ,'.csv','');  % save 1 experiment 
newStr = 'tracking_wells';
name = strcat(Header.FilePattern.folder,'\',newStr);
save(name,'Analyze','-v7.3')

clear('newStr','name')

%% Save Header
    
% newStr = regexprep(Header.Filename_csv{ii} ,'.csv','_Header'); % save 1 experiment  
newStr = 'tracking_wells_header';
name = strcat(Header.FilePattern.folder,'\',newStr);
save(name,'Header','-v7.3')
    
clear('newStr','name','ii') 

%% last update 10/02/2021 - 12:10
% listening: this will destroy you - they move on tracks...
