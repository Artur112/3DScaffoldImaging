        %4th Year Biomedical Engineering Final Year Project, written by Artur
%Jurgenson 2018-2019.

%Please run the first section first if you want to run any of the
%subsequent sections individually. All of the analyzed images, results and
%plots are saved in the folder that is specified in the first section. The
%folders that contain the scans from each channel for a week should all be
%placed in the same folder.

%Please note, this program requires the Matlab Image Processing Toolbox to
%be installed.

%% Get address of folder where data is stored and scan parameters
disp('Please select folder where the weekly scans are stored (ie folders that contain the Ch_x scans for each week)');
data_address = uigetdir('temp', 'Select folder where data is');
week_folders = dir(data_address);
dir_keep = [week_folders.isdir] & ~strcmp({week_folders.name},'.') & ~strcmp({week_folders.name},'..');
week_folders = week_folders(dir_keep); 
disp('Please select folder where to store the results');
store_address = uigetdir('temp','Select where to store results');

numweeks = length(week_folders);
nr_total_slices = 0;
for m = 2:numweeks
    nr_total_slices = nr_total_slices + numel(dir([strcat(week_folders(m).folder,'\',week_folders(m).name,'\','Ch2_Stitched_Sections') '\*.tif']));
end
    
disp('Please input slice thicknesses');
prompt = cell(1,numweeks);
for week = 1:numweeks
    prompt{week} = strcat("Week ", num2str(week), ' slice thickness / um:');
end
thicknesses = inputdlg(prompt);

%% Segment the images, convert scaffold scans to binary images and calculate metrics
bar = waitbar(0, {'Finding cells and calculating metrics', 'Estimated time remaining: ...'});
net = load('network1.mat'); %Load neural Network
net = net.net;
picsize = 300;  
%mkdir(strcat(store_address,'\Metrics')); %Create folder where to store calculated metrics
weeknr = 1;
slices_completed = 0;
time_per_slice = [];

for week = week_folders'
    mkdir(strcat(store_address,'\Cells\Week',num2str(weeknr))); %Create folder to store cell images in
    mkdir(strcat(store_address,'\Scaffolds\Week',num2str(weeknr))); %Create folder to store scaffold images in
    scafffiles = dir([strcat(week.folder,'\',week.name,'\','Ch1_Stitched_Sections') '\*.tif']);
    cellfiles = dir([strcat(week.folder,'\',week.name,'\','Ch2_Stitched_Sections') '\*.tif']);
    numslices = numel(cellfiles);
    metric_data = cell(1,numslices); %Pre allocating for speed
    
    for slice = 1:numslices
        tic;
    %%% Segment cell image by sliding a 300x300 section along the image that is segmented with neural network
        cell_img = imread(strcat(cellfiles(slice).folder,'\',cellfiles(slice).name));
        cell_img_segm = zeros(size(cell_img,1),size(cell_img,2));
        for ro = 1:floor(size(cell_img,1)/picsize)
            for co = 1:floor(size(cell_img,2)/picsize)
                section = cell_img(1+picsize*(ro-1) : picsize*ro , 1+picsize*(co-1) : picsize*co); %300x300 section to analyze
                if(any(section,'all')) %Skip if contains no cells          
                    %Find if touching any cells along edge,[left edge, right edge, top edge, bottom edge]
                    touching = [any(section(:,1)),any(section(:,picsize)),any(section(1,:)),any(section(picsize,:))];
                    C = semanticseg(section,net);
                    C = double(C)-1;
                    if(any(touching)) %If touching fill in the gaps in the edges
                        for m = 1:4
                            if(touching(m))
                                switch m
                                    case 1
                                        C(section(:,1)>0,1) = 1;
                                    case 2
                                        C(section(:,picsize) > 0,picsize) = 1;
                                    case 3
                                        C(1,section(1,:) > 0) = 1;
                                    case 4
                                        C(picsize,section(picsize,:)>0) = 1;
                                end
                            end
                        end                      
                    end
                    cell_img_segm(1+picsize*(ro-1) : picsize*ro , 1+picsize*(co-1) : picsize*co) = C;  
                end
            end
        end
        
    %%% Get scaffold image and convert it to binary  
        scaff_img = imread(scafffiles(slice).folder + "\" + scafffiles(slice).name);
        scaff_img = double(scaff_img > 0);
        
    %%% Calculate Metrics  
        cc = bwconncomp(cell_img_segm);
        slicedata.Centroids = cell2mat(struct2cell(regionprops(cc,'Centroid'))'); %Centroids of all cells
        slicedata.Areas = cell2mat(struct2cell(regionprops(cc,'Area'))'); %Areas of all cells
        %   slicedata.BBoxes = cell2mat(struct2cell(regionprops(cc,'BoundingBox'))'); %Bounding boxes of all cells
        MinAx = cell2mat(struct2cell(regionprops(cc,'MinorAxisLength'))');
        MajAx = cell2mat(struct2cell(regionprops(cc,'MajorAxisLength'))');
        slicedata.WHratios = MajAx ./ MinAx; %Width to height ratios of cells
        
        %Find areas that are smaller than 50 (these are cell fragments that count as noise or cells that partially show up from deeper layers.
        toremove = find(slicedata.Areas < 50);
        
        %Remove these areas from the cell image
        pxls = regionprops(cc,'PixelList');
        for ii = 1:length(toremove)
            cell_img_segm(pxls(toremove(ii)).PixelList(:,2), pxls(toremove(ii)).PixelList(:,1)) = 0;
        end
        
        % Remove these metrics from the metric_data
        slicedata.Centroids(toremove,:) = [];
        slicedata.Areas(toremove,:) = [];
        slicedata.WHratios(toremove,:) = [];
        
        % Calculate percentage of scaffold occupied by cells    
        totalscaff = sum(scaff_img,'all');
        slicedata.ScaffPerc = sum(cell_img_segm + scaff_img == 2,'all') / totalscaff * 100;

    %%% Store data and images for current slice
        metric_data{slice} = slicedata;
        if(slice < 10)
            imwrite(cell_img_segm,strcat(store_address,'\Cells\Week',num2str(weeknr),'\Slice_00',num2str(slice),'.png'));
           % imwrite(scaff_img,strcat(store_address,'\Scaffolds\Week',num2str(weeknr),'\Slice_00',num2str(slice),'.png'));
        elseif(slice<100)
            imwrite(cell_img_segm,strcat(store_address,'\Cells\Week',num2str(weeknr),'\Slice_0',num2str(slice),'.png'));
           % imwrite(scaff_img,strcat(store_address,'\Scaffolds\Week',num2str(weeknr),'\Slice_0',num2str(slice),'.png'));
        else
            imwrite(cell_img_segm,strcat(store_address,'\Cells\Week',num2str(weeknr),'\Slice_',num2str(slice),'.png'));
           % imwrite(scaff_img,strcat(store_address,'\Scaffolds\Week',num2str(weeknr),'\Slice_',num2str(slice),'.png'));
        end
    
    %%% Update progress bar
        slices_completed = slices_completed + 1;
        time_per_slice(slices_completed) = toc;
        time_left = mean(time_per_slice)*(nr_total_slices - slices_completed)/60/60;
        hours_left = floor(time_left);
        min_left = round((time_left - hours_left) * 60);
        waitbar(slices_completed/nr_total_slices,bar, {'Finding cells and calculating metrics', strcat('Estimated time remaining: ',' ', num2str(hours_left),' h',' ',num2str(min_left),' min')});
    end
    save(strcat(store_address,'\Metrics\Week',num2str(weeknr),'_metrics.mat'), 'metric_data');
    weeknr = weeknr + 1;
end
delete(bar);
        
%% Visualize Results
%mkdir(strcat(store_address,'\Plots'));
thicknesses{1} = '50';
thicknesses{2} = '70';
thicknesses{3} = '70';
numweeks = 3;
store_address = 'C:\users\Artur\Desktop\Project_data_new';
data = cell(1,numweeks);

for week = 1:numweeks
    weekly_data = load(strcat(store_address,'\Metrics\Week',num2str(week),'_metrics.mat'));
    data{week} = weekly_data.metric_data;
    cells = zeros(1,numweeks);
    alive = zeros(1,numweeks);
    perc = zeros(1,numweeks);
    for m = 1:length(data{week})
        cells(m) = length(data{week}{m}.Centroids);
        alive(m) = sum(data{week}{m}.WHratios > 2.5);
        perc(m) = data{week}{m}.ScaffPerc;
    end
    
    numcells{week} = cells;
    numalive{week} = alive;
    scaffperc{week} = perc;    
    numcells_lim(week) = max(cells);
    numalive_lim(week) = sum(data{week}{1}.WHratios > 2.5);
    scaffperc_lim(week) = max(perc);
end
YLim1 = max(numcells_lim)*1.05;
YLim2 = max(numalive_lim)*1.05;
YLim3 = max(scaffperc_lim)*1.05;
XLim = 0;
for week = 1:numweeks
    x_limit = length(numcells{week})*str2double(thicknesses{week});
    if(XLim < x_limit)
        XLim = x_limit;
    end
end
XLim_max = XLim;
numweeks = 3;
f = figure('Renderer', 'painters', 'Position', [300 100 900 600]);
for week = 1:numweeks  
    distance = (1:length(numcells{week})).*str2double(thicknesses{week});
    
    %Num of cells
    subplot(3,numweeks,week)
    bar(distance,numcells{week});
    cell_bar(week) = gca;
    if(week == 2)
        title({"Week " + num2str(week), 'Number of Cells'});
    else
        title({"Week " + num2str(week),''});
    end
    if(week == 1)
        ylabel('Number of Cells');
    end
    ylim([0 YLim1]);
    yticks([0:500:YLim1]);
    xlim([0 XLim]);
       
    %Num alive cells
    ff2 = subplot(3,numweeks,numweeks + week);
    ff2.Position = ff2.Position + [0, 0.02,0,0];%0.1300    0.42    0.775    0.06];  
    bar(distance,numalive{week});
    alive_bar(week) = gca;
    if(week == 2)
        title('Number of Alive Cells');
    end
    if(week == 1)
     ylabel('Number of alive cells');
    end
    ylim([1 YLim2]);
    xlim([0 XLim]);
    
    %Scaff perc
    ff3 = subplot(3,numweeks,numweeks*2 + week);
    bar(distance,scaffperc{week});
    ff3.Position = ff3.Position + [0, 0.04,0,0];%0.1300    0.42    0.775    0.06]; 
    scaff_bar(week) = gca;
    if(week == 2)
        title('% of Scaffold Occupied');
    end
    if(week == 1)
        ylabel('% of scaffold occupied');
    end
    xlabel('Distance into sample / \mu{\itm}');  
    ylim([1 YLim3]);
    xlim([0 XLim]);
    ax2(week) = axes('XLim',[0 length(numcells{week})],'XTick',0:length(numcells{week})/5:length(numcells{week}),'Position',get(scaff_bar(week),'Position') + [0 -0.06 0 0]...
       ,'XAxisLocation','bottom','color','none');
    ax2(week).YAxis.Visible = 'off';
    xlabel('Slice');
    
end
th = 2.5;
pos=get(alive_bar(end),'position');
t = annotation('textbox', [pos(1)+pos(3)-0.005 pos(2)+0.09 0.1 pos(4)-0.05], 'string', {'Width to Height Ratio Threshold:',th}, 'FontSize', 8, 'FontWeight','bold'...
    ,'EdgeColor','none','HorizontalAlignment','center');
S=['[th,alive_bar] = new_alive_th(get(gco,''value''),data,thicknesses,XLim,t,alive_bar);'];  
alive_slider=uicontrol('style','slider','units','normalized','position',[pos(1)+pos(3)+0.037 pos(2)-0.025 0.017 pos(4)-0.05],'callback',S,'min',0,'max',6,'value',th);

pos2 = get(scaff_bar(2),'position');
t2 = annotation('textbox', [pos2(1)-0.155 0.038 0.2 0.001], 'string','X-axis limit:', 'FontSize', 8, 'FontWeight','bold'...
     ,'EdgeColor','none','HorizontalAlignment','center');
t3 = annotation('textbox', [pos2(1)+ 0.157 0.038 0.2 0.001], 'string',num2str(XLim) + " \bf\mu\itm", 'FontSize', 8, 'FontWeight','bold'...
     ,'EdgeColor','none','HorizontalAlignment','center');
S2= ['XLim = change_xlimit(get(gco,''value''),cell_bar,alive_bar,scaff_bar,numweeks,thicknesses,ax2,t3, scaffperc);'];  
xlim_slider_pos = [pos2(1) pos(2)-0.423 pos2(3) 0.025];
xlim_slider=uicontrol('style','slider','units','normalized','position',xlim_slider_pos,...
    'callback',S2,'min',1,'max',XLim,'value',XLim);

S3= 'saveas(f, strcat(store_address,''\Plots\results_xlim_'',num2str(XLim),''_wtoh_'',num2str(th),''.jpg''));';  
save_button=uicontrol('style','pushbutton','String','Save Figure','position',[820 4 80 20],...
     'callback',S3);

%% 3D project data
command = 'C:\Users\artur\Desktop\4thYear\Project\Fiji.app\ImageJ-win64.exe -macro C:/Users/artur/Desktop/4thYear/Project/Project-code/3Dview.ijm';
system(command);

%% Scaled Numcells
h = figure('Renderer', 'painters', 'Position', [50 50 1200 600]);
a = 1;
for k = [0,5,10] 
    subplot(1,3,a);
    hold on;
    for week = 1:numweeks
        x_values = (1:length(numcells{week}))/length(numcells{week});
        y_values = numcells{week} / numcells_lim(week);
        if(k)
            y_values = movmean(y_values,k);
        end
        plot(x_values,y_values,'LineWidth',1.25);
        if(k==5)
            title({'Number of cell plots scaled and moving averaged filtered for k:',strcat("k = ",num2str(k))});
        else
            title(strcat("k = ",num2str(k)));
        end
        legend({'Week 1', 'Week 2', 'Week 3'});
        xlabel('Scaled distance into sample');
        ylabel('Scaled number of cells');
    end
    hold off;
    a = a + 1;
end
saveas(h, strcat(store_address,'\Plots\scaled_num_cells.jpg'));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Average Scaffold Pore Size

image = imread("D:/Artur/RawData/Fib_Week3/Ch1_Stitched_Sections/Stitched_Z001.tif");
image = image > 0;
se = strel('square',10);
image = imdilate(imerode(image,se),se);

props = regionprops(image,'Area');
image2 = zeros(size(image,1),size(image,2));
m = 700;
finished = false;
while(~finished)
    props = cell2mat(struct2cell(regionprops(image2,'Area')));
    se = strel('square',m);
    image2 = imerode(imdilate(image,se),se); %Closing operation
    subplot(1,2,1);
    imshow(image);
    subplot(1,2,2);
    imshow(image2);
    title("m = " + num2str(m) + ", numobjects = " + num2str(length(props)));
    m = m+20;
    if(sum(props>200000) == length(props))
        finished = true;
    end
    pause(0.1);  
end
  
%% Calculate metrics if images already processed
bar = waitbar(0, {'Calculating metrics', 'Estimated time remaining: ...'});
address = 'C:\users\Artur\Desktop\newtemp';
mkdir(strcat(address,'\Metrics')); %Create folder where to store calculated metrics
slices_completed = 0;
time_per_slice = [];

for week = 2:3
    scafffiles = dir([strcat(address,'\Scaffolds\Week',num2str(week)),'\*.png']);
    cellfiles = dir([strcat(address,'\Cells\Week',num2str(week)),'\*.png']);
    numslices = numel(cellfiles);
    metric_data = cell(1,numslices); %Pre allocating for speed   
    
    for slice = 1:numslices
        tic
        cell_img = imread(strcat(cellfiles(slice).folder,'\',cellfiles(slice).name));
        cell_img = double(cell_img>0);
        
        %%% Get scaffold image and convert it to binary  
        scaff_img = imread(scafffiles(slice).folder + "\" + scafffiles(slice).name);
        scaff_img = double(scaff_img > 0);
        
    %%% Calculate Metrics  
        cc = bwconncomp(cell_img);
        slicedata.Centroids = cell2mat(struct2cell(regionprops(cc,'Centroid'))'); %Centroids of all cells
        slicedata.Areas = cell2mat(struct2cell(regionprops(cc,'Area'))'); %Areas of all cells
        slicedata.BBoxes = cell2mat(struct2cell(regionprops(cc,'BoundingBox'))'); %Bounding boxes of all cells
        MinAx = cell2mat(struct2cell(regionprops(cc,'MinorAxisLength'))');
        MajAx = cell2mat(struct2cell(regionprops(cc,'MajorAxisLength'))');
        slicedata.WHratios = MajAx ./ MinAx; %Width to height ratios of cell
        
        % Calculate percentage of scaffold occupied by cells    
        totalscaff = sum(scaff_img,'all');
        slicedata.ScaffPerc = sum(cell_img + scaff_img == 2,'all') / totalscaff * 100;

        metric_data{slice} = slicedata;
    
    %%% Update progress bar
        slices_completed = slices_completed + 1;
        time_per_slice(slices_completed) = toc;
        time_left = mean(time_per_slice)*(nr_total_slices - slices_completed)/60/60;
        hours_left = floor(time_left);
        min_left = round((time_left - hours_left) * 60);
        waitbar(slices_completed/nr_total_slices,bar, {strcat('Calculating metrics for week',' ',num2str(week)), strcat('Estimated time remaining: ',' ', num2str(hours_left),' h',' ',num2str(min_left),' min')});
    end
    save(strcat(store_address,'\Metrics\Week',num2str(week),'_metrics.mat'), 'metric_data');
end
delete(bar);

%%
address = 'D:\ARTUR\RawData\Fib_Week2\Ch2_Stitched_Sections_JPEG';
files = dir([address '\*.jpg']);
a = 1;
for file = files'
    disp(a);
    img = imread(strcat(address,'\',file.name));
    img = imresize(img,0.1);
    imwrite(img, strcat('C:\users\Artur\Desktop\cell\',file.name));
    a = a + 1;
end

%%
%img_adr = 'C:\users\Artur\Desktop\Project_data\Cells\Week2\Slice_001.png';
img_adr = 'D:\ARTUR\RawData\Fib_Week2\Ch2_Stitched_Sections\Stitched_Z001.tif';
img = imread(img_adr);
img = img(9000:13500, 1200:4800);
imshow(img);
imwrite(img, 'C:\users\Artur\Desktop\4thYear\Project\images\sliding\small.tif');
img = img > 0;
imwrite(img, 'C:\users\Artur\Desktop\4thYear\Project\images\sliding\small_binary.png');
