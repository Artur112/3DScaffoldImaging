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
for m = 1:numweeks
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
mkdir(strcat(store_address,'\Metrics')); %Create folder where to store calculated metrics
weeknr = 1;
slices_completed = 0;

for week = week_folders'
    weeknr
    mkdir(strcat(store_address,'\Cells\Week',num2str(weeknr))); %Create folder to store cell images in
    mkdir(strcat(store_address,'\Scaffolds\Week',num2str(weeknr))); %Create folder to store scaffold images in
    scafffiles = dir([strcat(week.folder,'\',week.name,'\','Ch1_Stitched_Sections') '\*.tif']);
    cellfiles = dir([strcat(week.folder,'\',week.name,'\','Ch2_Stitched_Sections') '\*.tif']);
    numslices = numel(cellfiles);
    metric_data = cell(1,numslices); %Pre allocating for speed
    
    for slice = 1:1
        tic;
    %%% Segment cell image by adaptively sliding a 300x300 section along the image that is segmented with neural network
        cell_img = imread(strcat(cellfiles(slice).folder,'\',cellfiles(slice).name));
        cell_img_segm = zeros(size(cell_img,1),size(cell_img,2));
        for ro = 1:floor(size(cell_img,1)/picsize)
            for co = 1:floor(size(cell_img,2)/picsize)
                section = cell_img(1+picsize*(ro-1) : picsize*ro , 1+picsize*(co-1) : picsize*co);
                if(any(section,'all'))
                    shifts = [0,0,0,0];
                    if(ro > 1 && ro < 300 && co >1 && co <300)
                        while(any(section(:,1)))
                            shifts(1) = shifts(1) + 1;
                            section = cell_img(1+picsize*(ro-1) : picsize*ro , 1+picsize*(co-1) - 5*shifts(1): picsize*co - 5*shifts(1));
                        end
                        while(any(section(:,picsize)))
                            shifts(2) = shifts(2) + 1;
                            section = cell_img(1+picsize*(ro-1) : picsize*ro , 1+picsize*(co-1) + 5*shifts(2): picsize*co + 5*shifts(2));
                        end
                        while(any(section(1,:)))
                            shifts(3) = shifts(3) + 1;
                            section = cell_img(1+picsize*(ro-1) - 5*shifts(3) : picsize*ro - 5*shifts(3), 1+picsize*(co-1): picsize*co);
                        end
                        while(any(section(picsize,:)))
                            shifts(4) = shifts(4) + 1;
                            section = cell_img(1+picsize*(ro-1) + 5*shifts(4): picsize*ro + 5*shifts(4) , 1+picsize*(co-1): picsize*co);
                        end
                    end
                    C = semanticseg(section,net);
                    if(~any(shifts))
                        cell_img_segm(1+picsize*(ro-1) : picsize*ro , 1+picsize*(co-1) : picsize*co) = double(C) - 1;
                    else
                        cell_img_segm(1+picsize*(ro-1) + 5*(shifts(4)-shifts(3)): picsize*ro + 5*(shifts(4)-shifts(3)), 1+picsize*(co-1) +5*(shifts(2)-shifts(1)) : picsize*co +5*(shifts(2)-shifts(1))) = double(C) - 1;
                    end
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
        slicedata.BBoxes = cell2mat(struct2cell(regionprops(cc,'BoundingBox'))'); %Bounding boxes of all cells
        MinAx = cell2mat(struct2cell(regionprops(cc,'MinorAxisLength'))');
        MajAx = cell2mat(struct2cell(regionprops(cc,'MajorAxisLength'))');
        slicedata.WHratios = MajAx ./ MinAx; %Width to height ratios of cells
        
        toremove = find(slicedata.Areas < 50); %Find areas that are smaller than 50 (these are
        %are cell fragments that count as noise or cells that partially show up from deeper layers.
        
        %Remove these areas from the cell image
        pxls = regionprops(cc,'PixelList');
        for ii = 1:length(toremove)
            cell_img_segm(pxls(toremove(ii)).PixelList(:,2), pxls(toremove(ii)).PixelList(:,1)) = 0;
        end
        
        %Remove these metrics from the metric_data
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
            imwrite(scaff_img,strcat(store_address,'\Scaffolds\Week',num2str(weeknr),'\Slice_00',num2str(slice),'.png'));
        elseif(slice<100)
            imwrite(cell_img_segm,strcat(store_address,'\Cells\Week',num2str(weeknr),'\Slice_0',num2str(slice),'.png'));
            imwrite(scaff_img,strcat(store_address,'\Scaffolds\Week',num2str(weeknr),'\Slice_0',num2str(slice),'.png'));
        else
            imwrite(cell_img_segm,strcat(store_address,'\Cells\Week',num2str(weeknr),'\Slice_',num2str(slice),'.png'));
            imwrite(scaff_segm,strcat(store_address,'\Scaffolds\Week',num2str(weeknr),'\Slice_',num2str(slice),'.png'));
        end
    
    %%% Update bar
        slices_completed = slices_completed + 1;
        time_per_slice(slices_completed) = toc;
        time_left = mean(time_per_slice)*(nr_total_slices - slices_completed)/60/60;
        hours_left = floor(time_left);
        min_left = round((time_left - hours_left) * 60);
        waitbar(slices_completed/nr_total_slices,bar, {'Segmenting cells from images', strcat('Estimated time remaining: ', num2str(hours_left),' h',' ',num2str(min_left),' min')});
    end
    save(strcat(store_address,'\Metrics\Week',num2str(weeknr),'_metrics.mat'), 'metric_data');
    weeknr = weeknr + 1;
end
delete(bar);
        
%% Visualize Results
mkdir(strcat(store_address,'\Plots'));
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
   
    numcells_lim(week) = length(data{week}{1}.Centroids);
    numalive_lim(week) = sum(data{week}{1}.WHratios > 2.5);
end
     
YLim1 = max(numcells_lim)+100;
YLim2 = max(numalive_lim)+100;
XLim = 0;
for week = 1:numweeks
    x_limit = length(numcells{week})*str2double(thicknesses{week});
    if(XLim < x_limit)
        XLim = x_limit;
    end
end
   
f = figure('Renderer', 'painters', 'Position', [50 50 900 600]);
for week = 1:numweeks  

    distance = (1:length(numcells{week})).*str2double(thicknesses{week});
    
    %Num of cells
    subplot(3,3,week)
    plot(distance,numcells{week});
    if(week == 2)
        title({"Week " + num2str(week), 'Number of Cells'});
    else
        title({"Week " + num2str(week),''});
    end
    hold off;
    xlabel('Distance into sample / um');    
    ylabel('Number of Cells');
    ylim([0 YLim1]);
    xlim([0 XLim]);
       
    %Num alive cells
    ff2 = subplot(3,3,3 + week);
    %ff2.Position = ff2.Position + [0, 0.152,0,0];%0.1300    0.42    0.775    0.06];  
    plot(distance,numalive{week});
    if(week == 2)
        title('Number of Alive Cells');
    end
    ylim([1 YLim2]);
    xlim([0 XLim]);
    xlabel('Distance into sample / um');    
    ylabel('Number of Alive Cells');
    
    %Scaff perc
    ff3 = subplot(3,3,6 + week);
    plot(distance,scaffperc{week});
   % ff3.Position = ff3.Position + [0, 0.135,0,0];%0.1300    0.42    0.775    0.06]; 
    if(week == 2)
        title('% of Scaffold Occupied');
    end
    xlabel('Slice');    
    xlim([0 XLim]);
    ylabel('% of Scaffold covered');
end
saveas(f, strcat(store_address,'\Plots\results.jpg'));
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
    
  %%
  
  for m = 1:30
  img = imresize(imgs{m},0.005); 
  [row,col] = find(img == 0);
  alphmap = ones(size(img,1),size(img,2))*0.8;
  for n = 1:length(row)
  alphmap(row(n),col(n)) = 0;
  end
  H(m) = slice(repmat(double(img),[1 1 10]),[],[], m); %slice() requires at least 2x2x2
  set(H(m),'EdgeColor','none', 'AlphaData', alphmap,'FaceAlpha','flat');
 
  if(m == 1)
      hold on;
  end
  end

    