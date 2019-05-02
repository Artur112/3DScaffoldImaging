%4th Year Biomedical Engineering Final Year Project, written by Artur
%Jurgenson 2018-2019.


%% Get address of folder where data is stored
disp('Please select folder where the weekly scans are stored (ie folders that contain the Ch_x scans for each week');
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

%% Segment the images with Neural Network
bar = waitbar(0, {'Segmenting cells from images', 'Estimated time remaining: ...'});

net = load('network1.mat'); %Load neural Network
net = net.net;
picsize = 300;  
weeknr = 1;
slices_completed = 0;
for week = week_folders'
    files = dir([strcat(week.folder,'\',week.name,'\','Ch2_Stitched_Sections') '\*.tif']);
    mkdir(strcat(store_address,'\Segmented_Cells\Week',num2str(weeknr)));
    slice = 1;  
    for file = files'
        tic;
        img = imread(strcat(file.folder,'\',file.name));
        img_segm = zeros(size(img,1),size(img,2));
        for ro = 1:floor(size(img,1)/picsize)
            for co = 1:floor(size(img,2)/picsize)
                section = img(1+picsize*(ro-1) : picsize*ro , 1+picsize*(co-1) : picsize*co);
                if(any(section,'all'))
                    shifts = [0,0,0,0];
                    if(ro > 1 && ro < 300 && co >1 && co <300)
                        while(any(section(:,1)))
                            shifts(1) = shifts(1) + 1;
                            section = img(1+picsize*(ro-1) : picsize*ro , 1+picsize*(co-1) - 5*shifts(1): picsize*co - 5*shifts(1));
                        end
                        while(any(section(:,picsize)))
                            shifts(2) = shifts(2) + 1;
                            section = img(1+picsize*(ro-1) : picsize*ro , 1+picsize*(co-1) + 5*shifts(2): picsize*co + 5*shifts(2));
                        end
                        while(any(section(1,:)))
                            shifts(3) = shifts(3) + 1;
                            section = img(1+picsize*(ro-1) - 5*shifts(3) : picsize*ro - 5*shifts(3), 1+picsize*(co-1): picsize*co);
                        end
                        while(any(section(picsize,:)))
                            shifts(4) = shifts(4) + 1;
                            section = img(1+picsize*(ro-1) + 5*shifts(4): picsize*ro + 5*shifts(4) , 1+picsize*(co-1): picsize*co);
                        end
                    end
                    C = semanticseg(section,net);
                    if(~any(shifts))
                        img_segm(1+picsize*(ro-1) : picsize*ro , 1+picsize*(co-1) : picsize*co) = double(C) - 1;
                    else
                        img_segm(1+picsize*(ro-1) + 5*(shifts(4)-shifts(3)): picsize*ro + 5*(shifts(4)-shifts(3)), 1+picsize*(co-1) +5*(shifts(2)-shifts(1)) : picsize*co +5*(shifts(2)-shifts(1))) = double(C) - 1;
                    end
                end
            end
        end
        if(slice < 10)
            imwrite(img_segm,strcat(store_address,'\Segmented_Cells\Week',num2str(weeknr),'\Slice_00',num2str(slice),'.png'));
        elseif(slice<100)
            imwrite(img_segm,strcat(store_address,'\Segmented_Cells\Week',num2str(weeknr),'\Slice_0',num2str(slice),'.png'));
        else
            imwrite(img_segm,strcat(store_address,'\Segmented_Cells\Week',num2str(weeknr),'\Slice_',num2str(slice),'.png'));
        end
        time_per_slice = toc;
        time_left = time_per_slice*(nr_total_slices - slices_completed)/60/60;
        hours_left = floor(time_left);
        min_left = round((time_left - hours_left) * 60);
        waitbar(slices_completed/nr_total_slices,bar, {'Segmenting cells from images', strcat('Estimated time remaining: ', num2str(hours_left),' h',' ',num2str(min_left),' min')});
        slice = slice + 1;
        slices_completed = slices_completed + 1;
    end
    weeknr = weeknr + 1;
end
delete(bar);

%%  Convert Scaffold Images to Binary Images
bar = waitbar(0, 'Converting Scaffold Images to Binary');
weeknr = 1;
slices_completed = 0;
for week = week_folders'
    files = dir([strcat(week.folder,'\',week.name,'\','Ch1_Stitched_Sections') '\*.tif']);
    mkdir(strcat(store_address,'\Scaffold_Binary\Week',num2str(weeknr)));
    slice = 1;
    for file = files'
        disp(slice);
        image = imread(file.folder + "\" + file.name);
        image = double(image > 0);
        if(slice<10)
            imwrite(image,strcat(store_address,'\Scaffold_Binary\Week',num2str(weeknr),'\Slice_00', num2str(slice)+".png"));
        elseif(slice<100)
            imwrite(image,strcat(store_address,'\Scaffold_Binary\Week',num2str(weeknr),'\Slice_0', num2str(slice)+".png"));
        else
            imwrite(image,strcat(store_address,'\Scaffold_Binary\Week',num2str(weeknr),'\Slice_', num2str(slice)+".png"));
        end
        slice = slice + 1;
        slices_completed = slices_completed + 1;
        waitbar(slices_completed/nr_total_slices);
    end
    weeknr = weeknr + 1;
end
delete(bar);
%% Calculate Metrics
mkdir(strcat(store_address,'\Metrics'));
for week = 1:numweeks
    disp(week);
    cellfiles = dir([strcat(store_address,'\Segmented_Cells\Week',num2str(week)) '\*.png']);
    scafffiles = dir([strcat(store_address,'\Scaffold_Binary\Week',num2str(week)) '\*.png']);
    numslices = numel(cellfiles);
    metric_data = cell(1,numslices); %Pre allocating for speed
    
    for slice = 1:numslices
        % Calculate metrics of the cells
        cell_image = imread(cellfiles(slice).folder + "\" + cellfiles(slice).name);
        cell_image = im2double(cell_image);
        cc = bwconncomp(cell_image);
        slicedata.Centroids = cell2mat(struct2cell(regionprops(cc,'Centroid'))'); %Centroids of all cells
        slicedata.Areas = cell2mat(struct2cell(regionprops(cc,'Area'))'); %Areas of all cells
        slicedata.BBoxes = cell2mat(struct2cell(regionprops(cc,'BoundingBox'))'); %Bounding boxes of all cells
        MinAx = cell2mat(struct2cell(regionprops(cc,'MinorAxisLength'))');
        MajAx = cell2mat(struct2cell(regionprops(cc,'MajorAxisLength'))');
        slicedata.WHratios = MajAx ./ MinAx; %Width to height ratios of cells

        % Calculate percentage of scaffold occupied by cells    
        scaffold_image = imread(scafffiles(slice).folder + "\" + scafffiles(slice).name);
        scaffold_image = double(scaffold_image > 0);
        totalscaff = sum(scaffold_image,'all');
        slicedata.ScaffPerc = sum(cell_image + scaffold_image == 2,'all') / totalscaff * 100;

        %Store data for current slice and continue
        metric_data{slice} = slicedata;
    end
    save(strcat(store_address,'\Metrics\Week',num2str(week),'_metrics.mat'), 'metric_data');
end

%% Result Visualization
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
        
f = figure('Renderer', 'painters', 'Position', [50 50 900 600]);
%ha = tight_subplot(3,2,[.01 .03],[.1 .01],[.01 .01])
for week = 1:numweeks  
    if week == 1
        distz  = 50;
    else
        distz = 70;
    end
    distance = (1:length(numcells{week})).*distz;
    
    %Num of cells
    subplot(3,3,week)
    plot(distance,numcells{week});%length(metric_data{m}.Centroids));
    if(week == 2)
        title({"Week " + num2str(week), 'Number of Cells'});
    else
        title({"Week " + num2str(week),''});
    end
    hold off;
    xlabel('Distance into sample / um');    
    ylabel('Number of Cells');
    ylim([0 YLim1]);
    xlim([0 distance(end)]);
       
    %Num alive cells
    ff2 = subplot(3,3,3 + week);
    %ff2.Position = ff2.Position + [0, 0.152,0,0];%0.1300    0.42    0.775    0.06];  
    plot(distance,numalive{week});
    if(week == 2)
        title('Number of Alive Cells');
    end
    ylim([1 YLim2]);
    xlim([0 distance(end)]);
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
    xlim([0 distance(end)]);
    ylabel('% of Scaffold covered');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% ------------ Average Scaffold Pore Size

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
    numcells = [];
    for m = 1:50
        numcells = [numcells,length(metric_data{m}.Centroids)];    
    end
    f = figure;
    subplot(2,1,1);
    bar(numcells);
    hold on;
  %  plot(1:50,numcells/sum(numcells)*100,'LineWidth',1.5);
    matrix = numcells/sum(numcells)*100;
    ff = subplot(2,1,2, 'Parent',f);
    ff.Position = [0.1300    0.42    0.775    0.06];   
    imagesc(matrix); 
    colormap('autumn');
    set(gca, 'YTick', []);
    %%
    se = strel('square',500);
    image2 = imerode(imdilate(image,se),se); 
%     image3 = regionprops(image2,'FilledImage');
%     image3 = image3.FilledImage;
%     figure;
%     subplot(1,2,1);
%     imshow(image2);
%     subplot(1,2,2);
%     imshow(image3);
%     image4 = edge(image3);
%     figure;
%     i
%     imshow(image4);
    props = regionprops(image2,'Area');
    tic
    C = imfuse(image,image2,'ColorChannels',[1 2 0]);
    toc/60
    imshow(C);
    title("numobjects = " + num2str(length(props)));
    
    %%
    areascaffold = 0;
    for m = 1:size(image,1)
        for n = 1:size(image,2)
            if(image(m,n) && image2(m,n))
                areascaffold = areascaffold + 1;
            end
        end
    end
    props = regionprops(image2,'Area');
    pore_size = areascaffold
    
    %% 
   
            
%     if(sl<10)
%             movefile(address + "Slice_"+num2str(sl)+".png",address + "Slice_00"+num2str(sl)+".png");
%         elseif(sl<100)
%             movefile(address + "Slice_"+num2str(sl)+".png",address + "Slice_0"+num2str(sl)+".png");
%         end     

    %%
    filename = 'scaffold.gif';
    
    for n=1:50
        if(n<10)
            im = imread("C:\users\Artur\Desktop\temp\ScaffoldPics\Week3\Slice_00" + num2str(n) + ".png");
        else
            im = imread("C:\users\Artur\Desktop\temp\ScaffoldPics\Week3\Slice_0" + num2str(n) + ".png");
        end
        imgs{n} = im;
    end  
        %%
    address = 'C:\users\Artur\Desktop\temp\ScaffoldPics\Week2';
    files = dir([address '\*.png']);
    %all = zeros(20876*0.1,17250*0.1,60);
    sl = 1;
    tic;
    for img = files'
        imgg =  imresize(imread(strcat(strcat(address, "\"),img.name)),0.05) > 10;
        str = strel('diamond',10);
        imgg = imerode(imgg,str);
        imgg = imdilate(imgg,str);
        all(:,:,sl) = imresize(imread(strcat(strcat(address, "\"),img.name)),0.05) > 10;
        sl = sl + 1;
        display(sl);
    end
    toc/60
    img
%% 
%all(:,:,
  for n=1:50
        if(n<10)
            im = imread("C:\users\Artur\Desktop\temp\ScaffoldPics\Week3\Slice_00" + num2str(n) + ".png");
        else
            im = imread("C:\users\Artur\Desktop\temp\ScaffoldPics\Week3\Slice_0" + num2str(n) + ".png");
        end
        imgs{n} = im;
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
  
  %%
  address = 'C:\users\Artur\Desktop\temp2';
  files = dir([address '\*.png']);
  a = 1;
  for file = files'
      disp(a);
      img = imread(address + "\" + file.name);
      img = imresize(img,0.1);
      imwrite(img,"C:\users\Artur\Desktop\temp3\" + file.name);
      a = a + 1;
  end
 %%
 disp('c');
  for m = 1:5
      fprintf('aa');
  end
 % fprintf(repmat2)); %Deletes 'world' and '\n
  %%
  fprintf('Hello World\n');
pause
fprintf(repmat('\b', 1, 5)); %Deletes 'world' and '\n
%%
    
%     %Bars
%     ff1 = subplot(4,3,3 + week, 'Parent',f);
%     numcells = [];
%     for m = 1:length(metric_data)
%         numcells = [numcells,length(metric_data{m}.Centroids)];    
%     end
%     numcells = numcells/sum(numcells)*100;
%     imagesc(numcells); 
%     colormap('autumn');
%     ff1.Position = ff1.Position + [0,0.148,0,-0.14];%0.1300    0.42    0.775    0.06];   
%     set(gca, 'YTick', []);
%     set(gca,'XTick',[]);