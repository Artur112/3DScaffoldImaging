%4th Year Biomedical Engineering Final Year Project, written by Artur
%Jurgenson 2018-2019.

clc; clear all; close all;

%% Get address of folder where data is stored
disp('Please select folder where the data is (ie the Chx_Stitched_Sections folders)');
data_address = uigetdir('temp', 'Select folder where data is');

%% Segment the images with Neural Network
net = load('network1.mat'); %Load neural Network
net = net.net;
picsize = 300; 
files = dir([convertStringsToChars(strcat(data_address,"\Ch2_Stitched_Sections")) '\*.tif']);
slice = 1;  
tic;
for file = files'
    display("Segmenting Slice " + num2str(slice));
    img = imread(data_address + "\Ch2_Stitched_Sections\" + file.name);
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
    imwrite(img_segm,"C:\users\Artur\Desktop\temp\Slice_" + num2str(slice) + ".png");
    slice = slice + 1;
end
display("Segmenting all images took " + num2str(toc/60) + " minutes");

%% Calculate Metrics
segm_cell_imgs = dir(['C:\users\artur\Desktop\temp\Segmented' '\*.png']);
scaffold_imgs = dir(['D:\ARTUR\RawData\Fib_Week2\Ch1_Stitched_Sections' '\*.tif']);
numslices = numel(segm_cell_imgs);
metric_data = cell(numslices); %Pre allocating for speed
tic;
for slice = 1:numslices
    display("Analyzing Slice " + num2str(slice));
    
    % Calculate metrics of the cells
    cell_image = imread(segm_cell_imgs(1).folder + "\" + segm_cell_imgs(slice).name);
    cell_image = im2double(cell_image);
    cc = bwconncomp(cell_image);
    slicedata.Centroids = cell2mat(struct2cell(regionprops(cc,'Centroid'))'); %Centroids of all cells
    slicedata.Areas = cell2mat(struct2cell(regionprops(cc,'Area'))'); %Areas of all cells
    slicedata.BBoxes = cell2mat(struct2cell(regionprops(cc,'BoundingBox'))'); %Bounding boxes of all cells
    MinAx = cell2mat(struct2cell(regionprops(cc,'MinorAxisLength'))');
    MajAx = cell2mat(struct2cell(regionprops(cc,'MajorAxisLength'))');
    slicedata.WHRatios = MajAx ./ MinAx; %Width to height ratios of cells
    
    % Calculate percentage of scaffold occupied by cells    
    scaffold_image = imread(scaffold_imgs(1).folder + "\" + scaffold_imgs(slice).name);
    scaffold_image = double(scaffold_image > 0);
    totalscaff = sum(scaffold_image,'all');
    slicedata.scaffperc = sum(cell_image + scaffold_image == 2,'all') / totalscaff * 100
    
    %Store data for current slice and continue
    metric_data{slice} = slicedata;
end
display("Calculating metrics took " + num2str(toc/60) + " minutes");




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
    
%% Plot and Display Metrics
                
                
                
                