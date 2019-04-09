%4th Year Biomedical Engineering Final Year Project, written by Artur
%Jurgenson 2018-2019.

clc; clear all; close all;

%% Get address of folder where data is stored
disp('Please select folder where the data is (ie the Chx_Stitched_Sections folders)');
address = uigetdir('temp', 'Select folder where data is');
a = dir([address '\*.tif']);

%% Segment the images with Neural Network
net = load('network1.mat');
net = net.net;

%imDir = 'C:\users\artur\Desktop\temp\';
% img = imread('C:\users\Artur\Desktop\slice001pic0006ro006co044.tif');
% C = semanticseg(img,net);
% B = labeloverlay(img,C);
% imshow(img);
%%
k = double(C) - 1;
r = imfuse(img,k,'falsecolor','Scaling','joint','ColorChannels',[1 2 0]);
imshow(img);
%%
picsize = 300;
slice = 1;
for file = a'
    tic;
    img = imread(address + "\" + file.name);
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
    imwrite(img_segm,'C:\users\Artur\Desktop\temp2\segmented.png');
    toc/60
end

%% Calculate Metrics
address = uigetdir('temp', 'Select folder where segmented full images are');
a = dir([address '\*.tif']);
tic;
b = 1;
total = numel(a);

for file = a(1) 
    image = imread(address + "\" + file.name);
    cc = bwconncomp(image);
    celldata.Centroid = cell2mat(struct2cell(regionprops(cc,'Centroid'))'); %Centroids of all cells
    celldata.Area = cell2mat(struct2cell(regionprops(cc,'Area'))'); %Areas of all cells
    celldata.BBox = cell2mat(struct2cell(regionprops(cc,'BoundingBox'))'); %Bounding boxes of all cells
    MinAx = cell2mat(struct2cell(regionprops(cc,'MinorAxisLength'))');
    MajAx = cell2mat(struct2cell(regionprops(cc,'MajorAxisLength'))');
    celldata.WHRatio = MajAx / MinAx; %Width to height ratios of cells
    display(num2str(b/total*100) + "% done");
    all_cells{b} = celldata;
    b = b + 1;
end
toc/60

%%
%---------------Percentage of scaffold occupied by cells - first version
address_scaff = 'D:\ARTUR\RawData\Fib_Week2\Ch1_Stitched_Sections';% uigetdir('temp', 'Select folder where scaffold images are');
address_cells = 'D:\ARTUR\Nuclei-Cell\5_Stitched\Week2_Nuclei\GroundTruth';% uigetdir('temp', 'Select folder where cell images are');
a = dir([address_scaff '\*.tif']);
a2 = dir([address_cells '\*.tif']);
total = numel(a);
pic = 1;
for pic = 1:5   
    imgscaff = imread(address_scaff + "/" + a(pic).name);
    imgcells = imread(address_cells + "/" + a2(pic).name);
    
    totalscaff = sum(imgscaff>0,'all');
    imgscaffbin = imgscaff > 0;
    imgcellsbin = imgcells > 0;
    sumcells = 0;
    for m = 1:size(imgcells,1)
        for n = 1:size(imgcells,2)
            if(imgcellsbin(m,n) && imgscaffbin(m,n))
                sumcells = sumcells + 1;
            end
        end
    end
    percentage(pic) = sumcells / totalscaff * 100
    pic = pic + 1
end

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
                
                
                
                