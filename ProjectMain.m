clc; clear all; close all;
%% Section for loading in and visualizing individual images
%Load in the image - analysing week 2 files now
%img_scaffold = imread('D:\ARTUR\Fib_Week2\Ch1_Stitched_Sections\Stitched_Z030.tif');



img_cells = imread('D:\ARTUR\Fib_Week2\Ch2_Stitched_Sections\Stitched_Z030.tif');
picsize = 300;
img = 1;
for ro = 1:floor(size(img_cells,1)/picsize)
    for co = 1:floor(size(img_cells,2)/picsize)
        props = regionprops(img_cells(picsize*(ro-1)+1 : picsize*ro,picsize*(co-1)+1:picsize*co));
        if(max([props.Area]) > 50)
            smallerpics(:,:,img) = img_cells(picsize*(ro-1)+1 : picsize*ro,picsize*(co-1)+1:picsize*co);
            img = img + 1;
        end
    end
end

for img = 1:10
    subplot(2,5,img);
    imshow(smallerpics(:,:,200+img));
end
%img_nucleus = imread('D:\ARTUR\Fib_Week2\Ch3_Stitched_Sections\Stitched_Z030.tif');
imshow(img_cells)
%imgbin_scaffold = img_scaffold > 0;
%imgbin_cells = img_cells > 0;
%imgbin_nucleus = img_nucleus > 0;
%overlay_im = cat(3, imgbin_scaffold, imgbin_cells, imgbin_nucleus);
%C = imfuse(imgbin_scaffold,imgbin_nucleus);
%imshow(C);
%imagesc(image_scaffold); colormap('gray');

%D = -bwdist(~imgbin_cells);
%D(~imgbin_cells) = -Inf;
%I2 = imhmin(D,3);
%L = watershed(I2);
%imshow(label2rgb(L,'jet','w'));  title('Watershed');

%%
%Convert image to binary image
imgbin_scaffold = img_scaffold > 0;
imshow(imgbin_scaffold); title("Binary Scaffold Channel");
imgbin_cells = img_cells > 0;
figure; imshow(imgbin_cells); title("Binary Cell Channel");
%imagesc(image_scaffold_Binary); colormap('gray');
%clear img_scaffold;
%clear img_cells;

% %Split image in half since the JPEG scans are of two scaffolds
% image_scaffold_Binary1 = image_scaffold_Binary(:,1:round(length(image_scaffold_Binary)/2));
% image_scaffold_Binary2 = image_scaffold_Binary(:,round(length(image_scaffold_Binary)/2)+1 : length(image_scaffold_Binary));

%Gonna analyze the right side for now
%% Section for performing dilation and erosion for only considering scaffold area for segmentation
%Structuring elements
% se = strel('square',8);
% se2 = strel('square',30);
% se3 = strel('square',1300);
% 
% %Eroding to get rid of small things around
% erodeimg1 = imerode(image_scaffold_Binary,se);
% dilateimg1 = imdilate(erodeimg1,se2);
% erodeimg2 = imerode(dilateimg1,se3);
%finalimg = imerode(dilateimg2,se3);

% subplot(1,5,1); imshow(image_scaffold_Binary); title("Original Binary Image");
% subplot(1,5,2); imshow(erodeimg1); title("Eroded with 15x15 square");
% subplot(1,5,3); imshow(dilateimg1); title("Dilated with 7000 square");
% subplot(1,5,4); imshow(erodeimg2); title("Erode with 500 square");
%subplot(1,5,5); imshow(finalimg); title("Eroded with 300x300 square");

%subplot(1,4,1);
C = imfuse(imgbin_scaffold,imgbin_cells);
imshow(C);
title('Original Scaffold with cells overlayed');

% figure;
% C2 = imfuse(img_scaffold+50,imgbin_cells);
% imshow(C2);
% title('Modified Scaffold with cells overlayed');

%% Section for experimenting with segmentation methods

area = getfield(regionprops(finalimg,'BoundingBox'),'BoundingBox');
smallerimageglob = finalimg(round(area(2)):round(area(2)+area(4)),round(area(1)):round(area(1)+area(3)));
smallerimage = image_scaffold_Binary2(round(area(2)):round(area(2)+area(4)),round(area(1)):round(area(1)+area(3)));
subplot(2,5,6);imshow(smallerimageglob); title("Bounding box area");
subplot(2,5,7); imshow(smallerimage); title("Same area in original pic");


for m = 1:length(smallerimage(:,1))
    for n = 1:length(smallerimage(1,:))
        if(smallerimageglob(m,n) == 0)
            smallerimage(m,n) = 0;
        end
    end
end

subplot(2,5,8); imshow(smallerimage); title("Pixels values outside scaffold removed");
% justboundary = zeros(5327,9185);
% for row = 1:length(boundary)
%     justboundary(boundary(row,1),boundary(row,2)) = 1;
% end
% justboundary_thickened = imdilate(justboundary,[1, 1 ,1;1,1, 1;1,1,1]);
% subplot(1,4,3); imshow(justboundary_thickened); title("Boundary");
% 
% C = imfuse(finalimg,imageDataBinary,'blend');
% subplot(1,4,3); imshow(C); title("Images overlayed");
% B = bwboundaries(finalimg);
% boundary = B{1};

jpgFileName = '..\Raw Data\Ch2_Stitched_Sections_JPEG\Ch2_Stitched_Sections_JPEG\Stitched_Z001.jpg';
image_cells = imread(jpgFileName);

% figure;
% imshow(image_cells);

image_cellsHalf = image_cells(:,round(length(image_scaffold_Binary)/2)+1 : length(image_scaffold_Binary));
image_cellsSmall = image_cellsHalf(round(area(2)):round(area(2)+area(4)),round(area(1)):round(area(1)+area(3)));


for m = 1:length(smallerimage(:,1))
    for n = 1:length(smallerimage(1,:))
        if(smallerimageglob(m,n) == 0)
            image_cellsSmall(m,n) = 0;
        end
    end
end
% 
% for m = 1:length(smallerimage(:,1))
%     for n = 1:length(smallerimage(1,:))
%         if(image_cellsSmall(m,n) ~= 0)
%             image_cellsSmall(m,n) = image_cellsSmall(m,n) + 50;
%         end
%     end
% end

figure;
subplot(1,2,1);
imshow(image_cellsSmall); colormap('gray'); title("Cells in the scaffold area only Original");
%Convert image to binary image
image_cells_Binary = image_cells > 0;

%Split image in half since the JPE scans are of two scaffolds
image_cells_Binary1 = image_cells_Binary(:,1:round(length(image_cells_Binary)/2));
image_cells_Binary2 = image_cells_Binary(:,round(length(image_cells_Binary)/2)+1 : length(image_cells_Binary));
image_cells_Binary2 = image_cells_Binary2(round(area(2)):round(area(2)+area(4)),round(area(1)):round(area(1)+area(3)));
for m = 1:length(smallerimage(:,1))
    for n = 1:length(smallerimage(1,:))
        if(smallerimageglob(m,n) == 0)
            image_cells_Binary2(m,n) = 0;
        end
    end
end

subplot(1,2,2);
imshow(image_cells_Binary2); title("Cells in the scaffold area only Binary");
regioncells = regionprops(bwconncomp(image_cells_Binary2));
figure;
imshow(label2rgb(bwlabel(image_cells_Binary2))); title('Connected components');
cellareas = zeros(length(regioncells),1);
for n = 1: length(regioncells)
    cellareas(n) = getfield(regioncells(n),'Area');
end
figure;

hist(cellareas',10);


figure;
imshow(image_cells_Binary2); title('Bin');

figure;
D = -bwdist(~image_cells_Binary2);
D(~image_cells_Binary2) = -Inf;
I2 = imhmin(D,3);
L = watershed(I2);
imshow(label2rgb(L,'jet','w'));  title('Watershed');

regionprops(L);

%% Section for performing analysis on all images on all images at once.
%Reading in all the image files 
for k = 1:157
    if(k<10)
        jpgFileName = strcat('..\Raw Data\Ch2_Stitched_Sections_JPEG\Ch2_Stitched_Sections_JPEG\Stitched_Z00', num2str(k), '.jpg');
        if exist(jpgFileName, 'file')
            imageData(:,:,k) = imread(jpgFileName);
        else
            fprintf('File %s does not exist.\n', jpgFileName);
        end
    elseif(k>99)
        jpgFileName = strcat('..\Raw Data\Ch2_Stitched_Sections_JPEG\Ch2_Stitched_Sections_JPEG\Stitched_Z', num2str(k), '.jpg');
        if exist(jpgFileName, 'file')
            imageData(:,:,k) = imread(jpgFileName);
        else
            fprintf('File %s does not exist.\n', jpgFileName);
        end
    else
        jpgFileName = strcat('..\Raw Data\Ch2_Stitched_Sections_JPEG\Ch2_Stitched_Sections_JPEG\Stitched_Z0', num2str(k), '.jpg');
        if exist(jpgFileName, 'file')
            imageData(:,:,k) = imread(jpgFileName);
        else
            fprintf('File %s does not exist.\n', jpgFileName);
        end
    end
end

%Converting all the image files to binary
imageDataBinary = imageData > 0;

%Keeping only second half of each image - must preallocate first for speed
imageDataBinary1 = false(5327,4593,157);
imageDataBinary2 = false(5327,4592,157);
for pic = 1:length(imageDataBinary(1,1,:))
    imageDataBinary1(:,:,pic) = imageDataBinary(:,1:round(length(imageDataBinary)/2),pic);
    imageDataBinary2(:,:,pic) = imageDataBinary(:,round(length(imageDataBinary)/2)+1 : length(imageDataBinary),pic);      
end

for pic = 1:length(imageDataBinary(1,1,:))
  imageDataBinary1small(:,:,pic) = imageDataBinary1(round(area(2)):round(area(2)+area(4)),round(area(1)):round(area(1)+area(3)),pic);
  imageDataBinary2small(:,:,pic) = imageDataBinary2(round(area(2)):round(area(2)+area(4)),round(area(1)):round(area(1)+area(3)),pic);
end  
  
for m = 1:length(smallerimage(:,1))
    for n = 1:length(smallerimage(1,:))
        if(smallerimageglob(m,n) == 0)
            imageDataBinary2small(m,n,:) = 0;
        end
    end
end


% image1CC = bwconncomp(image1);
% image1CC_centroids = regionprops(image1CC,'Centroid');

% cellpoints = double(zeros(68344,3));
cellcounter = 1;
for pic = 1:157
    tempregionprops = regionprops(bwconncomp(imageDataBinary2small(:,:,pic)));
    for n = 1:length(tempregionprops)
        temp = getfield(tempregionprops(n),'Centroid');
        cellpoints(cellcounter,1) = temp(1);
        cellpoints(cellcounter,2) = temp(2);
        cellpoints(cellcounter,3) = pic;
        cellcounter = cellcounter + 1;
    end
end

figure;
plot3(cellpoints(:,1),cellpoints(:,2),cellpoints(:,3),'.'); 