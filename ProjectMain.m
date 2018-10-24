%% Section for analyzing one image only
clc; clear all; close all;
%Load in the image
jpgFileName = '..\Raw Data\Ch1_Stitched_Sections_JPEG\Ch1_Stitched_Sections_JPEG\Stitched_Z001.jpg';
image_scaffold = imread(jpgFileName);

%Convert image to binary image
image_scaffold_Binary = image_scaffold > 0;

%Split image in half since the JPEG scans are of two scaffolds
image_scaffold_Binary1 = image_scaffold_Binary(:,1:round(length(image_scaffold_Binary)/2));
image_scaffold_Binary2 = image_scaffold_Binary(:,round(length(image_scaffold_Binary)/2)+1 : length(image_scaffold_Binary));

%Gonna analyze the right side for now

%Structuring elements
se = strel('square',15);
se2 = strel('square',400);
se3 = strel('square',300);

%Eroding to get ri
erodeimg1 = imerode(image_scaffold_Binary2,se);
dilateimg1 = imdilate(erodeimg1,se);
dilateimg2 = imdilate(dilateimg1,se2);
finalimg = imerode(dilateimg2,se3);

subplot(2,5,1); imshow(image_scaffold_Binary2); title("Original Binary Image");
subplot(2,5,2); imshow(erodeimg1); title("Eroded with 15x15 square");
subplot(2,5,3); imshow(dilateimg1); title("Dilated with 15x15 square");
subplot(2,5,4); imshow(dilateimg2); title("Dilated with 400x400 square");
subplot(2,5,5); imshow(finalimg); title("Eroded with 300x300 square");


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

%Convert image to binary image
image_cells_Binary = image_cells > 0;

%Split image in half since the JPEG scans are of two scaffolds
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

subplot(2,5,8);
imshow(image_cells_Binary2); title("Cells in the scaffold area only");


%% Section doing same analysis on all images at once.
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