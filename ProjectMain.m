clc; clear all; close all;

%Iniatalize matrix to store pixel values into
%imageData = zeros(5327,9185,157,'uint8');
jpgFileName = '..\Raw Data\Ch1_Stitched_Sections_JPEG\Ch1_Stitched_Sections_JPEG\Stitched_Z001.jpg';
imageData = imread(jpgFileName);

imageDataBinary = imageData > 0;
% imagesc(imageData)
% title('Original Image');

%Structuring elements
se = strel('square',15);
se2 = strel('square',400);
se3 = strel('square',300);

%Eroding to get ri
erodeimg1 = imerode(imageDataBinary,se);
dilateimg1 = imdilate(erodeimg1,se);
dilateimg2 = imdilate(dilateimg1,se2);
finalimg = imerode(dilateimg2,se3);

subplot(1,4,1); imshow(imageDataBinary); title("Original BinaryJPEG");
subplot(1,4,2); imshow(finalimg); title("After erosion and dilation operations");

area = getfield(regionprops(finalimg,'BoundingBox'),'BoundingBox');
smallerimageglob = finalimg(round(area(2)):round(area(2)+area(4)),round(area(1)):round(area(1)+area(3)));
smallerimage = imageDataBinary(round(area(2)):round(area(2)+area(4)),round(area(1)):round(area(1)+area(3)));
subplot(1,4,3);imshow(smallerimage); title("Bounding box area");

for m = 1:length(smallerimage(:,1))
    for n = 1:length(smallerimage(1,:))
        if(smallerimageglob(m,n) == 0)
            smallerimage(m,n) = 0;
        end
    end
end

subplot(1,4,4); imshow(smallerimage); title("Pixels outside scaffold removed");
% justboundary = zeros(5327,9185);
% for row = 1:length(boundary)
%     justboundary(boundary(row,1),boundary(row,2)) = 1;
% end
% justboundary_thickened = imdilate(justboundary,[1, 1 ,1;1,1, 1;1,1,1]);
% subplot(1,4,3); imshow(justboundary_thickened); title("Boundary");

% C = imfuse(finalimg,imageDataBinary,'blend');
% subplot(1,4,3); imshow(C); title("Images overlayed");
% B = bwboundaries(finalimg);
% boundary = B{1};

%Reading in all the image files 
% for k = 1:157
%     if(k<10)
%         jpgFileName = strcat('..\Raw Data\Ch2_Stitched_Sections_JPEG\Ch2_Stitched_Sections_JPEG\Stitched_Z00', num2str(k), '.jpg');
%         if exist(jpgFileName, 'file')
%             imageData(:,:,k) = imread(jpgFileName);
%         else
%             fprintf('File %s does not exist.\n', jpgFileName);
%         end
%     elseif(k>99)
%         jpgFileName = strcat('..\Raw Data\Ch2_Stitched_Sections_JPEG\Ch2_Stitched_Sections_JPEG\Stitched_Z', num2str(k), '.jpg');
%         if exist(jpgFileName, 'file')
%             imageData(:,:,k) = imread(jpgFileName);
%         else
%             fprintf('File %s does not exist.\n', jpgFileName);
%         end
%     else
%         % Create an image filename, and read it in to a variable called imageData.
%         jpgFileName = strcat('..\Raw Data\Ch2_Stitched_Sections_JPEG\Ch2_Stitched_Sections_JPEG\Stitched_Z0', num2str(k), '.jpg');
%         if exist(jpgFileName, 'file')
%             imageData(:,:,k) = imread(jpgFileName);
%         else
%             fprintf('File %s does not exist.\n', jpgFileName);
%         end
%     end
% end

%%
% imageDataBinary = zeros(5327,9185,157,'uint8');
% for m = 1:length(imageData(1,1,:))
%     for i = 1:length(imageData(:,1,1))
%          for j = 1:length(imageData(1,:,1))
%              if(imageData(i,j,m) > 5)
%                  imageDataBinary(i,j,m) = 1;
%              else
%                  imageDataBinary(i,j,m) = 0;
%              end
%          end
%     end
% end

% image1CC = bwconncomp(image1);
% image1CC_centroids = regionprops(image1CC,'Centroid');

% cellpoints = double(zeros(68344,3));
% cellcounter = 1;
% for m = 1:157
%     tempregionprops = regionprops(bwconncomp(imageDataBinary(:,:,m)));
%     for n = 1:length(tempregionprops)
%         temp = getfield(tempregionprops(n),'Centroid');
%         cellpoints(cellcounter,1) = temp(1);
%         cellpoints(cellcounter,2) = temp(2);
%         cellpoints(cellcounter,3) = m;
%         cellcounter = cellcounter + 1;
%     end
% end
% 
% plot3(cellpoints(:,1),cellpoints(:,2),cellpoints(:,3),'.'); 