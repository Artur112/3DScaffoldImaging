clc;clear all;close all;
%%filePattern = fullfile('Raw Data\Ch2_Stitched_Sections_JPEG\Ch2_Stitched_Sections_JPEG\','*.jpg');

image1 = imread('Raw Data\Ch2_Stitched_Sections_JPEG\Ch2_Stitched_Sections_JPEG\Stitched_Z001.jpg');
for i = 1:length(image1(:,1))
    for j=1:length(image1(1,:))
        if(image1(i,j) > 0)
            image1binary(i,j) = 1;
        else
            image1binary(i,j) = 0;
        end
    end
end

image1CC = bwconncomp(image1);
image1CC_centroids = regionprops(image1CC,'Centroid');

hold on;
for n = 1:length(image1CC_centroids)
    temp = getfield(image1CC_centroids(n),'Centroid');
    plot(temp(1),temp(2),'*');    
end

% for k = 1:50
% 	% Create an image filename, and read it in to a variable called imageData.
% 	jpgFileName = strcat('Raw Data\Ch2_Stitched_Sections_JPEG\Ch2_Stitched_Sections_JPEG\Stitched_Z00', num2str(k), '.jpg');
% 	if exist(jpgFileName, 'file')
% 		imageData(:,:,k) = imread(jpgFileName);
% 	else
% 		fprintf('File %s does not exist.\n', jpgFileName);
%     end
% end
% 
% for m = 1:length(imageData(1,1,:))
%     for i = 1:length(imageData(:,1,1))
%          for j = 1:length(imageData(:,1,:))
%              if(imageData(i,j,m) > 0)
%                  imageDataBinary(i,j,m) = 1;
%              else
%                  imageDataBinary(i,j,m) = 0;
%              end
%          end
%     end
% end

