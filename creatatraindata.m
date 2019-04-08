% clear all; clc;
% 
% disp('Please select where the data to segment is');
% data_addr = uigetdir('temp', 'Select folder where data is');
% a = dir([data_addr '\*.tif']);
% numberofslices = numel(a);
% 
% picsize = 300; %what size images to be split into
% store_addr = "C:\users\Artur\Desktop\temp"; %where to store the images
% 
% if ~exist(store_addr, 'dir')
%     mkdir(store_addr);
%     mkdir(store_addr + "\1_Smallerpics");
%     mkdir(store_addr + "\2_Segmented");
%     mkdir(store_addr + "\3_Final");
% end
% 
% %% Split images into smaller images (for Neural Network training)
% disp('Splitting data into smaller images');
% tic;    
% for slice = 1:1%numberofslices
%     img = 1;
%     if(slice<10)
%         img_cells = imread(data_addr + "\Stitched_Z00" + num2str(slice)+".tif");
%     elseif(slice<100)
%         img_cells = imread(data_addr + "\Stitched_Z0" + num2str(slice)+".tif");
%     else
%         img_cells = imread(data_addr + "\Stitched_Z" + num2str(slice)+".tif");
%     end
% 
%     for ro = 1:floor(size(img_cells,1)/picsize)
%         for co = 1:floor(size(img_cells,2)/picsize)
%             props = regionprops(img_cells(picsize*(ro-1)+1 : picsize*ro,picsize*(co-1)+1:picsize*co));
%             if(max([props.Area]) > 50)
% 
%                 if(img<10)
%                     picnr = "000" + num2str(img);
%                 elseif(img<100)
%                     picnr = "00" + num2str(img);
%                 elseif(img<1000)
%                     picnr = "0" + num2str(img);
%                 else
%                     picnr = num2str(img);
%                 end
% 
%                 if(slice<10)
%                     slicenr = "00" + num2str(slice);
%                 elseif(slice<100)
%                     slicenr = "0" + num2str(slice);
%                 else
%                     slicenr = num2str(slice);
%                 end
% 
%                 if(ro<10)
%                     ronr = "00" + num2str(ro);
%                 elseif(ro<100)
%                     ronr = "0" + num2str(ro);
%                 else
%                     ronr = num2str(ro);
%                 end
% 
%                 if(co<10)
%                     conr = "00" + num2str(co);
%                 elseif(co<100)
%                     conr = "0" + num2str(co);
%                 else
%                     conr = numstr(co);
%                 end
% 
%                 imwrite(img_cells(picsize*(ro-1)+1 : picsize*ro,picsize*(co-1)+1:picsize*co),...
%                 store_addr + "\1_Smallerpics\slice"+slicenr+"pic"+picnr+"ro"+ronr+"co"+conr+".tif",'tif');
%                 display("Image "+ num2str(img)+" Slice " + num2str(slice)+ " Saved");
%                 img = img + 1;
%             end
%         end
%     end
%     display(num2str(slice) + " / " + num2str(numberofslices) + " slices done");
% end

%% Running ImageJ MorphoLibJ plugin to segment the smaller images
disp('Segmenting the smaller images with MorpholibJ');
a = dir([convertStringsToChars(store_addr + "\1_SmallerPics") '\*.tif']);

file = fopen('C:\Users\artur\Desktop\4thYear\Project\Project-Code\filename.txt','w');
fprintf(file,strrep(store_addr,'\','/') + "/2_Segmented/" +"\n");
for n = a'
    fprintf(file, strrep(store_addr,'\','/') + "/1_Smallerpics/" + n.name + "\n");
    fprintf(file, n.name(1:end-4)+"_segm" + "\n");
end
fclose(file);

command = 'vdesk on:2 run:"C:\Users\artur\Desktop\4thYear\Project\Fiji.app\ImageJ-win64.exe" -macro C:/Users/Artur/Desktop/4thYear/Project/Fiji.app/plugins/MorphSegmentation.ijm';
system(command);

%% Make segmented images suitable for neural network
disp('Correcting images to be suitable for Neural Network');
a = dir([convertStringsToChars(store_addr + "\2_Segmented") '\*.tif']);
total = numel(a);

%iter = 1;
for file = a'

    X = imread(store_addr + "\2_Segmented\" + file.name);  

    %Swap ones (background) and zeros(edges) in output from imagej to
    %ones for edges and zeros for background
    X2 = X;
    for m = 1:size(X,1)
        for n = 1:size(X,2)
            if(X(m,n) == 0)
                X2(m,n) = 1;
            elseif(X(m,n) == 1)
                X2(m,n) = 0;
            end
        end
    end

    %Remove edges from image
    X3 = X2;
    for m = 1:size(X2,1)
        for n = 1:size(X2,2)
            if(X2(m,n) == 1)
                X3(m,n) = 0;
            end
        end
    end

    %Find areas and positions of the cells / fragments
    %with edges removed
    ka = regionprops(X3);
    Area = [ka.Area];
    Centroids = {ka.Centroid};
    BoundingBoxes = {ka.BoundingBox};
    X4 = X2;

    %Loop through the areas found
    for blob = 1:length(Area)
        %If its area is smaller than 50 or bigger than 5000, discard it, by looping through all the cells
        %inside the bounding box and deleting those that are equal to
        %the class label (val ue of cell at centroid). Looping through
        %the bounding box instead of the whole image was done for
        %speed.
        if(Area(blob) < 50 && Area(blob) ~= 0 || Area(blob) > 5000 && Area(blob) ~= 0)
           labl = X2(floor(Centroids{blob}(2)), floor(Centroids{blob}(1)));
           for m = floor(Centroids{blob}(2)-BoundingBoxes{blob}(4)/2): round(Centroids{blob}(2) + BoundingBoxes{blob}(4)/2)
               for n = floor(Centroids{blob}(1)-BoundingBoxes{blob}(3)/2): round(Centroids{blob}(1) + BoundingBoxes{blob}(3)/2)
                   %Dont check outside image
                   if(m<301 && n<301 && m>0 && n>0)
                       if(X2(m,n) == labl)
                           X4(m,n) = 0;
                       end
                   end
               end
           end
        end
    end

    %Deleting the remaining edges of the globs that were removed. Done
    %by deleting every edge that isnt next to a class label.
    X5 = X4;
    for m = 1:size(X5,1)
        for n = 1:size(X5,2)
            if(X5(m,n) == 1)
                nums = 0;
                for a = -1:1
                    for b = -1:1
                        %dont check outside edges
                        if(m+a < 301 && n+b < 301 && m+a>0 && n+b > 0)
                            if(X5(m+a,n+b) > 1)
                                nums = nums + 1;
                            end
                        end
                    end
                end
                if(nums == 0)
                    X5(m,n) = 0;
                end
            end
        end
    end      

    % Remove class labels - make them all equal to 1
    for m = 1:size(X5,1)
        for n = 1:size(X5,2)
            if(X5(m,n)>1)
                X5(m,n) = 1;
            end
        end
    end
    name = file.name;
    imwrite(uint8(X5), store_addr + "\3_Final\" + name(1:end-8) + ".png");
    %display(num2str(iter/total*100)+"% complete");
    %iter = iter + 1;
end