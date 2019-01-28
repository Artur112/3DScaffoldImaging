function groundtruth()
    address = uigetdir('temp', 'Select folder where MorphoLibJ images are');
    a = dir([address '\*.tif']);
    total = numel(a);
    
    address2 = uigetdir('temp', 'Select folder where to store ground truth images');
    
    iter = 1;
    for file = a'

            X = imread(address + "\" + file.name);  
           
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
            %by deleting every edges that isnt next to a class label.
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

            % Remove class labels - make them equal to 1
            for m = 1:size(X5,1)
                for n = 1:size(X5,2)
                    if(X5(m,n)>1)
                        X5(m,n) = 1;
                    end
                end
            end
            name = file.name;
            imwrite(uint8(X5), address2 + "\" + name(1:end-8) + ".png");
            display(num2str(iter/total*100)+"% complete");
            iter = iter + 1;
    end
end

