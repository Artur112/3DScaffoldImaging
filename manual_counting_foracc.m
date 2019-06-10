%% Selecting random images to manually count cells from
%Run until a good suitable image is found, then quit program and write down
%details of the image section
while true
    week = randsample(3,1);
    files = dir([strcat('C:\users\Artur\Desktop\Project_data_new\Cells\Week',num2str(week)) '\*.png']);
    img_nr = randsample(10,1);
    img = imread(strcat(files(1).folder,'\',files(img_nr).name));
    y_pos = randsample(size(img,1)-5000,1);
    x_pos = randsample(size(img,2)-5000,1);
    img = img(y_pos:y_pos+1999,x_pos:x_pos+1999);
    if(~isempty(img))
    imshow(img>0);
    % files = dir([strcat('D:\Artur\RawData\Fib_week',num2str(week),'\Ch2_Stitched_Sections') '\*.tif']);
    pause;
    end
end

%% After image sections have been seleted take those sections pad with 0s and store
%Load in excel file first
for ii = 1:10
    disp(ii);
    weeknr = counting(ii,1);
    imgnr = counting(ii,2);
    xpos = counting(ii,3);
    ypos = counting(ii,4);
    if(imgnr<10)
        name = strcat('Stitched_Z00',num2str(imgnr),'.tif');
    else
        name = strcat('Stitched_Z0',num2str(imgnr),'.tif');
    end
    img = imread(strcat('D:\Artur\RawData\Fib_week',num2str(weeknr),'\Ch2_Stitched_Sections\',name));
    img = img(ypos:ypos+1999,xpos:xpos+1999);
    pad_img = uint16(zeros(2800,2800));
    pad_img(401:2400,401:2400) = img;
    imwrite(pad_img, "C:\Users\artur\Desktop\4thYear\Project\Images\To_Count\tocount" + num2str(ii) + ".tif");
end
    
%% After manual counting run Neural Network and compare
Num_Cells = [];
idx = 1;
net = load('network1.mat'); %Load neural Network
net = net.net;
picsize = 300;  
files = dir(['C:\users\Artur\Desktop\4thYear\Project\Images\To_Count' '\*.tif']);
for file = files'
    disp(idx);
    slicedata = [];
    %%% Segment cell image by sliding a 300x300 section along the image that is segmented with neural network
        cell_img = imread(strcat(file.folder,'\',file.name));
        cell_img_segm = zeros(size(cell_img,1),size(cell_img,2));
        for ro = 1:floor(size(cell_img,1)/picsize)
            for co = 1:floor(size(cell_img,2)/picsize)
                section = cell_img(1+picsize*(ro-1) : picsize*ro , 1+picsize*(co-1) : picsize*co); %300x300 section to analyze
                if(any(section,'all')) %Skip if contains no cells          
                    %Find if touching any cells along edge,[left edge, right edge, top edge, bottom edge]
                    touching = [any(section(:,1)),any(section(:,picsize)),any(section(1,:)),any(section(picsize,:))];
                    C = semanticseg(section,net);
                    C = double(C)-1;
                    if(any(touching)) %If touching fill in the gaps in the edges
                        for m = 1:4
                            if(touching(m))
                                switch m
                                    case 1
                                        C(section(:,1)>0,1) = 1;
                                    case 2
                                        C(section(:,picsize) > 0,picsize) = 1;
                                    case 3
                                        C(1,section(1,:) > 0) = 1;
                                    case 4
                                        C(picsize,section(picsize,:)>0) = 1;
                                end
                            end
                        end                      
                    end
                    cell_img_segm(1+picsize*(ro-1) : picsize*ro , 1+picsize*(co-1) : picsize*co) = C;  
                end
            end
        end
        
    %%% Calculate Metrics  
        cc = bwconncomp(cell_img_segm);
        slicedata.Centroids = cell2mat(struct2cell(regionprops(cc,'Centroid'))'); %Centroids of all cells
        slicedata.Areas = cell2mat(struct2cell(regionprops(cc,'Area'))'); %Areas of all cells
        
        %Find areas that are smaller than 50 (these are cell fragments that count as noise or cells that partially show up from deeper layers.
        toremove = find(slicedata.Areas < 50);
        
        %Remove these areas from the cell image
        pxls = regionprops(cc,'PixelList');
        for ii = 1:length(toremove)
            cell_img_segm(pxls(toremove(ii)).PixelList(:,2), pxls(toremove(ii)).PixelList(:,1)) = 0;
        end
        
        % Remove these metrics from the metric_data
        slicedata.Centroids(toremove,:) = [];
        slicedata.Areas(toremove,:) = [];
        
    %%% Store numceslls for image
        Num_Cells(idx) = length(slicedata.Centroids);
        idx = idx + 1;
end
    %%
counted = [40,46,49,63,28,41,58,45,58,50];
diff  = abs(counted - Num_Cells);
perc = diff./counted*100;
mean(perc);
rmse = sqrt(mean(counted-Num_Cells)^2);
    
    