function samplescans(picsize,numberofslices, address)
    tic;
    %Function to separate every tif image into smaller images, for only
    %those that contain cells in them. The data is then saved.
    for slice = 1:numberofslices
        if(slice<10)
            img_cells = imread(address + "\Stitched_Z00" + num2str(slice)+".tif");
        elseif(slice<100)
            img_cells = imread(address + "\Stitched_Z0" + num2str(slice)+".tif");
        else
            img_cells = imread(address + "\Stitched_Z" + num2str(slice)+".tif");
        end
        
        img = 1;
        for ro = 1:floor(size(img_cells,1)/picsize)
            for co = 1:floor(size(img_cells,2)/picsize)
                props = regionprops(img_cells(picsize*(ro-1)+1 : picsize*ro,picsize*(co-1)+1:picsize*co));
                if(max([props.Area]) > 50)
                    imwrite(img_cells(picsize*(ro-1)+1 : picsize*ro,picsize*(co-1)+1:picsize*co),...
                    "D:\ARTUR\Week2SmallerPics\slice"+num2str(slice)+"_img"+num2str(img)+".tif",'tif');
                    img = img + 1;
                    display("Image "+ num2str(img)+" Slice " + num2str(slice)+ " Saved");
                end
            end
        end
        display(num2str(slice/numberofslices*100)+"% complete");
    end
    
    display("Code took "+ num2str(toc/60) + " minutes");
end

