function storedlocation = smallerimages(picsize)
    address = uigetdir('temp', 'Select folder where images are');
    a = dir([address '\*.tif']);
    numberofslices = numel(a);
    
    address2 = uigetdir('temp', 'Select folder where to store smaller images');
    storedlocation = address2;
    
    tic;    
    %Function to separate every tif image into smaller images, for only
    %those that contain cells in them. The data is then saved.
    for slice = 1:numberofslices
        img = 1;
        
        if(slice<10)
            img_cells = imread(address + "\Stitched_Z00" + num2str(slice)+".tif");
        elseif(slice<100)
            img_cells = imread(address + "\Stitched_Z0" + num2str(slice)+".tif");
        else
            img_cells = imread(address + "\Stitched_Z" + num2str(slice)+".tif");
        end
        
        for ro = 1:floor(size(img_cells,1)/picsize)
            for co = 1:floor(size(img_cells,2)/picsize)
                props = regionprops(img_cells(picsize*(ro-1)+1 : picsize*ro,picsize*(co-1)+1:picsize*co));
                if(max([props.Area]) > 50)
                    
                    if(img<10)
                        picnr = "000" + num2str(img);
                    elseif(img<100)
                        picnr = "00" + num2str(img);
                    elseif(img<1000)
                        picnr = "0" + num2str(img);
                    else
                        picnr = num2str(img);
                    end
                    
                    if(slice<10)
                        slicenr = "00" + num2str(slice);
                    elseif(slice<100)
                        slicenr = "0" + num2str(slice);
                    else
                        slicenr = num2str(slice);
                    end
                    
                    if(ro<10)
                        ronr = "00" + num2str(ro);
                    elseif(ro<100)
                        ronr = "0" + num2str(ro);
                    else
                        ronr = num2str(ro);
                    end
                    
                    if(co<10)
                        conr = "00" + num2str(co);
                    elseif(co<100)
                        conr = "0" + num2str(co);
                    else
                        conr = numstr(co);
                    end
                    
                    imwrite(img_cells(picsize*(ro-1)+1 : picsize*ro,picsize*(co-1)+1:picsize*co),...
                    address2 + "\slice"+slicenr+"pic"+picnr+"ro"+ronr+"co"+conr+".tif",'tif');
                    img = img + 1;
                    display("Image "+ num2str(img)+" Slice " + num2str(slice)+ " Saved");
                end
            end
        end
        display(num2str(slice/numberofslices*100)+"% complete");
    end
    
    display("Code took "+ num2str(toc/60) + " minutes");
end

