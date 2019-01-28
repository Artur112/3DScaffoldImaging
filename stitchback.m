function stitchback(numrows, numcols, smallerimgsize)
    address = uigetdir('temp', 'Select folder where smaller images are');
    a = dir([address '\*.png']);
    
    address2 = uigetdir('temp', 'Select folder where to store bigger images');
        
    slice = 1;
    image = zeros(numrows,numcols);
    disp("Stitching slice 1");
    
    for file = a'
        name = file.name;
        slicenr = str2double(name(6:8));
        rownr = str2double(name(18:20));    
        colnr = str2double(name(23:25));
        
        if(slicenr == slice) 
            image(smallerimgsize*(rownr-1)+1:rownr*smallerimgsize, smallerimgsize*(colnr-1)+1:colnr*smallerimgsize) = imread(address+"\"+name);
        else
           imwrite(image,address2 + "\Slice" + num2str(slice) + ".tif");
           slice = slicenr;
           disp("Stitching slice " + num2str(slice));
           image = zeros(numrows,numcols);
           image(smallerimgsize*(rownr-1)+1:rownr*smallerimgsize, smallerimgsize*(colnr-1)+1:colnr*smallerimgsize) = imread(address+"\"+name);
        end    
    end
end

