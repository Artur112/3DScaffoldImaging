function stitchback(numrows, numcols, smallerimgsize)
    address = uigetdir('temp', 'Select folder where smaller images are');
    %address = 'D:/Artur/kaka';
    a = dir([address '\*.png']);
    
    address2 = uigetdir('temp', 'Select folder where to store bigger images');
    %address2 = 'D:/Artur/kaka2';
    total = numel(a);
    slice = 1;
    lastslicenr = str2double(a(total).name(6:8));
    image = zeros(numrows,numcols);
    disp("Stitching slice 1 of " + num2str(lastslicenr));
    iter = 1;
    for file = a'
        name = file.name;
        slicenr = str2double(name(6:8));
        rownr = str2double(name(18:20));    
        colnr = str2double(name(23:25));
        if(iter~= total)
            if(slicenr == slice) 
                image(smallerimgsize*(rownr-1)+1:rownr*smallerimgsize, smallerimgsize*(colnr-1)+1:colnr*smallerimgsize) = imread(address+"\"+name);
            else
               imwrite(image,address2 + "\Slice" + num2str(slice) + ".tif");
               slice = slicenr;
               disp("Stitching slice " + num2str(slice + " of " + num2str(lastslicenr)));
               image = zeros(numrows,numcols);
               image(smallerimgsize*(rownr-1)+1:rownr*smallerimgsize, smallerimgsize*(colnr-1)+1:colnr*smallerimgsize) = imread(address+"\"+name);
            end    
        else
            imwrite(image,address2 + "\Slice" + num2str(slice) + ".tif");
            disp('Stitching Complete');
        end
        iter = iter + 1;
    end
end

