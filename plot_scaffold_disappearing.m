a = 1;
ha = tight_subplot(3,5,[.01 .03],[.1 .01],[.01 .01]);       
for week = 1:3
    files = dir([strcat('D:\Artur\RawData\fib_week',num2str(week),'\Ch1_Stitched_Sections\') '*.tif']);
    for slice = [1,2,10,20,40]
        disp(a);
        img = imread(strcat(files(slice).folder,'\',files(slice).name));
        img = imresize(img,0.1);
        %subplot(3,5,a);
        axes(ha(a));
        imagesc(img);
        colormap('gray');
        set(gca,'visible','off');
        if(slice == 1)
            ylabel("Week " + num2str(week));
        end
        title("Slice " + num2str(slice));
        a = a+1; 
    end
end