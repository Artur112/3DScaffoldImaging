run("Image Sequence...", "open=D:/ARTUR/RawData/Fib_Week3/Ch1_Stitched_Sections/Stitched_Z001.tif scale=10 convert sort");
selectWindow("Ch1_Stitched_Sections");

//run("Brightness/Contrast...");
run("Enhance Contrast", "saturated=0.35");
run("Apply LUT", "stack");
run("Close");
run("Image Sequence...", "open=C:/Users/artur/Desktop/Project_data_new/Cells/Week3/Slice_001.png");
selectWindow("Week3");
run("Histogram", "stack");
close();
run("Close");
run("Merge Channels...", "c2=Week3 c4=Ch1_Stitched_Sections create");
saveAs("Tiff", "C:/Users/artur/Desktop/Project_data_new/3D/Week3_Composite.tif");
//setTool("hand");
run("3D Viewer");
call("ij3d.ImageJ3DViewer.setCoordinateSystem", "false");
call("ij3d.ImageJ3DViewer.add", "Week3_Composite.tif", "None", "Week3_Composite.tif", "0", "true", "true", "true", "2", "0");
