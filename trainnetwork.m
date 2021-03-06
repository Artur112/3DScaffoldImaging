%% Train Segmentation Network 
picsize = 300; %How big the images' width and heights are

dataDir = fullfile(toolboxdir('vision'),'visiondata');
imDir = 'C:\users\Artur\Desktop\temp\1_Smallerpics\'; %Directory where raw images are
pxDir = 'C:\users\Artur\Desktop\temp\3_Final\'; %Directory where pixel labeled images are

imds = imageDatastore(imDir); %Imagedatastore for all images
classNames = ["background", "cell"];
pxlLabels = [0, 1];
pxds = pixelLabelDatastore(pxDir,classNames,pxlLabels); %Pixel labeldatastore for all images
pximds = pixelLabelImageDatastore(imds,pxds);

splitratio = 0.75; %Training and testing data split ratio
train_pximds = partitionByIndex(pximds,1:floor(splitratio*length(imds.Files)));
test_pximds = partitionByIndex(pximds,ceil(splitratio*length(imds.Files)):length(imds.Files));

%Specifying the neural network architecture
numFilters = 64;
filterSize = 3;
numClasses = 2;
layers = [
    imageInputLayer([picsize picsize 1])
    convolution2dLayer(filterSize,numFilters,'Padding',1)
    reluLayer()
    maxPooling2dLayer(2,'Stride',2)
    convolution2dLayer(filterSize,numFilters,'Padding',1)
    reluLayer()
    transposedConv2dLayer(4,numFilters,'Stride',2,'Cropping',1);
    convolution2dLayer(1,numClasses);
    softmaxLayer()
    pixelClassificationLayer()
    ];

opts = trainingOptions('sgdm', ...
    'InitialLearnRate',1e-3, ...
    'MaxEpochs',100, ...
    'MiniBatchSize',15, ...
    'Plots', 'training-progress');

net = trainNetwork(train_pximds,layers,opts);

save('../networktemp.mat','net');

%%

imgs_test = imageDatastore(test_pximds.Images);
classNames = ["background", "cell"];
pxlLabels = [0, 1];
pxdsTruth = pixelLabelDatastore(test_pximds.Images,classNames,pxlLabels);

for n = 1:length(imgs_test.Files)
    imagetotest = partition(imgs_test,'Files',n);
    imagetruth = partition(pxdsTruth,'Files',n);
    name = cell2mat(imagetotest.Files);
    pxdsResults = semanticseg(imagetotest,net,"WriteLocation",'D:\Artur\NetworkOutput','NamePrefix', name(end-28:end-4));
   % metrics = evaluateSemanticSegmentation(pxdsResults,imagetruth);
    display(n);
end
