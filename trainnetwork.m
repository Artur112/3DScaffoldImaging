%% Train Segmentation Network 
dataDir = fullfile(toolboxdir('vision'),'visiondata');
imDir = 'D:\ARTUR\Week3SmallerPics-Nuclei\'; %Directory where raw images are
pxDir = 'D:\ARTUR\Week3Final-Nuclei\'; %Directory where pixel labeled images are

imds = imageDatastore(imDir); %Imagedatastore for all images

classNames = ["background", "cell"];
pxlLabels = [0, 1];
pxds = pixelLabelDatastore(pxDir,classNames,pxlLabels); %Pixel labeldatastore for all images

pximds = pixelLabelImageDatastore(imds,pxds);

splitratio = 0.75;

train_pximds = partitionByIndex(pximds,1:floor(splitratio*length(imds.Files)));
test_pximds = partitionByIndex(pximds,ceil(splitratio*length(imds.Files)):length(imds.Files));

numFilters = 64;
filterSize = 3;
numClasses = 2;
layers = [
    imageInputLayer([300 300 1])
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

save('../network.mat','net');

%%

imgs_test = imageDatastore(test_pximds.Images);
classNames = ["background", "cell"];
pxlLabels = [0, 1];
pxdsTruth = pixelLabelDatastore(test_pximds.Images,classNames,pxlLabels);

for n = 1:length(imgs_test.Files)
    imagetotest = partition(imgs_test,'Files',n);
    imagetruth = partition(pxdsTruth,'Files',n);
%net = load('net.mat');
%net = SeriesNetwork(net.layers);
    name = cell2mat(imagetotest.Files);
    pxdsResults = semanticseg(imagetotest,net,"WriteLocation",'D:\Artur\NetworkOutput','NamePrefix', name(end-28:end-4));
   % metrics = evaluateSemanticSegmentation(pxdsResults,imagetruth);
    display(n);
end
