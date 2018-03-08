%% Load data:
addpath(genpath('D:\GitHub\other\NMF-Source-Extraction'));
load('D:\Data\scratch\MM102_NMF_results\MM102_160727_main_slice01_patchResults_v0913.mat');

%% Selmaans clustering:
nClusters = 4;
sourceProps = clusterSourceTypes(A, nClusters);

% Plot high-probability ROIs:
clusterProb = max(sourceProps.gmProbs, [], 2);
thresh = 0.90;

cIm = zeros(512,512,3);
cIm(:,:,1) = reshape(sum(A(:,clusterProb>thresh & sourceProps.idx==sourceProps.sigRank(1)),2),512,512);
cIm(:,:,2) = reshape(sum(A(:,clusterProb>thresh & sourceProps.idx==sourceProps.sigRank(2)),2),512,512);
cIm(:,:,3) = reshape(sum(A(:,clusterProb>thresh & sourceProps.idx==sourceProps.sigRank(3)),2),512,512);
cIm(:,:,3) = cIm(:,:,3) + reshape(sum(A(:,clusterProb>thresh & sourceProps.idx==sourceProps.sigRank(4)),2),512,512);
cIm(:,:,2) = cIm(:,:,2) + reshape(sum(A(:,clusterProb>thresh & sourceProps.idx==sourceProps.sigRank(4)),2),512,512);

figure(4)
imshow(cIm*4, 'parent', gca)

sourceClass = sourceProps.idx;
sourceClass(clusterProb<thresh) = nan;
% save('MM102_160725_goodTrainingLabels.mat', 'sourceClass');

%% Create source patches:
% Align Sources
winRad = 12;
winWidth = 2*winRad+1;
nROIs = size(A,2);
alignedMasks = zeros(winWidth,winWidth,nROIs);
for nROI = 1:nROIs
    thisMask = padarray(reshape(A(:,nROI),512,512),[winRad,winRad],0);
    maskProps = regionprops(thisMask>0, thisMask, 'WeightedCentroid','Area');
    if length(maskProps)>1
        [~,thisComp] = max([maskProps.Area]);
        maskProps = maskProps(thisComp);
    end
    thisCent = round(maskProps.WeightedCentroid([2,1]));
    tempMask = thisMask(thisCent(1)-winRad:thisCent(1)+winRad,...
        thisCent(2)-winRad:thisCent(2)+winRad);
    tempMask = medfilt2(tempMask,[3 3]);
    alignedMasks(:,:,nROI) = tempMask;
end

% save('MM102_160725_rawTrainingData.mat', 'alignedMasks');

%% Cluster with current NN:
Xtest = permute(alignedMasks, [1 2 4 3]);
YTest = classify(convnet, Xtest);
% YTest = addcats(YTest, 'cutEdge');
% YTest = renamecats(YTest, {'doughnutSoma', 'outOfPlaneSoma', ...
%                 'smallRoundProcess', 'complexProcess', ...
%                 'messyAndMultiples', 'cutEdge'});
classes = categories(YTest);
classHere = 1;

cIm = mat2gray(reshape(full(sum(A(:,YTest==classes(classHere)),2)),512,512));
figure(44563456)
imshow(cIm, 'initialmag', 'fit')
title(classes{classHere})

% Refine manually:
% o = classifierGui(A, double(YTest), sourceProps.allCentroids);

%% Load training data from file:
fileNames = {'MM102_160725_manualSourceLabels', ...
             'MM102_160726_manualSourceLabels', ...
             'MM102_160727_manualSourceLabels'};
         
for i = 1:numel(fileNames)
    if i == 1
        manLabels = load('MM102_160725_manualSourceLabels');
    else
        manLabelsHere = load('MM102_160725_manualSourceLabels');
        assert(isequal(manLabels.classNames, manLabelsHere.classNames), 'saved class names don''t match -- check!');
        manLabels.labels = cat(1, manLabels.labels, manLabelsHere.labels);
        manLabels.sources = cat(2, manLabels.sources, manLabelsHere.sources);
    end
end

% Align Sources
winRad = 12;
winWidth = 2*winRad+1;
nROIs = size(manLabels.sources,2);
alignedMasks = zeros(winWidth,winWidth,nROIs);
for nROI = 1:nROIs
    thisMask = padarray(reshape(manLabels.sources(:,nROI),512,512),[winRad,winRad],0);
    maskProps = regionprops(thisMask>0, thisMask, 'WeightedCentroid','Area');
    if length(maskProps)>1
        [~,thisComp] = max([maskProps.Area]);
        maskProps = maskProps(thisComp);
    end
    thisCent = round(maskProps.WeightedCentroid([2,1]));
    tempMask = thisMask(thisCent(1)-winRad:thisCent(1)+winRad,...
        thisCent(2)-winRad:thisCent(2)+winRad);
    tempMask = medfilt2(tempMask,[3 3]);
    alignedMasks(:,:,nROI) = tempMask;
end
imgSize = [winWidth winWidth];

%% Augment data:
realLabels = manLabels.labels;
augmentedMasks = alignedMasks;
augmentedClass = realLabels;
nClasses = numel(unique(augmentedClass));
classSize = histcounts(augmentedClass);
classSizeReal = classSize;
nRealExamples = sum(~isnan(sourceClass));

realInd = {};
for i = 1:nClasses
    realInd{i} = find(augmentedClass==i);
end

% Pre-allocate:
nTotal = 200000;
augmentedMasks(:,:,nTotal) = 0;
augmentedClass(nTotal) = 0;

tformTemplate = projective2d(eye(3));
maxTranslPix = 2;
outRef = imref2d(size(augmentedMasks(:,:,1)));
outRef.XWorldLimits = outRef.XWorldLimits-imgSize(1)/2-0.5;
outRef.YWorldLimits = outRef.YWorldLimits-imgSize(2)/2-0.5;

for i = nRealExamples + 1:nTotal
    i
    
%     % Pick real example from smallest class at random:
    [~, classHere] = min(classSize);
    classSize(classHere) = classSize(classHere)+1;
    realIndHere = realInd{classHere}(randi(classSizeReal(classHere)));

    % Pick real example (maintaining class size bias):
%     realIndHere = randi(nRealExamples);
    
    % Create random rotation and translation:
    augmentedClass(i) = realLabels(realIndHere);
    x = (rand*2-1) * maxTranslPix;
    y = (rand*2-1) * maxTranslPix;
    
    if augmentedClass(i)==6 
        % Don't rotate cut-edge sources!
        th = 0;
    else
        th = rand * 360;
    end
    
    % Matrix representing a rotation followed by a translation:
    % http://planning.cs.uiuc.edu/node99.html
    tform = tformTemplate;
    tform.T = [cosd(th) -sind(th) x;
        sind(th)  cosd(th) y;
        0         0 1]';
    
    augmentedMasks(:,:,i) = imwarp(alignedMasks(:,:,realIndHere), ...
        outRef, tform, 'FillValues', 0, ...
        'OutputView', outRef);
end
% save('MM102_160725-27_augmentedTrainingData_2pxTransl.mat', 'augmentedMasks', 'augmentedClass', '-v7.3');
% load('MM102_160725-27_augmentedTrainingData_2pxTransl.mat')

%% Adjust data format:

% Input for NN, with dimensions h-w-ch-n:
X = permute(augmentedMasks, [1 2 4 3]);
Y = categorical(augmentedClass);
nCat = numel(unique(Y));
classSizes = histcounts(Y);

% % Equalize class sizes:
% nMin = min(classSizes);
% Ydouble = double(Y);
% for i = 1:nCat
%     iThis = find(Ydouble==i);
%     if numel(iThis) > nMin
%         X(:,:,:, iThis(nMin+1:end)) = [];
%         Y(iThis(nMin+1:end)) = [];
%         Ydouble(iThis(nMin+1:end)) = [];
%     end
% end

% Randomize order:
ord = randperm(numel(Y));
Y = Y(ord);
X = X(:,:,:,ord);

fprintf('Class sizes: ')
fprintf('%d ', histcounts(Y));
fprintf('\n')

%% Display some:
% figure
cats = categories(Y);
for i = 2000:10000
    if Y(i)==cats(6)
        imagesc(X(:,:,:,i));
        drawnow
        pause(0.2);
    else
        continue
    end
end
    

%% Neural network

% Architecture:
% layers = [imageInputLayer(imgSize,'Normalization','zerocenter','Name','inputl')
%     convolution2dLayer([4 3],12,'NumChannels',1,'Name','conv1')
%     reluLayer('Name','relu1')
%     crossChannelNormalizationLayer(4,'Name','cross1')
%     maxPooling2dLayer(2,'Stride',2,'Name','max1')
%     fullyConnectedLayer(nCat,'Name','full2')
%     softmaxLayer('Name','softm')
%     classificationLayer('Name','out')];

% layers = [imageInputLayer(imgSize,'Normalization','zerocenter','Name','inputl')
%     convolution2dLayer(15, 32, 'Stride', 1, 'Name', 'conv1')
%     reluLayer('Name','relu1')
%     crossChannelNormalizationLayer(5,'Name','cross1')
%     convolution2dLayer(10, 4, 'Stride', 1, 'Name', 'conv2')
%     reluLayer('Name','relu1')
%     crossChannelNormalizationLayer(5,'Name','cross1')
%     reluLayer('Name','relu1')
%     fullyConnectedLayer(256, 'Name', 'full1')
%     reluLayer('Name', 'relu4')
%     fullyConnectedLayer(nCat, 'Name', 'full2')
%     softmaxLayer('Name', 'softm')
%     classificationLayer('Name', 'out')];

layers = [imageInputLayer(imgSize,'Normalization','zerocenter','Name','inputl')
    convolution2dLayer(18, 32, 'Stride', 1, 'Name', 'conv1', 'Padding', 0)
    reluLayer('Name','relu1')
    crossChannelNormalizationLayer(5,'Name','cross1')
    convolution2dLayer(5, 3, 'Stride', 1, 'Name', 'conv1', 'Padding', 0)
    reluLayer('Name','relu1')
    crossChannelNormalizationLayer(5,'Name','cross1')
    fullyConnectedLayer(256, 'Name', 'full1')
    reluLayer('Name', 'relu4')
    fullyConnectedLayer(nCat, 'Name', 'full2')
    softmaxLayer('Name', 'softm')
    classificationLayer('Name', 'out')];

% layers = [imageInputLayer(imgSize,'Normalization','zerocenter','Name','inputl')
%     fullyConnectedLayer(256, 'Name', 'full1')
%     reluLayer('Name', 'relu4')
%     fullyConnectedLayer(256, 'Name', 'full1')
%     reluLayer('Name', 'relu4')
%     fullyConnectedLayer(nCat, 'Name', 'full2')
%     softmaxLayer('Name', 'softm')
%     classificationLayer('Name', 'out')];

% Initialize with sqrt(n):
% for i = 1:numel(layers)
%     if isprop(layers(i), 'Weights')
%         if i==2
%             fanIn = prod(imgSize);
%         else
%             % Fan-in is previous number of weights:
%             fanIn = prod(weightDims);
%         end
%         switch class(layers(i))
%             case 'nnet.cnn.layer.Convolution2DLayer'
%                 weightDims = [layers(i).FilterSize, 1, layers(i).NumFilters];
%                 
%             case 'nnet.cnn.layer.FullyConnectedLayer'
%                 weightDims = 
%         layers(i).Weights = gpuArray(randn(weightDims) / sqrt(fanIn));
%     end
% end
            

%% Train:
% options = trainingOptions('sgdm','MaxEpochs',10000,...
%     'MiniBatchSize', 512, ...
%     'InitialLearnRate',0.0001, ...
%     'CheckpointPath', 'D:\GitHub\other\NMF-Source-Extraction\MM_util\clusterSources');
options = trainingOptions('sgdm','MaxEpochs',20,...
    'MiniBatchSize', 512, ...
    'L2Regularization', 0.0001, ...
    'InitialLearnRate', 0.1, ...
    'LearnRateSchedule', 'piecewise', ...
    'LearnRateDropFactor', 1, ...
    'LearnRateDropPeriod', 1);

nSamples = numel(Y);
% nSamples = 100;
convnet = trainNetwork(X(:,:,:,1:nSamples), Y(1:nSamples), layers, options);


% Refinement:
options = trainingOptions('sgdm','MaxEpochs',30,...
    'MiniBatchSize', 512, ...
    'L2Regularization', 0.0001, ...
    'InitialLearnRate', 0.03, ...
    'LearnRateSchedule', 'piecewise', ...
    'LearnRateDropFactor', 0.2, ...
    'LearnRateDropPeriod', 5, ...
    'CheckpointPath', 'D:\GitHub\other\NMF-Source-Extraction\MM_util\clusterSources');

nSamples = numel(Y);
convnet = trainNetwork(X(:,:,:,1:nSamples), Y(1:nSamples), convnet.Layers, options);

%% Check weights:
w = gather(convnet.Layers(2).Weights);
w = permute(w, [1 2 4 3]);
w = reshape(w, size(w, 1), size(w, 2)*4, []);
w = permute(w, [2, 1, 3]);
w = reshape(w, size(w, 1), []);
imagesc(w)
axis equal

%% Check performance on all examples:
% checkPoint = load('convnet_checkpoint__127140__2016_10_03__11_37_49.mat');
% net = checkPoint.net;
net = convnet;

Xtest = permute(alignedMasks, [1 2 4 3]);
[YTest, score] = classify(net, Xtest);
% YTest = categorical(manLabels.labels);
classes = categories(YTest);
classHere = 1;

classScore = max(score, [], 2);
cIm = zeros(512,512,3);
sourcesHere = full(manLabels.sources(:,YTest==classes(classHere)));
scoresHere = classScore(YTest==classes(classHere))';
cIm(:,:,1) = mat2gray(reshape(sum(bsxfun(@times, scoresHere*0+1, sourcesHere) , 2), 512, 512));
cIm(:,:,2) = mat2gray(reshape(sum(bsxfun(@times, scoresHere, sourcesHere) , 2), 512, 512));
cIm(:,:,3) = mat2gray(reshape(sum(bsxfun(@times, scoresHere, sourcesHere) , 2), 512, 512));
% cIm(:,:,2) = reshape(sum(A(:,clusterProb>thresh & sourceProps.idx==sourceProps.sigRank(2)),2),512,512);
% cIm(:,:,3) = reshape(sum(A(:,clusterProb>thresh & sourceProps.idx==sourceProps.sigRank(3)),2),512,512);
% cIm(:,:,3) = cIm(:,:,3) + reshape(sum(A(:,clusterProb>thresh & sourceProps.idx==sourceProps.sigRank(4)),2),512,512);
% cIm(:,:,2) = cIm(:,:,2) + reshape(sum(A(:,clusterProb>thresh & sourceProps.idx==sourceProps.sigRank(4)),2),512,512);

figure(44563456)
imshow(cIm, 'initialmag', 'fit')

% Refine manually:
% o = classifierGui(A, double(YTest), centroids);

%% Confusion matrix:
nY = numel(Y(1:nSamples));
targets = false(nCat, nY);
targets(double(Y(1:nSamples))' + (0:(nY-1))*nCat) = true;

outputs = false(nCat, nY);
outputs(double(YTest)' + (0:(nY-1))*nCat) = true;

plotconfusion(double(targets), double(outputs))
set(gcf, 'position', [236 241 500 500]);
shg


