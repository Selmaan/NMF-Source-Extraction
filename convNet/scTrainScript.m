manLabels.labels = [];
manLabels.sources = [];
c = parcluster('local');
c.NumWorkers = 10;

%% Load in Data
includeID = find(labels<6);
labels = labels(includeID);
sources = sources(:,includeID);

manLabels.labels = cat(1,manLabels.labels,labels(:));
manLabels.sources = cat(2,manLabels.sources,sources);

%% Align and assign category

manLabels.labels(manLabels.labels<3) = 1;
manLabels.labels(manLabels.labels==3) = 2;
manLabels.labels(manLabels.labels>3) = 3;

% Align Sources
winRad = 12;
winWidth = 2*winRad+1;
nROIs = size(manLabels.sources,2);
alignedMasks = zeros(winWidth,winWidth,nROIs);
parfor nROI = 1:nROIs
    thisMask = padarray(reshape(manLabels.sources(:,nROI),512,512),[winRad,winRad],0);
    maskProps = regionprops(thisMask>0, thisMask, 'WeightedCentroid','Area');
    if length(maskProps)>1
        [~,thisComp] = max([maskProps.Area]);
        maskProps = maskProps(thisComp);
    end
    thisCent = round(maskProps.WeightedCentroid([2,1]));
    tempMask = thisMask(thisCent(1)-winRad:thisCent(1)+winRad,...
        thisCent(2)-winRad:thisCent(2)+winRad);
    alignedMasks(:,:,nROI) = tempMask;
end
imgSize = [winWidth winWidth];

%% Augment data:
realLabels = manLabels.labels;
augmentedMasks = alignedMasks;
nClasses = numel(unique(realLabels));
classSize = histcounts(realLabels);
classSizeReal = classSize;
% nRealExamples = sum(~isnan(sourceClass));
nRealExamples = length(realLabels);

realInd = {};
for i = 1:nClasses
    realInd{i} = find(realLabels==i);
end

% Determine needed class types and indeces and pre-allocate:
nTotal = size(alignedMasks,3)*25;
augmentedClass = realLabels;
augmentedInd = (1:length(realLabels))';
nPerClass = round(nTotal/nClasses);
nTotal = nPerClass*nClasses;
augmentedMasks(:,:,nTotal) = 0;
augPerClass = bsxfun(@minus,nPerClass,classSizeReal);
for iClass = 1:nClasses
    augmentedClass = cat(1,augmentedClass,ones(augPerClass(iClass),1)*iClass);
    augmentedInd = cat(1,augmentedInd,realInd{iClass}(randi(classSizeReal(iClass),augPerClass(iClass),1)));
end

% Create transformed masks
tformTemplate = affine2d(eye(3));
maxTranslPix = 3;
maxScaleFac = .2;
outRef = imref2d(size(augmentedMasks(:,:,1)));
outRef.XWorldLimits = outRef.XWorldLimits-imgSize(1)/2-0.5;
outRef.YWorldLimits = outRef.YWorldLimits-imgSize(2)/2-0.5;

augmentedMasks = alignedMasks(:,:,augmentedInd);
parfor_progress(nTotal);
parfor i = (nRealExamples + 1):nTotal
    parfor_progress;
    
    % Pick real example
%     classHere = augmentedClass(i);
%     realIndHere = augmentedInd(i);
%     realIndHere = realInd{classHere}(randi(classSizeReal(classHere)));    
    thisMask = augmentedMasks(:,:,i);

    % Create random rotation, translation, and scaling:
%     augmentedClass(i) = realLabels(realIndHere);
    x = (rand*2-1) * maxTranslPix;
    y = (rand*2-1) * maxTranslPix;
    scFac = (rand*2-1) * maxScaleFac + 1;
    
    if augmentedClass(i)==6 
        % Don't rotate cut-edge sources!
        th = 0;
    else
        th = rand * 360;
    end
    
    % Matrix representing a rotation followed by scaling followed by a translation:
    % http://planning.cs.uiuc.edu/node99.html
    tform = tformTemplate;
    tform.T = [cosd(th) -sind(th) x;
        sind(th)  cosd(th) y;
        0         0 1]';
    tform.T(1:2,1:2) = tform.T(1:2,1:2) .* scFac;
    
    augmentedMasks(:,:,i) = imwarp(thisMask, ...
        outRef, tform, 'FillValues', 0, ...
        'OutputView', outRef);
end
parfor_progress(0);

%% Adjust data format:

% Input for NN, with dimensions h-w-ch-n:
X = permute(augmentedMasks, [1 2 4 3]);
Y = categorical(augmentedClass);
nCat = numel(unique(Y));
classSizes = histcounts(Y);

% Randomize order:
ord = randperm(numel(Y));
Y = Y(ord);
X = X(:,:,:,ord);

fprintf('Class sizes: ')
fprintf('%d ', histcounts(Y));
fprintf('\n')

%% Neural network

layers = [imageInputLayer(imgSize,'Normalization','zerocenter','Name','inputl')
    convolution2dLayer(18, 32, 'Stride', 1, 'Name', 'conv1', 'Padding', 0)
    reluLayer('Name','relu1')
    crossChannelNormalizationLayer(5,'Name','cross1')
    convolution2dLayer(5, 3, 'Stride', 1, 'Name', 'conv2', 'Padding', 0)
    reluLayer('Name','relu2')
    crossChannelNormalizationLayer(5,'Name','cross2')
    fullyConnectedLayer(256, 'Name', 'full1')
    reluLayer('Name', 'relu3')
    fullyConnectedLayer(nCat, 'Name', 'full2')
    softmaxLayer('Name', 'softm')
    classificationLayer('Name', 'out')];

%% Train:
options = trainingOptions('sgdm','MaxEpochs',20,...
    'MiniBatchSize', 1e3, ...
    'L2Regularization', 0.0001, ...
    'InitialLearnRate', 0.1, ...
    'LearnRateSchedule', 'piecewise', ...
    'LearnRateDropFactor', 1, ...
    'LearnRateDropPeriod', 1);

nSamples = numel(Y);
convnet = trainNetwork(X(:,:,:,1:nSamples), Y(1:nSamples), layers, options);


% Refinement:
options = trainingOptions('sgdm','MaxEpochs',40,...
    'MiniBatchSize', 1e3, ...
    'L2Regularization', 0.0001, ...
    'InitialLearnRate', 0.03, ...
    'LearnRateSchedule', 'piecewise', ...
    'LearnRateDropFactor', 0.2, ...
    'LearnRateDropPeriod', 5, ...
    'CheckpointPath', 'C:\Users\Selmaan\Documents\GitHub\NMF-Source-Extraction\convNet');

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

%% Confusion matrix:
Xtest = permute(alignedMasks, [1 2 4 3]);
[YTest, score] = classify(convnet, Xtest);

nY = numel(realLabels);
targets = false(nCat, nY);
targets(realLabels' + (0:(nY-1))*nCat) = true;

outputs = false(nCat, nY);
outputs(double(YTest)' + (0:(nY-1))*nCat) = true;

plotconfusion(double(targets), double(outputs))