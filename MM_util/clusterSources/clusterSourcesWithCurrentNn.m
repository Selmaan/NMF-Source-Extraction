function [labels, centroids] = clusterSourcesWithCurrentNn(A)
% Quick and dirty convenience function to cluster sources using the CNN.

% Currently best convnet:
load('D:\GitHub\other\NMF-Source-Extraction\MM_util\clusterSources\conventThatCanRecorgnizeCutSources.mat')

%% Create source patches:
% Align Sources
winRad = 12;
winWidth = 2*winRad+1;
nRois = size(A,2);
alignedMasks = zeros(winWidth,winWidth,nRois);
centroids = zeros(nRois, 2);
for iRoi = 1:nRois
    thisMask = padarray(reshape(A(:,iRoi),512,512),[winRad,winRad],0);
    maskProps = regionprops(thisMask>0, thisMask, 'WeightedCentroid','Area');
    if length(maskProps)>1
        [~,thisComp] = max([maskProps.Area]);
        maskProps = maskProps(thisComp);
    end
    centroids(iRoi, :) = round(maskProps.WeightedCentroid([2,1]));
    tempMask = thisMask(centroids(iRoi, 1)-winRad:centroids(iRoi, 1)+winRad,...
        centroids(iRoi, 2)-winRad:centroids(iRoi, 2)+winRad);
    tempMask = medfilt2(tempMask,[3 3]);
    alignedMasks(:,:,iRoi) = tempMask;
end



%% Cluster with current NN:
Xtest = permute(alignedMasks, [1 2 4 3]);
labels = classify(convnet, Xtest);
% YTest = addcats(YTest, 'cutEdge');
% YTest = renamecats(YTest, {'doughnutSoma', 'outOfPlaneSoma', ...
%                 'smallRoundProcess', 'complexProcess', ...
%                 'messyAndMultiples', 'cutEdge'});
labels = double(labels);
classes = unique(labels);

cIm = zeros(512,512,3);
cIm(:,:,1) = mat2gray(reshape(full(sum(A(:,labels==classes(2)),2)),512,512));
cIm(:,:,2) = mat2gray(reshape(full(sum(A(:,labels==classes(1)),2)),512,512));
cIm(:,:,3) = mat2gray(reshape(full(sum(A(:,labels==classes(3)|labels==classes(4)),2)),512,512), [0, 0.3]);

figure(44563456)
imshow(cIm, 'initialmag', 'fit')
title('Good sources')

cIm = zeros(512,512,3);
cIm(:,:,1) = mat2gray(reshape(full(sum(A(:,labels==classes(5)),2)),512,512), [0, 0.3]);
cIm(:,:,3) = mat2gray(reshape(full(sum(A(:,labels==classes(6)),2)),512,512), [0, 0.3]);

figure(44563457)
imshow(cIm, 'initialmag', 'fit')
title('Bad sources')