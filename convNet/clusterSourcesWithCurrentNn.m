function [labels, centroids, alignedMasks] = clusterSourcesWithCurrentNn(A,convnetFn)
% Quick and dirty convenience function to cluster sources using the CNN.

if nargin<2
    convnetFn = 'SC_convNetCutSources_3class_allLayers.mat';
end

% Currently best convnet:
load(convnetFn)

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
    alignedMasks(:,:,iRoi) = tempMask;
end
    
    
%% Cluster with current NN:
Xtest = permute(alignedMasks, [1 2 4 3]);
labels = classify(convnet, Xtest);
labels = double(labels);