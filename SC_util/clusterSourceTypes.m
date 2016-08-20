function sourceProps = clusterSourceTypes(A,nClusters)

if ~exist('nClusters','var') || isempty(nClusters)
    nClusters = 3;
end

%% Parameters
winRad = 15;

%% Normalize Sources
nA = full(sqrt(sum(A.^2)));
nr = length(nA);
A = full(A/spdiags(nA(:),0,nr,nr));

%% Align Sources
winWidth = 2*winRad+1;
nROIs = size(A,2);
alignedMasks = zeros(winWidth,winWidth,nROIs);
allCentroids = nan(nROIs,2);
% allSkews = nan(nROIs,1);
for nROI = 1:nROIs
    thisMask = padarray(reshape(A(:,nROI),512,512),[winRad,winRad],0);
    maskProps = regionprops(thisMask>0, thisMask, 'WeightedCentroid','Area');
    if length(maskProps)>1
        [~,thisComp] = max([maskProps.Area]);
        maskProps = maskProps(thisComp);
    end
    allCentroids(nROI,:) = maskProps.WeightedCentroid;
    thisCent = round(maskProps.WeightedCentroid([2,1]));
    tempMask = thisMask(thisCent(1)-winRad:thisCent(1)+winRad,...
        thisCent(2)-winRad:thisCent(2)+winRad);
    tempMask = medfilt2(tempMask,[3 3]);
    alignedMasks(:,:,nROI) = tempMask;
end

%% PCA + gmm clustering
pcaDims = 2;
[mEig,mSco,mLat] = pca(reshape(alignedMasks,winWidth^2,nROIs)');
clusterScores = mSco(:,1:pcaDims);
% clusterScores = tsne(reshape(alignedMasks,31^2,nROIs)',[],pcaDims);
gmModel = fitgmdist(clusterScores,nClusters,'Start','plus','Replicates',10);
[~,pID] = max(gmModel.mu(:,1));
remIDs = setdiff(1:nClusters,pID);
[~,cID] = max(gmModel.mu(remIDs,2));
cID = remIDs(cID);
jID = setdiff(1:nClusters,[pID,cID]);
[idx,~,gmProbs] = cluster(gmModel,clusterScores);

cIm = zeros(512,512,3);
cIm(:,:,1) = reshape(sum(A(:,idx==pID),2),512,512);
cIm(:,:,2) = reshape(sum(A(:,idx==cID),2),512,512);
for i=1:length(jID)
    cIm(:,:,3) = cIm(:,:,3) + reshape(sum(A(:,idx==jID(i)),2),512,512);
end

%% Output
sourceProps.allCentroids = allCentroids;
sourceProps.mEig = mEig(:,1:pcaDims);
sourceProps.mSco = mSco(:,1:pcaDims);
sourceProps.mLat = mLat(1:pcaDims)/sum(mLat);
sourceProps.gmModel = gmModel;
sourceProps.cID = cID;
sourceProps.jID = jID;
sourceProps.pID = pID;
sourceProps.idx = idx;
sourceProps.gmProbs = gmProbs;
sourceProps.cIm = cIm;

%% Display

figure,imshow(cIm*5),
for i=1:3
    figure,imshow(cIm(:,:,i)*5),
end