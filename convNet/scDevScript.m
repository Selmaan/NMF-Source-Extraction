%% Load Data
load('Z:\HarveyLab\Tier1\Shih_Yi\Imaging\9\170810\Slice02_patchResults_v170311.mat','A')
memMap = matfile('Z:\HarveyLab\Tier1\Shih_Yi\Imaging\9\170810\Slice2_MemMapReconstructed.mat');
Y = memMap.Y(:,:,1:1000);

%% Prepare GUI and Labels

A = A ./ sqrt(sum(A.^2));

winRad = 12;
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
end

l = clusterSourcesWithCurrentNn(A);
o = classifierGui(A,l,allCentroids);

%%
