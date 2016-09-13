%% Load Data
D(1) = load('patchResults_slice01_v0820.mat', 'A','C');
D(2) = load('patchResults_slice02_v0820.mat', 'A','C');
D(3) = load('patchResults_slice03_v0820.mat', 'A','C');
D(4) = load('patchResults_slice04_v0820.mat', 'A','C');

As = []; Cs = [];
for nSlice=1:4
    As = [As, D(nSlice).A];
    thisC = [];
    for nBlock = 1:size(syncObj.sliceFrames,1)
        if nBlock > 1
            offsetInd = syncObj.sliceFrames(nBlock-1,nSlice);
        else
            offsetInd = 0;
        end
        blockLength = syncObj.validFrameCount(nBlock)/5;
        thisC = [thisC, D(nSlice).C(:,offsetInd+1:offsetInd+blockLength)];
    end
    Cs = [Cs; thisC];
end

% sourceProps = clusterSourceTypes(As);
sourceClusterGUI(As,3);
clear D
%% Cluster Cell Sources
% cIdx = sourceProps.gmProbs(:,sourceProps.pID)>0.1;
cIdx = guiClusterIDs;
cellFilts = As(:,cIdx);
imshow(reshape(full(sum(cellFilts,2)),512,512)*3),
%% Detrend and Deconvolve Traces, save output
acqBlocks = syncObj.sliceFrames(:,end);
cellCalcium = Cs(cIdx,:);
cellDenoised = nan(size(cellCalcium));
cellDeconv = nan(size(cellCalcium));
cellBaseline = nan(size(cellCalcium,1),length(acqBlocks));
cellG = nan(size(cellCalcium,1),length(acqBlocks));
cellNoise = nan(size(cellCalcium,1),length(acqBlocks));
parfor_progress(size(cellCalcium,1)*length(acqBlocks));
for nBlock = 1:length(acqBlocks)
        if nBlock > 1
            inds = 1+acqBlocks(nBlock-1):acqBlocks(nBlock);
        else
            inds = 1:acqBlocks(nBlock);
        end
        
        parfor nCell = 1:size(cellCalcium,1)
            thisTrace = removeSourceBaseline(cellCalcium(nCell,inds));
            cellCalcium(nCell,inds) = thisTrace;
            [cDe,bs,~,g,sn,sp] = constrained_foopsi(thisTrace);
            cellDenoised(nCell,inds) = cDe+bs;
            cellDeconv(nCell,inds) = sp;
            cellBaseline(nCell,nBlock) = bs;
            cellG(nCell,nBlock) = g(1);
            cellNoise(nCell,nBlock) = sn;
            parfor_progress;
        end
end
parfor_progress(0);

save('cellData_0902','cellDeconv','cellCalcium',...
    'syncObj','cellDenoised','cellFilts','cellG','cellNoise')