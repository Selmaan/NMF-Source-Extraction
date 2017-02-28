function stimExpt = stimCNMF(acqObj,stimBlocks)

if nargin<2
    stimBlocks = logical([0 1 1 1 0]);
end

cd(acqObj.defaultDir),
if exist(acqObj.indexedMovie.slice(1).channel(1).memMap,'file')
    memMap = matfile(acqObj.indexedMovie.slice(1).channel(1).memMap);
else
    error('Could not find memMap File'),
end

%% Get stimulation frames and target IDs
stimExpt = struct;
stimExpt.stimBlocks = stimBlocks;

syncFiles = dir(fullfile(acqObj.defaultDir,'*.abf'));
for i=1:length(syncFiles)
    stimExpt.syncFns{i} = fullfile(acqObj.defaultDir,syncFiles(i).name);
    stimExpt.resFns{i} = fullfile(acqObj.defaultDir,...
        sprintf('resFOV1_0000%d_00001.tif',i));
end

for nBlock = 1:length(stimExpt.syncFns)
    syncDat = abfload(stimExpt.syncFns{nBlock});
    stimExpt.frameTimes{nBlock} = find(syncDat(2:end,1)>1 & syncDat(1:end-1,1)<1);
    stimExpt.psychTimes{nBlock} = find(syncDat(2:end,2)>1 & syncDat(1:end-1,2)<1);
    stimExpt.psych2frame{nBlock} = interp1(stimExpt.frameTimes{nBlock},...
        1:length(stimExpt.frameTimes{nBlock}),stimExpt.psychTimes{nBlock},'next');
    if stimExpt.stimBlocks(nBlock)
        thisHeader = scanimage.util.opentif(stimExpt.resFns{nBlock});
        stimExpt.stimOrder{nBlock} = thisHeader.SI.hPhotostim.sequenceSelectedStimuli-1;
    else
        stimExpt.stimOrder{nBlock} = [];
    end
end
clear syncDat

%% Get Stimulation Frames for all Targets

nTargs = length(unique(stimExpt.stimOrder{find(stimExpt.stimBlocks,1)}));
allStim = cell(nTargs,1);

for nBlock = find(stimExpt.stimBlocks)
    for nTarg = 1:nTargs
        blockOffsetFrame = length(cat(1,stimExpt.frameTimes{1:nBlock-1}));
        theseStim = blockOffsetFrame+stimExpt.psych2frame{nBlock}(stimExpt.stimOrder{nBlock}==nTarg);
        allStim{nTarg} = [allStim{nTarg};ceil(theseStim/memMap.dsRatio)];
    end
end

%% Load Dataset into memory and create mean images
fprintf('Loading down-sampled data into memory...'),
Y = memMap.Yr;
fprintf('Done! \n'),

allStimIm = nan(512^2,nTargs);
for nTarg = 1:nTargs
    s1f = allStim{nTarg};
    allStimIm(:,nTarg) = mean(Y(:,[s1f; s1f+1]),2)-mean(Y(:,[s1f-1; s1f-2]),2);
end
allStimIm = reshape(allStimIm,512,512,length(allStim));
%% Get Scanfield and Target Coordinates
[hRoiGroup,stimGroups] = scanimage.util.readTiffRoiData(stimExpt.resFns{find(stimExpt.stimBlocks,1)});
scanfield = hRoiGroup.rois(1).scanfields(1);
resRA = imref2d([512 512],[scanfield.rect(1),scanfield.rect(1)+scanfield.rect(3)],...
    [scanfield.rect(2),scanfield.rect(2)+scanfield.rect(4)]);
for nROI = 2:length(stimGroups)
    roiCentroid(nROI-1,:) = stimGroups(nROI).rois(2).scanfields.centerXY;
end

%% Post process stimAvg images

normStimIm = bsxfun(@rdivide,allStimIm,sqrt(meanRef(acqObj)));

[xIntrinsic,yIntrinsic] = worldToIntrinsic(resRA,roiCentroid(:,1),roiCentroid(:,2));
winStimIm = zeros(size(normStimIm));
for nTarg = 1:nTargs
    xWin = min(512,max(1,round((-60:60)+xIntrinsic(nTarg))));
    yWin = min(512,max(1,round((-60:60)+yIntrinsic(nTarg))));
    thisWin = normStimIm(yWin,xWin,nTarg);
    thisWin(thisWin<prctile(thisWin(:),90)) = 0;
    winStimIm(yWin,xWin,nTarg) = thisWin;
end
    
stimExpt.rawStimIm = allStimIm;
stimExpt.procStimIm = winStimIm;

%% Extract Sources w/ avgStim Image initializations
acqObj.extractSources(1,reshape(Y,512,512,size(Y,2)),winStimIm),

