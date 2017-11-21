function stimExpt = stimCNMF(acqObj,stimBlocks,copyMemMap)


%% Parameters and Arguments

if nargin< 3
    copyMemMap = true;
end

if nargin<2
%     stimBlocks = logical([0 1 1 1 0]);
    stimBlocks = logical([1 1 1 0]);
end

cd(acqObj.defaultDir),
if exist(acqObj.indexedMovie.slice(1).channel(1).memMap,'file')
    memMap = matfile(acqObj.indexedMovie.slice(1).channel(1).memMap);
else
    error('Could not find memMap File'),
end

linFOVum = [500 500];
fprintf('Select Linear Reference Image \n'),
[linMovNames, linMovPath] = uigetfile([cd,'\*.tif'],...
    'MultiSelect','off');

if copyMemMap
    fprintf('\n Copying dsMemMap to local Drive...'),
    origPath = acqObj.indexedMovie.slice.channel.memMap;
    newPath = 'G:\tmpFiles\Slice1_dsMemMap.mat';
    copyfile(origPath,newPath),
    acqObj.indexedMovie.slice.channel.memMap = newPath;
    fprintf('Done \n'),
end
%% Get stimulation frames and target IDs
stimExpt = struct;
stimExpt.stimBlocks = stimBlocks;

syncFiles = dir(fullfile(acqObj.defaultDir,'*.abf'));
if length(syncFiles) ~= length(stimBlocks)
    error('num stimBlocks does not equal num syncFiles'),
end
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
    try
        tmp = medfilt1(syncDat(:,3:5),100) - 1.517;
    catch
        warning('This Experiment lacks Yaw Velocity'),
        tmp = medfilt1(syncDat(:,3:4),100) - 1.517;
    end
    stimExpt.ballVel{nBlock} = tmp(stimExpt.frameTimes{nBlock},:);
end

if sum(acqObj.correctedMovies.slice.channel.size(:,3))...
        == sum(cellfun(@length,stimExpt.frameTimes))
    fprintf('\n Logged Frames equals Detected Frame Times \n'),
else
    warning('Logged Frames Differs From Detected Frame Times'),
    fprintf('Manually fix frame times field: \n'),
    keyboard,
    fprintf('Correcting psych2frames: \n'),
    for nBlock = 1:length(stimExpt.syncFns)
        stimExpt.psych2frame{nBlock} = interp1(stimExpt.frameTimes{nBlock},...
            1:length(stimExpt.frameTimes{nBlock}),stimExpt.psychTimes{nBlock},'next');
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
% minStimFrame = min(cat(1,allStim{:}))-5;
% maxStimFrame = max(cat(1,allStim{:}))+5;
Y = memMap.Yr;
% Y = memMap.Yr(:,minStimFrame:maxStimFrame);
fprintf('Done! \n'),

allStimIm = nan(512^2,nTargs);
for nTarg = 1:nTargs
    s1f = allStim{nTarg};
%     s1f = allStim{nTarg}-minStimFrame+1;
    allStimIm(:,nTarg) = mean(Y(:,[s1f; s1f+1]),2)-mean(Y(:,[s1f-1; s1f-2]),2);
end
allStimIm = reshape(allStimIm,512,512,length(allStim));
% clear Y,
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

%% Create Alignment and reference structure
stimExpt.fnLin = fullfile(linMovPath,linMovNames);
[stimExpt.gLin,stimExpt.rLin,stimExpt.linRA,...
    stimExpt.gRes,stimExpt.resRA,stimExpt.resHeader,...
    stimExpt.roiCentroid,stimExpt.stimGroups] = ...
    alignStimExpt(stimExpt.resFns{2},stimExpt.fnLin);
stimExpt.xConvFactor = linFOVum(1)/stimExpt.linRA.ImageExtentInWorldX;
stimExpt.yConvFactor = linFOVum(2)/stimExpt.linRA.ImageExtentInWorldY;


%% Extract Sources w/ avgStim Image initializations
clear Y
acqObj.extractSources(1,[],winStimIm),
% acqObj.extractSources(1,reshape(Y,512,512,size(Y,2)),winStimIm),
% clear Y,
% acqObj.extractSources(1,[],winStimIm),
update_temporal_components_fromTiff(acqObj);

if copyMemMap
    acqObj.indexedMovie.slice.channel.memMap = origPath;
    delete(newPath),
end

%% Identify cells and other potential stim-sources and get traces

% Manually verify stimulation targets
[targetID,targetLabel] = confirmStimTargets(stimExpt,acqObj);
stimExpt.stimSources = targetID';
stimExpt.targetLabel = targetLabel';

A = load(acqObj.roiInfo.slice.NMF.filename,'A');
A = A.A;
l = clusterSourcesWithCurrentNn(A);
cellSources = find(l==1);
validSources{1} = [stimExpt.stimSources(targetLabel>=1);...
    setdiff(cellSources,stimExpt.stimSources(targetLabel>=1))];

stimFrames = cell(0);
for nBlock = find(stimExpt.stimBlocks)
    blockOffsetFrame = length(cat(1,stimExpt.frameTimes{1:nBlock-1}));
    stimFrames{nBlock} = blockOffsetFrame + stimExpt.psych2frame{nBlock}(1:length(stimExpt.stimOrder{nBlock}));
end

stimFrames = cat(1,stimFrames{stimExpt.stimBlocks});
stimFrames = repmat(stimFrames,1,4) + repmat(0:2:6,size(stimFrames,1),1);
stimFrames = stimFrames(:);
interpFrames = cell(0);
interpFrames{1} = stimFrames;
interpFrames{2} = stimFrames+1;
interpFrames{3} = stimFrames-1;

[dF,deconv,denoised,Gs,Lams,A,b,f] = extractTraces_NMF(acqObj,validSources,interpFrames);
dF = cell2mat(dF);
A = cell2mat(A);
denoised = cell2mat(denoised);
deconv = cell2mat(deconv);
%% 

cellFilts = bsxfun(@rdivide,A,sum(A,1));
[i,j] = ind2sub([512,512],1:512^2);
cellCentroids = ([j;i]*cellFilts)';
[xWorld,yWorld] = intrinsicToWorld(stimExpt.resRA,...
    cellCentroids(:,1),cellCentroids(:,2));
cellCentroids = [xWorld,yWorld];
nStim = size(stimExpt.roiCentroid,1);
nCell = size(cellCentroids,1);
distMat = sqrt(sum((repmat(reshape(cellCentroids,nCell,1,2),[1 nStim 1]) -...
    repmat(reshape(stimExpt.roiCentroid,1,nStim,2),[nCell 1 1])).^2,3));

stimExpt.dF_deconv = deconv;
stimExpt.dF = dF;
stimExpt.cIds = validSources{1};
stimExpt.cellFilts = cellFilts;
stimExpt.cellStimDistMat = distMat;
stimExpt.cellCentroids = cellCentroids;

tmpDir = acqObj.defaultDir;
expName = strrep(tmpDir(end-10:end-1),'\','_');
for i=1:length(stimExpt.syncFns)
    stimExpt.stimFns{i} = fullfile(acqObj.defaultDir,sprintf('%s_%d',expName,i));
    stimFID = fopen(stimExpt.stimFns{i});
    stimExpt.stimInfo{i} = fread(stimFID,'float64');
    fclose(stimFID);
end

save('stimExpt','stimExpt'),
acqObj.syncInfo.stimExpt = stimExpt;
acqObj.save,