function syncObj = bSyncIm(acqObj,wsFile,vrFile)

%% Find files if unspecified and load data

% get wavesurfer file
if ~exist('wsFile','var') || isempty(wsFile)
    [fName,fDir] = uigetfile('*.h5','Locate wavesurfer file');
    wsFile = fullfile(fDir,fName);
end

% get vermin iteration data file
if ~exist('vrFile','var') || isempty(vrFile)
    [fName,fDir] = uigetfile('*.mat','Locate vermin file');
    vrFile = fullfile(fDir,fName);
end

% load files
wsInfo = h5info(wsFile);
syncGroup = wsInfo.Groups(end).Name;
wvData = h5read(wsFile,[syncGroup '/analogScans']);
load(vrFile),

syncObj.wsFile = wsFile;
syncObj.vrFile = vrFile;
%% Identify contiguous acquisition blocks and corresponding frame times

% find all imaging frames from frame clock rising edge
frameTrig = diff(wvData(:,4)>8e3);
frameOnsets = find(frameTrig==1);
% focus blocks between acquisitions give a pair of long inter-frame-intervals
try
    focusBlocks = reshape(find(diff(frameOnsets)>1e3),2,[])';
catch
    warning(['Problem with Automatic Focus Block Detection: Using 1st and Last Points']),
    focusPoints = find(diff(frameOnsets)>1e3);
    focusBlocks = reshape(focusPoints([1,end]),2,[])';
end
nBlocks = size(focusBlocks,1) + 1;
% For each imaging block, remove the frame after first entry in pair
% through till the second entry in pair from valid frame counter
excludeInd = [];
for nBlock = 1:nBlocks-1
    excludeInd = [excludeInd, (focusBlocks(nBlock,1)+1):focusBlocks(nBlock,2)];
end
frameOnsets(excludeInd) = [];

% For each acquisition, for each slice, find contiguous acquisitions.
% Because when abort pressed frame clock continues after disk logging ends and
% different slices may have logged different # of frames, find minimum # of
% logged frames for each acquisition that are valid, and exclude all other
% frames from synchronization counter

% This loop attempts to identify continuous acquisitions from file names,
% and the largest number of frames consistent over all slices
nSlices = length(acqObj.correctedMovies.slice);
acqPref = acqObj.acqName;
for nBlock = 1:nBlocks
    thisPref = sprintf('%s_%0.5d',acqPref,nBlock);
    theseMovs =cellfun(@(x)strfind(x,thisPref),acqObj.Movies,'UniformOutput',0);
    acqBlockMembers(nBlock,:) = ~cellfun(@isempty, theseMovs);
    blockEnd = find(acqBlockMembers(nBlock,:),1,'last');
    for nSlice = 1:nSlices
        sliceFrames(nBlock,nSlice) = sum(acqObj.correctedMovies.slice(nSlice).channel(1).size(1:blockEnd,3));
    end
    if nBlock > 1
        validFrameCount(nBlock) = min(sliceFrames(nBlock,:)-sliceFrames(nBlock-1,:));
    else
        validFrameCount(nBlock) = min(sliceFrames(nBlock,:));
    end
end
% Convert from slice to unrolled volume acquisition frame number
validFrameCount = validFrameCount*(nSlices+1);

% This loop starts from the beginning of each acquisition block and gets
% rid of extra frame ticks that do not correspond to logged data, or that
% correspond to additional frames dropped for slice alignment purposes
blockFrameOffsets = [0; find(diff(frameOnsets)>1e3)];
validFrameInd = [];
for i=1:length(blockFrameOffsets)
    validFrameInd = [validFrameInd,...
        blockFrameOffsets(i)+1:blockFrameOffsets(i)+validFrameCount(i)];
end
frameOnsets(setdiff(1:length(frameOnsets),validFrameInd)) = [];

% Mark last frame for each acquisition
lastBlockFrames = find(diff(frameOnsets)>1e3);
frameTrig = zeros(size(wvData,1)-1,1);
frameTrig(frameOnsets) = 1;
IFI = diff(frameOnsets);
IFI(lastBlockFrames) = [];
fprintf('Frame Intervals: min %03d, max %03d\n',min(IFI),max(IFI)),
fprintf('%02d Imaging Blocks and %d Complete Volumes Detected\n',nBlocks,length(frameOnsets)/(nSlices+1)),
frameMidTimes = frameOnsets + round(mean(IFI)/2);

syncObj.sliceFrames = sliceFrames;
syncObj.validFrameCount = validFrameCount;
syncObj.validFrameInd = validFrameInd;
syncObj.frameOnsets = frameOnsets;
syncObj.focusBlocks = focusBlocks;
syncObj.lastBlockFrames = lastBlockFrames;
syncObj.frameMidTimes = frameMidTimes;
%% Find vermin iteration times
vermSync = wvData(:,5);
vermSessionStart = find(vermSync<-1e4,1,'last');
vermSync(1:vermSessionStart) = 0; %Set values before onset to 0
vermTrig = diff(vermSync>8e3);
vermOnsets = find(vermTrig == 1 | vermTrig == -1);
fprintf('Vermin Intervals: min: %03d, max: %03d\n',min(diff(vermOnsets)),max(diff(vermOnsets)))
vermMidTimes = (vermOnsets(2:end)+vermOnsets(1:end-1))/2;
vermTrig = zeros(size(wvData,1)-1,1);
vermTrig(vermOnsets) = 1;

syncObj.vermOnsets = vermOnsets;
syncObj.vermMidTimes = vermMidTimes;

%% Create verm2frame and frame2verm conversions

syncObj.verm2frame = interp1(1:length(frameTrig),cumsum(frameTrig==1),vermMidTimes,'nearest');
syncObj.frame2verm = interp1(1:length(vermTrig),cumsum(vermTrig==1),frameMidTimes,'nearest');

for i=1:length(syncObj.lastBlockFrames)
    syncObj.verm2frame(syncObj.verm2frame == syncObj.lastBlockFrames(i)) = nan; %Iterations lost during focus/alignment blocks
end
% if ~isempty(syncObj.lastBlockFrames)
%     syncObj.verm2frame(syncObj.verm2frame == syncObj.lastBlockFrames) = nan; %Iterations lost during focus/alignment blocks
% end
%%

nTrials = max(sessionData(end,:));
tReward = nan(nTrials,1);
tWorld = nan(nTrials,1);
trialStarts = nan(nTrials,1);
trialChoices = nan(nTrials,1);
trialEnds = nan(nTrials,1);

for nTrial = 1:nTrials
    trialOffset = find(sessionData(11,:)==nTrial,1)-1;
    
    trialIts = sessionData(:,sessionData(11,:)==nTrial);
    trialStarts(nTrial) = trialOffset + find(trialIts(8,:)==0,1,'first');
    trialChoices(nTrial) = trialOffset + find(trialIts(8,:)==-1,1,'first');
    if nTrial > 1
        preTrialInds = find(trialIts(8,:)==1);
        if tReward(nTrial-1)==1
            itiDur = 2;
        elseif tReward(nTrial-1) == 0
            itiDur = 4;
        end
        itiStart = find(cumsum(trialIts(10,preTrialInds(2:end)))>itiDur,1,'first');
        if ~isempty(itiStart)
            trialEnds(nTrial-1) = trialOffset + itiStart;
        else
            trialEnds(nTrial-1) = nan;
        end
    end  
    tReward(nTrial) = sum(trialIts(9,:));
    tWorld(nTrial)= trialIts(1,find(trialIts(8,:)==-1,1));
end

syncObj.tReward = tReward;
syncObj.tWorld = tWorld;
syncObj.trialStarts = trialStarts;
syncObj.trialChoices = trialChoices;
syncObj.trialEnds = trialEnds;
