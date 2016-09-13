function [binData,trialInfo, vrData] = binAlignNeurons(neurData, syncObj)

%% Binning Parameters
% For now these are set to default values
% nPosBins = 30;
% posBins = linspace(5,205,nPosBins+1); %Edges for 30 spatial bins
% posBins(1)=0;
preTrialShifts = -3:-1;
postTrialShifts = 0:23;

trialInfo.preTrialShifts = preTrialShifts;
trialInfo.postTrialShifts = postTrialShifts;
%% Create trial info structure
trialInfo.missingTrials = union(find(isnan(syncObj.verm2frame(syncObj.trialStarts))),find(isnan(syncObj.verm2frame(syncObj.trialChoices))));
nTrials = min(length(syncObj.trialStarts),length(syncObj.trialChoices));
trialInfo.validTrials = setdiff(1:nTrials,trialInfo.missingTrials);
trialInfo.worldVec = syncObj.tWorld(trialInfo.validTrials);

trialInfo.wTrials = find(mod(trialInfo.worldVec,2)==1); %white cue
trialInfo.bTrials = find(mod(trialInfo.worldVec,2)==0); %black cue
trialInfo.kTrials = find(trialInfo.worldVec<5);         %checker trial
trialInfo.nTrials = find(trialInfo.worldVec>=5);        %no checker trial
trialInfo.lTrials = find(mod(trialInfo.worldVec,4)<2);  %left chamber
trialInfo.rTrials = find(mod(trialInfo.worldVec,4)>=2); %right chamber
trialInfo.cTrials = find(syncObj.tReward(trialInfo.validTrials) == 1); %Correct Trials
trialInfo.iTrials = find(syncObj.tReward(trialInfo.validTrials) == 0); %Incorrect Trials

%%
trialDurs = (syncObj.vermOnsets(syncObj.trialChoices)-syncObj.vermOnsets(syncObj.trialStarts))/1e3;
nPosBins = round(median(trialDurs)*6);
load(syncObj.vrFile),
positionIterations = sessionData([6,11],sessionData(8,:)==0);

posPrctiles = nan(nTrials,nPosBins+1);
for nTrial = 1:nTrials
    posPrctiles(nTrial,:) = prctile(positionIterations(1,positionIterations(2,:)==nTrial),...
        linspace(0,100,nPosBins+1));
end

posBins = median(posPrctiles);
posBins(1) = 0;
trialInfo.posBins = posBins;

%%
vT = trialInfo.validTrials;
posData = nan(size(neurData,1),length(posBins)-1,length(vT));
preData = nan(size(neurData,1),length(preTrialShifts),length(vT));
postData = nan(size(neurData,1),length(postTrialShifts),length(vT));

posVR = nan(size(sessionData,1),length(posBins)-1,length(vT));
preVR = nan(size(sessionData,1),length(preTrialShifts),length(vT));
postVR = nan(size(sessionData,1),length(postTrialShifts),length(vT));

for iTrial = 1:length(vT)
    nTrial = vT(iTrial);
    trialOffset = syncObj.trialStarts(nTrial)-1;
    trialData = sessionData(:,syncObj.trialStarts(nTrial):syncObj.trialChoices(nTrial));
    
    for posBin = 1:length(posBins)-1
        binInd = find(trialData(6,:)>posBins(posBin) & trialData(6,:)<= posBins(posBin+1));
        if isempty(binInd)
            binInd = find(trialData(6,:)>posBins(posBin),1,'first');
        end
        posVR(:,posBin,iTrial) = mean(sessionData(:,binInd + trialOffset),2);
%         posVR(10,posBin,iTrial) = sum(sessionData(10,binInd + trialOffset),2);
        posVR(10,posBin,iTrial) = syncObj.vermOnsets(binInd(end)+trialOffset+1) ...
            - syncObj.vermOnsets(binInd(1) + trialOffset);
        frameInds = ceil(syncObj.verm2frame(binInd + trialOffset)/5);
        posData(:,posBin,iTrial) = mean(neurData(:,frameInds),2);
    end
    
    for preBin = 1:length(preTrialShifts)
        thisShift = preTrialShifts(preBin);
        frameShifted = thisShift+ceil(syncObj.verm2frame(syncObj.trialStarts(nTrial))/5);
        preData(:,preBin,iTrial) = neurData(:,frameShifted);
        vermInd = syncObj.frame2verm(frameShifted*5)+(-1:1);
        vermInd(vermInd<1) = 1;
        preVR(:,preBin,iTrial) = mean(sessionData(:,vermInd),2);
        preVR(10,preBin,iTrial) = syncObj.frameOnsets(5*(frameShifted+1))-syncObj.frameOnsets(5*frameShifted);
    end
    
%     preData(:,:,iTrial) = neurData(:,preTrialShifts+ceil(syncObj.verm2frame(syncObj.trialStarts(nTrial))/5));
    for postBin = 1:length(postTrialShifts)
        thisShift = postTrialShifts(postBin);
        frameShifted = thisShift+ceil(syncObj.verm2frame(syncObj.trialChoices(nTrial))/5);
        postData(:,postBin,iTrial) = neurData(:,frameShifted);
        vermInd = syncObj.frame2verm(frameShifted*5)+(-1:1);
        vermInd(vermInd>size(sessionData,2)) = size(sessionData,2);
        postVR(:,postBin,iTrial) = mean(sessionData(:,vermInd),2);
        postVR(10,postBin,iTrial) = syncObj.frameOnsets(5*(frameShifted+1))-syncObj.frameOnsets(5*frameShifted);
    end

% postData(:,:,iTrial) = neurData(:,postTrialShifts+ceil(syncObj.verm2frame(syncObj.trialChoices(nTrial))/5));
end

binData = cat(2,preData,posData,postData);
vrData = cat(2,preVR,posVR,postVR);    

%%
choiceBin = length(trialInfo.preTrialShifts)+length(trialInfo.posBins)+3;
trialInfo.rChoice = find(squeeze(vrData(5,choiceBin,:))<0);
trialInfo.lChoice = find(squeeze(vrData(5,choiceBin,:))>0);