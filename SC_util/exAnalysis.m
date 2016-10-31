%%
% normTraces = bsxfun(@rdivide,abs(de_),mean(abs(cellDeconv),2));
[bD,tI,vD] = binAlignNeurons(de_cell,FOV1.syncInfo);
[maxVal,maxInd] = max(mean(bD,3),[],2);
[~,pkSrt] = sort(maxInd,'ascend');
figure,plot(squeeze(median(vD(10,:,:),3)))

if mod(tI.worldVec(1),4) == 1 | mod(tI.worldVec(1),4) == 2
    switchPoints(1) = find(abs(mod(tI.worldVec,4)-1.5)>1,1,'first');
    switchPoints(2) = find(abs(mod(tI.worldVec,4)-1.5)>1,1,'last');
else
    switchPoints(1) = find(abs(mod(tI.worldVec,4)-1.5)<1,1,'first');
    switchPoints(2) = find(abs(mod(tI.worldVec,4)-1.5)<1,1,'last');
end

X=zeros(length(tI.validTrials),7)-1;
X(tI.wTrials,1) = 1;
% X(tI.rTrials,2) = 1/2; %use 'correct' choice
X(tI.rChoice,2) = 1; %use actual choice
X(tI.cTrials,3) = 1;
X(tI.kTrials,4) = 1;
X(switchPoints(1):switchPoints(2),5) = 1;
X(:,6) = 1;
X(:,7) = linspace(-1,1,length(tI.validTrials));

invX = pinv(X);
bMat = nan(size(bD,1),size(bD,2),size(X,2));
for nNeuron = 1:size(bD,1)
    neurCoefs = (invX*squeeze(bD(nNeuron,:,:))')';
    bMat(nNeuron,:,:) = neurCoefs;
%     parfor ind = 1:size(bD,2)
%         bMat(nNeuron,ind,:) = regress(squeeze(bD(nNeuron,ind,:)),X);
%     end
end

figure,hold on
for i=1:5
plot(mean(abs(bMat(:,:,i)))./mean(abs(bMat(:,:,6))))
end
axis tight
ay = ylim;
line((length(tI.preTrialShifts))*[1 1], [ay(1) ay(2)],'color','k'),
line((length(tI.preTrialShifts)+find(tI.posBins>150,1))*[1 1],...
    [ay(1) ay(2)],'color','k','linestyle','--'),
line((length(tI.preTrialShifts)+length(tI.posBins-1))*[1 1],...
    [ay(1) ay(2)],'color','k'),
line((length(tI.preTrialShifts)+length(tI.posBins-1)+12)*[1 1],...
    [ay(1) ay(2)],'color','k','linestyle','--'),
line((length(tI.preTrialShifts)+find(tI.posBins>90,1))*[1 1],...
    [ay(1) ay(2)],'color','g','linestyle','--'),
figure,plot(mean(abs(bMat(:,:,6))))

fNorm = sum(f)./mean(sum(f,1),2);
fD = binAlignNeurons(fNorm,syncObj);
fMat = nan(size(bD,2),size(X,2));
for ind = 1:size(bD,2)
    fMat(ind,:) = regress(squeeze(fD(1,ind,:)),X);
end
figure,plot(fMat(:,1:5)./std(fMat(:,6))),
figure,plot(fMat(:,6)),

%% Trial Modulation Analysis
normTraces = bsxfun(@rdivide,cellDeconv,mean(cellDeconv,2));
onOffsets = -3:18;
choiceOffsets = -12:24;

sliceOnsets = ...
    ceil(syncObj.verm2frame(syncObj.trialStarts)/5);
sliceOnsets(isnan(sliceOnsets)) = [];
sliceChoices = ...
    ceil(syncObj.verm2frame(syncObj.trialChoices)/5);
sliceChoices(isnan(sliceChoices)) = [];

neurAvg = nan(size(normTraces,1),...
    length(onOffsets)+length(choiceOffsets));

for ind = 1:length(onOffsets)
    neurAvg(:,ind) = ...
        mean(normTraces(:,onOffsets(ind)+sliceOnsets),2);
end

for ind = 1:length(choiceOffsets)
    neurAvg(:,ind+length(onOffsets)) = ...
        mean(normTraces(:,choiceOffsets(ind)+sliceChoices),2);
end


numShufs = 1e3;
neurShufs = nan(size(normTraces,1),numShufs);
for nShuf = 1:numShufs
    neurShufs(:,nShuf) = mean(normTraces(:,...
        randi(size(normTraces,2),1,length(sliceOnsets))),2);
end
normAvg = bsxfun(@minus,neurAvg,mean(neurShufs,2));
normAvg = bsxfun(@rdivide,normAvg,std(neurShufs,[],2));
[pkVal,pkInd] = max(normAvg,[],2);
[~,pkSrt] = sort(pkInd,'ascend');

%% Trial Modulation Plots
sigPks = pkSrt(pkVal(pkSrt)>10);
figure,
imagesc(normAvg(sigPks,:),[0 40]),
figure,
imagesc(bsxfun(@rdivide,normAvg(sigPks,:),max(normAvg(sigPks,:),[],2))),

figure,plot(mean(normAvg)),
axis tight,
ay = ylim;
tOnset = find(onOffsets==0);
alignBoundary = length(onOffsets)+1;
choiceOnset = find(choiceOffsets==0) + length(onOffsets);
itiOnset = find(choiceOffsets==13) + length(onOffsets);
line(tOnset*[1 1],[ay(1) ay(2)],'color','k'),
line(alignBoundary*[1 1],[ay(1) ay(2)],'color','k','linestyle',':'),
line(choiceOnset*[1 1],[ay(1) ay(2)],'color','k'),
line(itiOnset*[1 1],[ay(1) ay(2)],'color','k','linestyle',':'),

figure,imagesc(corrcoef(normAvg)),
% figure,imagesc(cov(normAvg)),
%% Condition Modulation Analysis
normTraces = bsxfun(@rdivide,cellDeconv,mean(cellDeconv,2));
onOffsets = -3:18;
choiceOffsets = -12:24;
neurAvg = nan(size(normTraces,1),...
    length(onOffsets)+length(choiceOffsets),4);
normAvg = neurAvg;
nShuffles = 1e3;
shufAvg = nan(size(normTraces,1),length(onOffsets)+length(choiceOffsets),nShuffles);

for cond = 1:4
    theseTrials = find((syncObj.tWorld == cond+4) & (syncObj.tReward == 1));
    
    sliceOnsets = ...
        ceil(syncObj.verm2frame(syncObj.trialStarts(theseTrials))/5);
    sliceOnsets(isnan(sliceOnsets)) = [];
    sliceChoices = ...
        ceil(syncObj.verm2frame(syncObj.trialChoices(theseTrials))/5);
    sliceChoices(isnan(sliceChoices)) = [];
    
    for ind = 1:length(onOffsets)
        neurAvg(:,ind,cond) = ...
            mean(normTraces(:,onOffsets(ind)+sliceOnsets),2);
    end

    for ind = 1:length(choiceOffsets)
        neurAvg(:,ind+length(onOffsets),cond) = ...
            mean(normTraces(:,choiceOffsets(ind)+sliceChoices),2);
    end
    
%     numShufs = 1e3;
%     neurShufs = nan(size(normTraces,1),numShufs);
%     for nShuf = 1:numShufs
%         neurShufs(:,nShuf) = mean(normTraces(:,...
%             randi(size(normTraces,2),1,length(sliceOnsets))),2);
%     end
%     
%     normAvg(:,:,cond) = bsxfun(@minus,neurAvg(:,:,cond),mean(neurShufs,2));
%     normAvg(:,:,cond) = bsxfun(@rdivide,normAvg(:,:,cond),std(neurShufs,[],2));
end

for nShuffle = 1:nShuffles
    theseTrials = randperm(length(syncObj.trialStarts),round(length(syncObj.trialStarts)/4));
    
    sliceOnsets = ...
        ceil(syncObj.verm2frame(syncObj.trialStarts(theseTrials))/5);
    sliceOnsets(isnan(sliceOnsets)) = [];
    sliceChoices = ...
        ceil(syncObj.verm2frame(syncObj.trialChoices(theseTrials))/5);
    sliceChoices(isnan(sliceChoices)) = [];
    
    for ind = 1:length(onOffsets)
        shufAvg(:,ind,nShuffle) = ...
            mean(normTraces(:,onOffsets(ind)+sliceOnsets),2);
    end

    for ind = 1:length(choiceOffsets)
        shufAvg(:,ind+length(onOffsets),nShuffle) = ...
            mean(normTraces(:,choiceOffsets(ind)+sliceChoices),2);
    end
end

normAvg = bsxfun(@minus,neurAvg,mean(shufAvg,3));
normAvg = bsxfun(@rdivide,normAvg,std(shufAvg,[],3));
normAvg(isnan(normAvg))=0;
[pkVal,pkInd] = max(mean(shufAvg,3),[],2);
[~,pkSrt] = sort(pkInd,'ascend');

%% Condition Modulation Plots
cueDiff = mean(normAvg(:,:,[1 3]),3)-mean(normAvg(:,:,[2 4]),3);
turnDiff = mean(normAvg(:,:,[1 4]),3)-mean(normAvg(:,:,[2 3]),3);
figure,imagesc(cueDiff(pkSrt,:),[-15 15])
figure,imagesc(turnDiff(pkSrt,:),[-15 15])
% figure,imagesc(cov(cueDiff))
% figure,imagesc(cov(turnDiff))
figure,imagesc(sqrtCov(cueDiff)),
figure,imagesc(sqrtCov(turnDiff)),
figure,hold on
plot(diag(sqrtCov(cueDiff)),'linewidth',2),
plot(diag(sqrtCov(turnDiff)),'linewidth',2),
axis tight,
ay = ylim;
tOnset = find(onOffsets==0);
alignBoundary = length(onOffsets)+1;
choiceOnset = find(choiceOffsets==0) + length(onOffsets);
itiOnset = find(choiceOffsets==13) + length(onOffsets);
line(tOnset*[1 1],[ay(1) ay(2)],'color','k'),
line(alignBoundary*[1 1],[ay(1) ay(2)],'color','k','linestyle',':'),
line(choiceOnset*[1 1],[ay(1) ay(2)],'color','k'),
line(itiOnset*[1 1],[ay(1) ay(2)],'color','k','linestyle',':'),
%% NMF Visualization
[w,h] = nnmf(normTraces',10);
[~,maxH] = max(h);
[~,displayInd] = sort(maxH,'ascend');
figure,
imagesc(matConv(normTraces(displayInd,:),2),[0 6]),

%% Single-Trial Analyses
normTraces = bsxfun(@rdivide,cellDeconv,mean(cellDeconv,2));

sliceOnsets = ...
    ceil(syncObj.verm2frame(syncObj.trialStarts)/5);
onsetMissing = find(isnan(sliceOnsets));
sliceChoices = ...
    ceil(syncObj.verm2frame(syncObj.trialChoices)/5);
choiceMissing = find(isnan(sliceChoices));
missingTrials = union(onsetMissing,choiceMissing);
validTrials = setdiff(1:length(syncObj.trialStarts),missingTrials);

neurResp = nan(size(normTraces,1),...
    length(validTrials),3);

for iTrial = 1:length(validTrials)
    nTrial = validTrials(iTrial);
    neurResp(:,iTrial,1) = mean(normTraces(:,sliceOnsets(nTrial)+(0:6)),2);
    neurResp(:,iTrial,2) = mean(normTraces(:,sliceChoices(nTrial)+(-6:0)),2);
    neurResp(:,iTrial,3) = mean(normTraces(:,sliceChoices(nTrial)+(2:14)),2);
end

%% Single Trial Plots

wTrialsC = find(syncObj.tReward(validTrials) == 1 & mod(syncObj.tWorld(validTrials),2)==0);
bTrialsC = find(syncObj.tReward(validTrials) == 1 & mod(syncObj.tWorld(validTrials),2)==1);
wTrialsI = find(syncObj.tReward(validTrials) == 0 & mod(syncObj.tWorld(validTrials),2)==0);
bTrialsI = find(syncObj.tReward(validTrials) == 0 & mod(syncObj.tWorld(validTrials),2)==1);
cueTrials = {wTrialsC,bTrialsC,wTrialsI,bTrialsI};
rTrialsC = find(syncObj.tReward(validTrials) == 1 & abs(syncObj.tWorld(validTrials)-6.5) == 3/2);
rTrialsI = find(syncObj.tReward(validTrials) == 0 & abs(syncObj.tWorld(validTrials)-6.5) == 1/2);
lTrialsC = find(syncObj.tReward(validTrials) == 1 & abs(syncObj.tWorld(validTrials)-6.5) == 1/2);
lTrialsI = find(syncObj.tReward(validTrials) == 0 & abs(syncObj.tWorld(validTrials)-6.5) == 3/2);
turnTrials = {rTrialsC,lTrialsC,rTrialsI,lTrialsI};

Y = mdscale(pdist(neurResp(:,:,3)','spearman'),2);
figure,hold on
for i=1:length(turnTrials)
    ind = turnTrials{i};
    plot(Y(ind,1),Y(ind,2),'.','markersize',10),
end