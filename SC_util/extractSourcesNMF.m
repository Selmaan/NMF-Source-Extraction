function extractSourcesNMF(acqObj,nSlice,data,initImages)

syncObj = acqObj.syncInfo;
acqBlocks = [1 syncObj.sliceFrames(1,nSlice)];
for blockNum = 2:size(syncObj.sliceFrames,1)
    acqBlocks(blockNum,:) = ...
        [1+syncObj.sliceFrames(blockNum-1,nSlice), syncObj.sliceFrames(blockNum,nSlice)];
end

memMap = matfile(acqObj.indexedMovie.slice(nSlice).channel.memMap);  %useful to construct memmap even w/ data
if isempty(data) %using memory map
    imSize = acqObj.correctedMovies.slice(nSlice).channel.size(1,1:2);
    nFramesDS = size(memMap,'Y',3);
else
    imSize = [size(data,1),size(data,2)];
    nFramesDS = size(data,3);
end

options = CNMFSetParms(...
    'd1',imSize(1),'d2',imSize(2),...
    'spatial_method','constrained',...
    'search_method','dilate','se',strel('disk',2,0),...
    'deconv_method','constrained_foopsi',...    % activity deconvolution method
    'merge_thr',0.8,...                    % merging threshold
    'nB',3,...
    'thr_method','nrg',...
    'nrgthr',0.99,...
    'clos_op',strel('square',1),...
    'medw',[1 1]);

% Use a separate options structure for when thresholding/postprocessing is
% desired, as opposed to when it's automatically called by spatial updates
optionsThresh = options;
optionsThresh.medw = [3 3];
optionsThresh.clos_op = strel('square',3);
optionsThresh.nrgthr = 0.95;

%% Get Patches Results
patch_size = [52,52];                   % size of each patch along each dimension (optional, default: [32,32])
overlap = [6,6];                        % amount of overlap in each dimension (optional, default: [4,4])
nFactors = 12;
patches = construct_patches(imSize,patch_size,overlap);

% parfor_progress(length(patches));
% parfor_progress;
fprintf('Initializing patches...');

% Commented out code uses initImages during initialization phase
% if isempty(data) %memory mapped
%     RESULTS = patchInitNMF(acqObj,nSlice,patches,1,nFactors,initImages);
%     parfor patchNum = 2:length(patches)
%     %     parfor_progress;
%         RESULTS(patchNum) = patchInitNMF(acqObj,nSlice,patches,patchNum,nFactors,initImages);
%     end
% else
%     for patchNum = 1:length(patches)
%         RESULTS(patchNum) = patchInitNMF(data,nSlice,patches,patchNum,nFactors,initImages);
%     end
% end

% Code used here does not use initImages, it'll be added after initialization
if isempty(data) %memory mapped
    RESULTS = patchInitNMF(acqObj,nSlice,patches,1,nFactors);
    parfor patchNum = 2:length(patches)
    %     parfor_progress;
        RESULTS(patchNum) = patchInitNMF(acqObj,nSlice,patches,patchNum,nFactors);
    end
else
    for patchNum = 1:length(patches)
        RESULTS(patchNum) = patchInitNMF(data,nSlice,patches,patchNum,nFactors);
    end
end
fprintf(' done. \n');
%% combine results into one structure
fprintf('Combining results from different patches...');
C = double(cell2mat({RESULTS(:).C}'));
d = prod(imSize);
A = sparse(d,size(C,1));

cnt = 0;
B = sparse(prod(imSize),length(patches));
MASK = zeros(imSize);
F = zeros(length(patches),nFramesDS);
for i = 1:length(patches)
    for k = 1:size(RESULTS(i).A,2)
            cnt = cnt + 1;
            Atemp = zeros(imSize);
            Atemp(patches{i}(1):patches{i}(2),patches{i}(3):patches{i}(4)) = ...
                reshape(RESULTS(i).A(:,k),patches{i}(2)-patches{i}(1)+1,patches{i}(4)-patches{i}(3)+1);            
            A(:,cnt) = sparse(double(Atemp(:)));
    end
    
    b_temp = sparse(imSize(1),imSize(2));
    b_temp(patches{i}(1):patches{i}(2),patches{i}(3):patches{i}(4)) = reshape(RESULTS(i).b,patches{i}(2)-patches{i}(1)+1,patches{i}(4)-patches{i}(3)+1);  
    B(:,i) = double(b_temp(:));
    F(i,:) = double(RESULTS(i).f);
    MASK(patches{i}(1):patches{i}(2),patches{i}(3):patches{i}(4)) = MASK(patches{i}(1):patches{i}(2),patches{i}(3):patches{i}(4)) + 1;
%     IND(patches{i}(1):patches{i}(2),patches{i}(3):patches{i}(4)) = IND(patches{i}(1):patches{i}(2),patches{i}(3):patches{i}(4)) + 1;
end

A = spdiags(1./MASK(:),0,prod(imSize),prod(imSize))*A;
B = spdiags(1./MASK(:),0,prod(imSize),prod(imSize))*B;
% A(A<0) = 0;
% B(B<0) = 0;
% C(C<0) = 0;
% F(F<0) = 0;

fprintf(' done. \n');
% clear RESULTS

%% load subset of high temporal resolution data to get imaging noise
minNoisePrctile = 5;
nMovs = 15;
movNums = round(linspace(2,length(acqObj.correctedMovies.slice.channel.fileName)-1,nMovs));
allSN = nan(512^2,nMovs);
for nMov = 1:nMovs
    tempMov = single(readCor(acqObj,movNums(nMov)));
    thisSN = get_noise_fft(tempMov,options);
    allSN(:,nMov) = thisSN(:);    
end
sn = mean(allSN,2);
clear P
P.snRaw = sn;
% Enforce minimum noise
snThresh = prctile(sn,minNoisePrctile);
sn(sn<snThresh) = snThresh;
P.snThresh = sn;
% Adjust for downsampling
sn = sn ./ sqrt(memMap.dsRatio);
P.snDS = sn;
%% Eliminate Empty Sources and Initialize common background
emptyROIs = find(sum(A,1)==0);
if ~isempty(emptyROIs)
        A(:,emptyROIs) = [];
        C(emptyROIs,:) = [];
end

f1 = nanmedian(F);
fSplit(1,:) = medfilt1(f1,1e3,'truncate');
fSplit(2,:) = f1-fSplit(1,:);
fBase = prctile(fSplit(2,:),10);
fSplit(1,:) = fSplit(1,:)+fBase;
fSplit(2,:) = fSplit(2,:)-fBase;
f(1,:) = fSplit(1,:).*linspace(1,0,length(f1));
f(2,:) = fSplit(1,:).*linspace(0,1,length(f1));
f(3,:) = fSplit(2,:);
% bReps = options.nb/2;
% fStd = std(f1);
% f(1:bReps,:) = repmat(fSplit(1,:),bReps,1) + fStd*1/2*randn(bReps,size(fSplit,2));
% f(1+bReps:options.nb,:) = repmat(fSplit(2,:),bReps,1) + fStd*1/2*randn(bReps,size(fSplit,2));

%% Merge sources and eliminate very weak sources
fracRetain = 1;
nSources = size(A,2);
P.b = cell(nSources,1);
P.c1 = cell(nSources,1);
P.gn = cell(nSources,1);
P.neuron_sn = cell(nSources,1);

Km = 0;
Kn = nSources;

while Km < Kn
    Kn = size(A,2);
    [A,C] = merge_components([],A,[],C,f,P,C,options);
    Km = size(A,2)
end

[A,C] = order_ROIs(A,C);

numRetain = ceil(fracRetain*size(A,2));
A = A(:,1:numRetain);
C = C(1:numRetain,:);

%% Add initImages as Sources

if ~isempty(initImages)
    fprintf('Adding initImages to Source Initialization...'),
    % Add initImages to source matrix
    A = cat(2,A,sparse(reshape(initImages,prod(imSize),size(initImages,3))));
    % Normalize spatial components
    A = bsxfun(@rdivide,A,sqrt(sum(A.^2)));
    % use temporary meanRef as background (note, do NOT replace f)
    b = reshape(meanRef(acqObj),prod(imSize),1);
    b = size(A,2).*b./sqrt(sum(b.^2));
    % Ensure components are finite
    A(~isfinite(A)) = 0;
    b(~isfinite(b)) = 0;
    if ~isempty(data)
        [C,~,A] = pinv_temporal_components(data,nSlice,A,b);
    else
        [C,~,A] = pinv_temporal_components(acqObj,nSlice,A,b);
    end
    fprintf(' done. \n');
end

%% Initialize dummy variables
nSources = size(A,2);
P.b = cell(nSources,1);
P.c1 = cell(nSources,1);
P.gn = cell(nSources,1);
P.neuron_sn = cell(nSources,1);

% label initImage sources, hacky solution here
if ~isempty(initImages)
    for i = 1:size(initImages,3)
        sourceInd = i - size(initImages,3) + nSources;
        P.b{sourceInd} = i;
    end
end 

%% Baseline and normalize traces

binSize = 500 / memMap.dsRatio;
basePrct = 2;

for nAcq = 1:size(acqBlocks,1)
    thisBlock = acqBlocks(nAcq,:);
    acqInd = ceil(thisBlock(1)/ memMap.dsRatio):floor(thisBlock(2)/ memMap.dsRatio);
    parfor nSource = 1:size(C,1)
        C(nSource,acqInd) = removeSourceBaseline(C(nSource,acqInd),binSize,basePrct);
    end
end

C(C<0) = 0;
cNorm = sqrt(sum(C.^2,2));
C = bsxfun(@rdivide,C,cNorm);
fNorm = sqrt(sum(f.^2,2));
f = bsxfun(@rdivide,f,fNorm);
scaleFactor = 1e3;
f = f * scaleFactor;

%% First pass to clean up initialization
P.sn = P.snDS;
[A,b,C,f,P,options] = updateCNMF_all...
    (A,C,f,P,options,acqObj,data,acqBlocks,memMap,nSlice);
%% Enforce robustness with noise inflation
noiseTolerance = 1.1;
P.sn = P.snDS * noiseTolerance;
[A,b,C,f,P,options] = updateCNMF_all...
    (A,C,f,P,options,acqObj,data,acqBlocks,memMap,nSlice);

%% Clean up robust results

cNorm = sqrt(sum(C.^2,2));
C = bsxfun(@rdivide,C,cNorm);
fNorm = sqrt(sum(f.^2,2));
f = bsxfun(@rdivide,f,fNorm);
scaleFactor = 1e3;
f = f * scaleFactor;

P.sn = P.snDS;
[A,b,C,f,P,options] = updateCNMF_all...
    (A,C,f,P,options,acqObj,data,acqBlocks,memMap,nSlice);

%% Save Results

saveFile = fullfile(acqObj.defaultDir,sprintf('Slice%0.2d_patchResults_v170306.mat',nSlice));
acqObj.roiInfo.slice(nSlice).NMF.filename = saveFile;
save(saveFile,'A','b','C','f','P','options','optionsThresh'),

end

function [A,b,C,f,P,options] = updateCNMF_all...
    (A,C,f,P,options,acqObj,data,acqBlocks,memMap,nSlice)

fprintf('Updating spatial components...');
pctRunOnAll warning('off','MATLAB:nargchk:deprecated'),
if ~isempty(data) %data in memory
    [A,b,C] = update_spatial_components(reshape(data,size(A,1),size(C,2)),C,f,A,P,options);
else %memMapped
    [A,b,C] = update_spatial_components(memMap,C,f,A,P,options);
end
fprintf(' done. \n');

fprintf('Updating temporal components... (1)'),
if ~isempty(data)
    [C,f] = coDesc_temporal_components(data,nSlice,A,b,C,f);
else
    [C,f] = coDesc_temporal_components(acqObj,nSlice,A,b,C,f);
end
fprintf(' done. \n'),

% % Re-normalize spatial components
% A = bsxfun(@rdivide,A,sqrt(sum(A.^2)));
% % scale b as well? maybe not, since it is supposed to be much bigger
% % compromise might be scale b so its larger by factor of nSources 
% b = bsxfun(@rdivide,b,sqrt(sum(b.^2)))*size(A,2);

% fprintf('Updating temporal components... (1)')
% if ~isempty(data)
%     [C,f,A] = pinv_temporal_components(data,nSlice,A,b);
% else
%     [C,f,A] = pinv_temporal_components(acqObj,nSlice,A,b);
% end
% fprintf(' done. \n');

fprintf('Merging overlaping components... (1)')
Km = 0;
Kn = size(A,2);
while Km < Kn
    Kn = size(A,2);
    [A,C,~,~,P] = merge_components([],A,b,C,f,P,C,options);
    Km = size(A,2),
end
fprintf(' done. \n');

end

% Remove Source Baseline and neuropil (rough estimate w robust fit)
% This is just to help background + neuropil signals be absorbed into the
% background component by taking (most of) it out of the neural sources
% warnState = warning('off', 'stats:statrobustfit:IterationLimit');
% for nAcq = 1:size(acqBlocks,1)
%     acqInd = acqBlocks(nAcq,1):acqBlocks(nAcq,2);
%     fSub = removeSourceBaseline(sum(f(:,acqInd),1));
%     parfor nSource = 1:size(C,1)
%         warning('off', 'stats:statrobustfit:IterationLimit');
%         cSub = removeSourceBaseline(C(nSource,acqInd));
%         pilFit = robustfit(fSub,cSub,'bisquare',2);
%         C(nSource,acqInd) = cSub - pilFit(1) - pilFit(2)*fSub;
%     end
% end
% warning(warnState);
