function extractSourcesNMF(acqObj,nSlice)

syncObj = acqObj.derivedData(1).syncObj;
acqBlocks = [1 syncObj.sliceFrames(1,nSlice)];
for blockNum = 2:size(syncObj.sliceFrames,1)
    acqBlocks(blockNum,:) = ...
        [1+syncObj.sliceFrames(blockNum-1,nSlice), syncObj.sliceFrames(blockNum,nSlice)];
end

imSize = acqObj.correctedMovies.slice(nSlice).channel.size(1,1:2);
nFrames = sum(acqObj.correctedMovies.slice(nSlice).channel.size(:,3));

if acqObj.metaDataSI.SI.hFastZ.enable
    frameRate = round(acqObj.metaDataSI.SI.hRoiManager.scanFrameRate...
        /acqObj.metaDataSI.SI.hFastZ.numFramesPerVolume);
else
    frameRate = round(acqObj.metaDataSI.SI.hRoiManager.scanFrameRate);
end

%% Set parameters
patch_size = [52,52];                   % size of each patch along each dimension (optional, default: [32,32])
overlap = [6,6];                        % amount of overlap in each dimension (optional, default: [4,4])
nFactors = 15;
tBin = 1 * frameRate;
p = 1;                                            % order of autoregressive system (p = 0 no dynamics, p=1 just decay, p = 2, both rise and decay)
merge_thr = 0.8;                                  % merging threshold

options = CNMFSetParms(...
    'd1',imSize(1),'d2',imSize(2),...
    'search_method','dilate','se',strel('disk',2,0),...
    'deconv_method','constrained_foopsi',...    % activity deconvolution method
    'merge_thr',merge_thr,...                    % merging threshold
    'nB',2, ...
    'medw',[1 1]);

%% Get Patches Results
patches = construct_patches(imSize,patch_size,overlap);

parfor_progress(length(patches));
parfor_progress;
RESULTS = patchInitNMF(acqObj,nSlice,patches,1,nFactors,tBin);
parfor patchNum = 2:length(patches)
    parfor_progress;
    RESULTS(patchNum) = patchInitNMF(acqObj,nSlice,patches,patchNum,nFactors,tBin);
end
parfor_progress(0);

%% combine results into one structure
fprintf('Combining results from different patches...');
C = double(cell2mat({RESULTS(:).C}'));
d = prod(imSize);
A = sparse(d,size(C,1));
P.sn = zeros(imSize);
P.active_pixels = zeros(imSize);
IND = zeros(imSize);
P.b = {};
P.c1 = {};
P.gn = {};
P.neuron_sn = {};

cnt = 0;
B = sparse(prod(imSize),length(patches));
MASK = zeros(imSize);
F = zeros(length(patches),nFrames);
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
    P.sn(patches{i}(1):patches{i}(2),patches{i}(3):patches{i}(4)) = reshape(RESULTS(i).P.sn,patches{i}(2)-patches{i}(1)+1,patches{i}(4)-patches{i}(3)+1);
    IND(patches{i}(1):patches{i}(2),patches{i}(3):patches{i}(4)) = IND(patches{i}(1):patches{i}(2),patches{i}(3):patches{i}(4)) + 1;
end

A = spdiags(1./MASK(:),0,prod(imSize),prod(imSize))*A;
B = spdiags(1./MASK(:),0,prod(imSize),prod(imSize))*B;
% A(A<0) = 0;
% B(B<0) = 0;
% C(C<0) = 0;
% F(F<0) = 0;

fprintf(' done. \n');
clear RESULTS

%% Eliminate Empty Sources and Initialize common background
emptyROIs = find(prod(A==0));
if ~isempty(emptyROIs)
        A(:,emptyROIs) = [];
        C(emptyROIs,:) = [];
end

f = nanmedian(F);
fSplit(1,:) = medfilt1(f,1e3,'truncate');
fSplit(2,:) = f-fSplit(1,:);
fBase = prctile(fSplit(2,:),10);
fSplit(1,:) = fSplit(1,:)+fBase;
fSplit(2,:) = fSplit(2,:)-fBase;
f = fSplit;

%% Merge sources and eliminate very weak sources
fracRetain = 0.95;
nSources = size(A,2);
P.b = cell(nSources,1);
P.c1 = cell(nSources,1);
P.gn = cell(nSources,1);
P.neuron_sn = cell(nSources,1);

Km = 0;
Kn = nSources;

while Km < Kn
    Kn = size(A,2);
    [A,C] = merge_components([],A,[],C,f,P,[],options);
    Km = size(A,2),
end

[A,C] = order_ROIs(A,C);

numRetain = ceil(fracRetain*size(A,2));
A = A(:,1:numRetain);
C = C(1:numRetain,:);
P.b = cell(numRetain,1);
P.c1 = cell(numRetain,1);
P.gn = cell(numRetain,1);
P.neuron_sn = cell(numRetain,1);

%% 1-pass on data to get traces w/out dynamics
% P.p = 0;
% options.temporal_iter = 10;
% [C,f,P] = update_temporal_components(data,A,b,C,f,P,options);


% Remove Source Baseline and neuropil (rough estimate w robust fit)
% This is just to help background + neuropil signals be absorbed into the
% background component by taking (most of) it out of the neural sources
warnState = warning('off', 'stats:statrobustfit:IterationLimit');
for nAcq = 1:size(acqBlocks,1)
    acqInd = acqBlocks(nAcq,1):acqBlocks(nAcq,2);
    fSub = removeSourceBaseline(sum(f(:,acqInd),1));
    parfor nSource = 1:size(C,1)
        warning('off', 'stats:statrobustfit:IterationLimit');
        cSub = removeSourceBaseline(C(nSource,acqInd));
        pilFit = robustfit(fSub,cSub,'bisquare',2);
        C(nSource,acqInd) = ...
            removeSourceBaseline(C(nSource,acqInd) - pilFit(1) - pilFit(2)*fSub);
    end
end
warning(warnState);
 
% Set negative values of C to 0 and rescale C + f
% Rescale because the original lars problem for each pixel minimizes the
% total weight over all sources including background! Don't want to
% penalize the background weight dramatically. This 'shrunk' estimate for
% the background is replaced by OLS at end of algorithm
C(C<0) = 0;
cNorm = sqrt(sum(C.^2,2));
C = bsxfun(@rdivide,C,cNorm);
fNorm = sqrt(sum(f.^2,2));
f = bsxfun(@rdivide,f,fNorm);

fprintf('Updating spatial components... (1)');
warning('off','MATLAB:nargchk:deprecated'),
binFile = acqObj.indexedMovie.slice(nSlice).channel(1).fileName;
[A,b,C] = update_spatial_components(binFile,C,f,A,P,options);
fprintf(' done. \n');

% Re-normalize spatial components
A = bsxfun(@rdivide,A,sqrt(sum(A.^2)));
% scale b as well? maybe not, since it is supposed to be much bigger
% compromise might be scale b so its larger by factor of nSources 
b = bsxfun(@rdivide,b,sqrt(sum(b.^2)))*size(A,2);

fprintf('Updating temporal components... (1)')
[C,f,A] = pinv_temporal_components(acqObj,nSlice,A,b);
fprintf(' done. \n');

fprintf('Merging overlaping components... (1)')
Km = 0;
Kn = size(A,2);
while Km < Kn
    Kn = size(A,2);
    [A,C,~,~,P] = merge_components([],A,b,C,f,P,[],options);
    Km = size(A,2),
end
fprintf(' done. \n');

% if size(YrA,1) > size(C,1) %If we merged components in line above
%     [C,f,P] = update_temporal_components(data,A,b,C,f,P,options);
% end

%% Save Results

saveFile = fullfile(acqObj.defaultDir,sprintf('Slice%0.2d_patchResults_v0912',nSlice));
save(saveFile,'A','b','C','f','P'),
% if ~isempty(saveFile)
%     save(saveFile,'A','b','C','f','P'),
% end

%% update spatial and temporal components

% P.p = 0; % Don't model temporal dynamics on first pass
% options.temporal_iter = 1; %Don't use multiple iterations on first pass
% fprintf('Updating spatial components... (1)');
% warning('off','MATLAB:nargchk:deprecated'),
% [A,b,C] = update_spatial_components(data,C,f,A,P,options);
% fprintf(' done. \n');
% fprintf('Updating temporal components... (1)')
% [C,f,P,S,YrA] = update_temporal_components(data,A,b,C,f,P,options);
% fprintf(' done. \n');
% 
% %% Merge Components
% fprintf('Merging overlaping components... (1)')
% Km = 0;
% Kn = size(A,2);
% 
% while Km < Kn
%     Kn = size(A,2);
%     [A,C,~,~,P,S] = merge_components([],A,b,C,f,P,S,options);
%     Km = size(A,2),
% end
% fprintf(' done. \n');
% 
% %% Update components again and merge
% P.p = P.pFinal; % P.p was set to 0 for first spatio-temporal update
% fprintf('Updating Spatial components (2)... ')
% [A,b,C] = update_spatial_components(data,C,f,A,P,options);
% fprintf(' done. \n');
% fprintf('Updating temporal components (2)... ')
% [C,f,P,S,YrA] = update_temporal_components(data,A,b,C,f,P,options);
% fprintf(' done. \n');
% fprintf('Merging overlaping components... (2)')
% Km = 0;
% Kn = size(A,2);
% while Km < Kn
%     Kn = size(A,2);
%     [A,C,~,~,P,S] = merge_components([],A,b,C,f,P,S,options);
%     Km = size(A,2),
% end
% fprintf(' done. \n');
% %% Final spatio-temporal update
% fprintf('Updating spatial components (3)...');
% options.temporal_iter = options.temporal_iterFinal; %Use multiple iterations for final pass
% [A,b,C] = update_spatial_components(data,C,f,A,P,options);
% fprintf(' done. \n');
% fprintf('Updating temporal components (3)... ')
% [C,f,P,S,YrA] = update_temporal_components(data,A,b,C,f,P,options);
% fprintf(' done. \n');

% %% Save Results
% if ~isempty(saveFile)
%     save(saveFile,'A','b','C','f','P','S','YrA'),
% end