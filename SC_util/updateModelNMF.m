function [A,b,C,f,P] = updateModelNMF(acqObj,nSlice,A,b,C,f,P)

% load(FOV1.roiInfo.slice(nSlice).NMF.filename),

% Global Time-Constant
gTau = 1.5;

if acqObj.metaDataSI.SI.hFastZ.enable
    frameRate = round(acqObj.metaDataSI.SI.hRoiManager.scanFrameRate...
        /acqObj.metaDataSI.SI.hFastZ.numFramesPerVolume);
else
    frameRate = round(acqObj.metaDataSI.SI.hRoiManager.scanFrameRate);
end

%global gValue
gVal = 1 - 1/(frameRate*gTau);

syncObj = acqObj.syncInfo;
acqBlocks = [1 syncObj.sliceFrames(1,nSlice)];
for blockNum = 2:size(syncObj.sliceFrames,1)
    acqBlocks(blockNum,:) = ...
        [1+syncObj.sliceFrames(blockNum-1,nSlice), syncObj.sliceFrames(blockNum,nSlice)];
end

imSize = acqObj.correctedMovies.slice(nSlice).channel.size(1,1:2);
nFrames = sum(acqObj.correctedMovies.slice(nSlice).channel.size(:,3));

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

%%
fprintf('\n Projecting Traces with Global Tau: %0.2d',(1/(1-gVal))/frameRate),
warnState = warning('off', 'stats:statrobustfit:IterationLimit');
parfor_progress(size(acqBlocks,1)*size(C,1));
for nAcq = 1:size(acqBlocks,1)
    acqInd = acqBlocks(nAcq,1):acqBlocks(nAcq,2);
%     gMat = make_G_matrix(length(acqInd),gVal);
    parfor nSource = 1:size(C,1)
        parfor_progress;
        warning('off', 'stats:statrobustfit:IterationLimit');
        C(nSource,acqInd) = removeSourceBaseline(C(nSource,acqInd));
        C(nSource,acqInd) = constrained_foopsi(C(nSource,acqInd),[],[],gVal);
    end
end
warning(warnState);
parfor_progress(0);
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

% Trying alternate change here, where f is scaled way up compared to C, so
% that weights on background signals are effectively free
scaleFactor = 1e3;
f = f * scaleFactor;

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