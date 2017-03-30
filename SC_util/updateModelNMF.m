function [A,b,C,f,P] = updateModelNMF(acqObj,nSlice,A,b,C,f,P)

syncObj = acqObj.syncInfo;
acqBlocks = [1 syncObj.sliceFrames(1,nSlice)];
for blockNum = 2:size(syncObj.sliceFrames,1)
    acqBlocks(blockNum,:) = ...
        [1+syncObj.sliceFrames(blockNum-1,nSlice), syncObj.sliceFrames(blockNum,nSlice)];
end

imSize = acqObj.correctedMovies.slice(nSlice).channel(1).size(1,1:2);
nFrames = sum(acqObj.correctedMovies.slice(nSlice).channel(1).size(:,3));

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

%%
warnState = warning('off', 'stats:statrobustfit:IterationLimit');
for nAcq = 1:size(acqBlocks,1)
    acqInd = acqBlocks(nAcq,1):acqBlocks(nAcq,2);
    parfor nSource = 1:size(C,1)
        warning('off', 'stats:statrobustfit:IterationLimit');
        C(nSource,acqInd) = removeSourceBaseline_lowpassfilter(C(nSource,acqInd));
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