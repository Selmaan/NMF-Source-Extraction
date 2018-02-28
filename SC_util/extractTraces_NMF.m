function [dF,deconv,denoised,Gs,Lams,A,b,f] = extractTraces_NMF(acqObj,cellLabels,interpFrames)

if nargin < 2 || isempty(cellLabels)
    cellLabels = cell(length(acqObj.roiInfo.slice),1);
end

if nargin < 3
    interpFrames = [];
end

% Get frame rate and conversion factor for volume imaging
if isfield(acqObj.metaDataSI,'SI')
    metaData = acqObj.metaDataSI.SI;
else
    metaData = acqObj.metaDataSI;
end

frameRate = metaData.hRoiManager.scanVolumeRate;
gDecay = exp(-1/frameRate); %initialize deconv at 1s tau value
fprintf('Deconvolution initialized at tau-1s, g is %0.3f \n',gDecay);

if metaData.hFastZ.enable
    frame2slice = 1/metaData.hFastZ.numFramesPerVolume;
else
    frame2slice = 1/1;
end

%% Get slice and acq info and load data
nSlices = length(acqObj.roiInfo.slice);
nBlocks = size(acqObj.syncInfo.sliceFrames,1);

A = cell(1,nSlices);
b = cell(1,nSlices);
C = cell(nSlices,1);
f = cell(nSlices,1);

% Load NMF sources and concatenate data
for nSlice = 1:nSlices
    D = load(acqObj.roiInfo.slice(nSlice).NMF.filename);
    load(acqObj.roiInfo.slice(nSlice).NMF.traceFn,'Cf');
    A{nSlice} = D.A;
    b{nSlice} = D.b;
    nSources = size(A{nSlice},2);
    nBkgd = size(b{nSlice},2);
    if ~isempty(cellLabels{nSlice})
        validSources = cellLabels{nSlice};
        A{nSlice} = D.A(:,validSources);
    else
        validSources = 1:nSources;
    end
    thisC = cell(1,nBlocks);
    thisF = cell(1,nBlocks);
    for nBlock = 1:nBlocks
        % Clip 'extra' frames from slices without complete volume
        if nBlock > 1
            offsetInd = acqObj.syncInfo.sliceFrames(nBlock-1,nSlice);
        else
            offsetInd = 0;
        end
        blockLength = acqObj.syncInfo.validFrameCount(nBlock) * frame2slice;
        thisC{nBlock} = Cf(validSources,offsetInd+1:offsetInd+blockLength);
        thisF{nBlock} =  Cf(nSources+1:nSources+nBkgd,offsetInd+1:offsetInd+blockLength);
%         C{nSlice,nBlock} = Cf(1:nSources,offsetInd+1:offsetInd+blockLength);
%         f{nSlice,nBlock} = Cf(nSources+1:nSources+nBkgd,offsetInd+1:offsetInd+blockLength);
    end
    C{nSlice} = cell2mat(thisC);
    f{nSlice} = cell2mat(thisF);
end

clear Cf D
%% Get Traces
dF = cell(nSlices,1);
deconv = cell(nSlices,1);
denoised = cell(nSlices,1);
Gs = cell(nSlices,1);
Lams = cell(nSlices,1);

for nSlice = 1:nSlices
    % get F_baseline from background component
    thisA = A{nSlice};
    thisB = b{nSlice};
    thisF = f{nSlice};
    if size(thisB,2)==3 && size(thisF,1)==3
        % Use only 'tonic' components of baseline, ignoring phasic/neuropil
        thisB = thisB(:,1:2);
        thisF = thisF(1:2,:);
    else
        warning('Non-Standard Baseline Used'),
    end
    sumA = full(sum(thisA));
    normA = bsxfun(@rdivide,thisA,sumA);
    baseF = normA'*thisB * median(thisF,2);
    thisC = bsxfun(@times,C{nSlice},full(sum(normA.*thisA))');
%     thisC = bsxfun(@times,C{nSlice},sumA');
%     thisC = bsxfun(@times,C{nSlice},full(diag(normA'*thisA)));
    
    % Interpolate Data btw Blank Frames
    if ~isempty(interpFrames)
        nInterp = length(interpFrames)-1;
        blankFrames = interpFrames{1};
        interpVals = 0*thisC(:,blankFrames);
        interpF = 0*thisF(:,blankFrames);
        for i=1:nInterp
            interpVals = interpVals + thisC(:,interpFrames{i+1})/nInterp;
            interpF = interpF + thisF(:,interpFrames{i+1})/nInterp;
        end
        thisC(:,blankFrames) = interpVals;
        thisF(:,blankFrames) = interpF;
        f{nSlice} = thisF;
    end

    % Deconvolve all Signals
    thisDenoised = nan(size(thisC));
    thisDeconv = nan(size(thisC));
    thisB = nan(size(thisC,1),1);
    thisG = nan(size(thisC,1),1);
    thisLam = nan(size(thisC,1),1);
    parfor_progress(size(thisC,1));
    parfor nSig = 1:size(thisC,1)
        parfor_progress;
        [c, s, b, g, lam] = sc_constrained_oasisAR1(double(thisC(nSig,:)), gDecay);
        thisDenoised(nSig,:) = b+c;
        thisDeconv(nSig,:) = s;
        thisB(nSig) = b;
        thisG(nSig) = g;
        thisLam(nSig) = lam;
    end
    parfor_progress(0);
    
    % Use inferred baseline and background to get dF/F, then scale trace,
    % denoised and spiking data
    thisB(thisB<0) = 0;
    baseF(baseF<0) = 0;
    F_ = double(baseF) + thisB;
    thisDF = bsxfun(@rdivide,double(thisC),F_);
    thisDenoised = bsxfun(@rdivide,thisDenoised,F_);
    thisDeconv = bsxfun(@rdivide,thisDeconv,F_);
    
    % Store data for each slice
    dF{nSlice} = thisDF;
    deconv{nSlice} = thisDeconv;
    denoised{nSlice} = thisDenoised;
    Gs{nSlice} = thisG;
    Lams{nSlice} = thisLam;
end
