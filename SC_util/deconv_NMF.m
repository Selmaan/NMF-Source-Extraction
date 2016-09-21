function [dF,dF_denoised,dF_deconv,...
    traceBs,traceGs,traceSNs,A,b,f] = deconv_NMF(acqObj)

%% Get slice and acq info and load data
nSlices = length(acqObj.roiInfo.slice);
nBlocks = size(acqObj.syncInfo.sliceFrames,1);

A = cell(1,nSlices);
b = cell(1,nSlices);
C = cell(nSlices,nBlocks);
f = cell(nSlices,nBlocks);

% Load NMF sources and concatenate data
for nSlice = 1:nSlices
    D = load(acqObj.roiInfo.slice(nSlice).NMF.filename);
    A{nSlice} = D.A;
    b{nSlice} = D.b;
    for nBlock = 1:nBlocks
        % Clip 'extra' frames from slices without complete volume
        if nBlock > 1
            offsetInd = acqObj.syncInfo.sliceFrames(nBlock-1,nSlice);
        else
            offsetInd = 0;
        end
        blockLength = acqObj.syncInfo.validFrameCount(nBlock)/5;
        C{nSlice,nBlock} = D.C(:,offsetInd+1:offsetInd+blockLength);
        f{nSlice,nBlock} = D.f(:,offsetInd+1:offsetInd+blockLength);
    end
end

%% Deconvolve all traces

dF = cell(nSlices,1);
dF_denoised = cell(nSlices,1);
dF_deconv = cell(nSlices,1);
traceBs = cell(nSlices,1);
traceGs = cell(nSlices,1);
traceSNs = cell(nSlices,1);

for nSlice = 1:nSlices
    % get F_baseline from 1st (tonic) background component
    thisA = A{nSlice};
    thisB = b{nSlice};
    thisF = cell2mat(f(nSlice,:));
    normA = bsxfun(@rdivide,thisA,sum(thisA));
    baseF = median(thisF(1,:)) * normA'*thisB(:,1);
    
    % remove trace baseline for each acquisition block independently
    for nBlock = 1:nBlocks
        thisC = C{nSlice,nBlock};
        parfor nSig = 1:size(thisC,1)
            thisC(nSig,:) = removeSourceBaseline(thisC(nSig,:)/baseF(nSig));
        end
        C{nSlice,nBlock} = thisC;
    end
    
    % Pre-allocate memory for all data from this slice
    thisDF = cell2mat(C(nSlice,:));
    thisDF_denoised = nan*thisDF;
    thisDF_deconv = nan*thisDF;
    thisBs = nan*thisDF(:,1);
    thisGs = nan*thisDF(:,1);
    thisSNs = nan*thisDF(:,1);
    
    % Deconvolve each source
    parfor nSig = 1:size(thisDF,1)
        [cDe,bs,c1,g,sn,sp] = constrained_foopsi(thisDF(nSig,:));
        thisDF_denoised(nSig,:) = cDe + bs;
        thisDF_deconv(nSig,:) = sp;
        thisBs(nSig) = bs;
        thisGs(nSig) = g(1);
        thisSNs(nSig) = sn(1);
    end
    
    % Store data for each slice
    dF{nSlice} = thisDF;
    dF_denoised{nSlice} = thisDF_denoised;
    dF_deconv{nSlice} = thisDF_deconv;
    traceBs{nSlice} = thisBs;
    traceGs{nSlice} = thisGs;
    traceSNs{nSlice} = thisSNs;
end
