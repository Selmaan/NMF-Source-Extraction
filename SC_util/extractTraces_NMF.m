function [dF,A,b,f] = extractTraces_NMF(acqObj)

% Get frame rate and conversion factor for volume imaging
if acqObj.metaDataSI.SI.hFastZ.enable
    frame2slice = 1/acqObj.metaDataSI.SI.hFastZ.numFramesPerVolume;
else
    frame2slice = 1/1;
end

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
        blockLength = acqObj.syncInfo.validFrameCount(nBlock) * frame2slice;
        C{nSlice,nBlock} = D.C(:,offsetInd+1:offsetInd+blockLength);
        f{nSlice,nBlock} = D.f(:,offsetInd+1:offsetInd+blockLength);
    end
end

%% Get Traces
dF = cell(nSlices,1);
for nSlice = 1:nSlices
    % get F_baseline from 1st (tonic) background component
    thisA = A{nSlice};
    thisB = b{nSlice};
    thisF = cell2mat(f(nSlice,:));
    normA = bsxfun(@rdivide,thisA,sum(thisA));
    baseF = normA'*thisB * median(thisF,2);
%     baseF = median(thisF(1,:),2) * normA'*thisB(:,1);
    

    % remove trace baseline for each acquisition block independently
    for nBlock = 1:nBlocks
        thisC = C{nSlice,nBlock};
        for nSig = 1:size(thisC,1)
            thisC(nSig,:) = removeSourceBaseline(thisC(nSig,:)/baseF(nSig));
        end
        C{nSlice,nBlock} = thisC;
    end
    
    thisDF = cell2mat(C(nSlice,:));
    
    % Store data for each slice
    dF{nSlice} = thisDF;
end
