function [dF_cells,dF_denoised,dF_deconv,...
    traceBs,traceGs,traceSNs,traceSnScales,A,b,f,cIds] = deconvCells_NMF(acqObj,cellLabel,interpFrames)

if nargin < 2 || isempty(cellLabel)
    cellLabel = 1;
end

if nargin < 3 || isempty(interpFrames)
    interpFrames = [];
end

[dF,A,b,f] = extractTraces_NMF(acqObj);
nSlices = length(dF);

if ~isempty(interpFrames)
    for nSlice = 1:nSlices
        thisDF = dF{nSlice};
        thisF = f{nSlice};
        nInterp = length(interpFrames)-1;
        blankFrames = interpFrames{1};
        interpVals = 0*thisDF(:,blankFrames);
        interpF = 0*thisF(:,blankFrames);
        for i=1:nInterp
            interpVals = interpVals + thisDF(:,interpFrames{i+1})/nInterp;
            interpF = interpF + thisF(:,interpFrames{i+1})/nInterp;
        end
        thisDF(:,blankFrames) = interpVals;
        thisF(:,blankFrames) = interpF;
        dF{nSlices} = thisDF;
        f{nSlice} = thisF;
    end
end

fR = acqObj.metaDataSI.SI.hRoiManager.scanVolumeRate;
fprintf('Frame Rate: %.2f \n',fR);

cIds = cell(nSlices,1);
dF_cells = cell(nSlices,1);
dF_denoised = cell(nSlices,1);
dF_deconv = cell(nSlices,1);
for nSlice = 1:nSlices
    thisA = A{nSlice};
    sourceLabels = clusterSourcesWithCurrentNn(thisA);
    cId = find(sourceLabels==cellLabel);
    fprintf('%d Cells Found in Slice %d \n',length(cId),nSlice);
    
    thisDF = dF{nSlice}(cId,:);
    thisDF_denoised = nan*thisDF;
    thisDF_deconv = nan*thisDF;
    thisBs = nan*thisDF(:,1);
    thisGs = nan*thisDF(:,1);
    thisSNs = nan*thisDF(:,1);
    thisSnScale = nan*thisDF(:,1);
    
    % Deconvolve each source
    fprintf('Solving Deconvolution for Slice %0.2d \n',nSlice),
    parfor_progress(size(thisDF,1));
    parfor nSig = 1:size(thisDF,1)
        parfor_progress;
        try
            [cDe,bs,c1,g,sn,sp,snScale] = constrained_foopsi...
                (thisDF(nSig,:),[],[],[],[],[],fR);
            thisDF_denoised(nSig,:) = cDe + bs;
            thisDF_deconv(nSig,:) = sp;
            thisBs(nSig) = bs;
            thisGs(nSig) = g(1);
            thisSNs(nSig) = sn;
            thisSnScale(nSig) = snScale;
        catch
            warning('Deconvolution Failure in Slice %d Sig %d',nSlice,nSig),
        end
    end
    parfor_progress(0);
    
    % Store data for each slice
    cIds{nSlice} = cId;
    dF_cells{nSlice} = thisDF;
    dF_denoised{nSlice} = thisDF_denoised;
    dF_deconv{nSlice} = thisDF_deconv;
    traceBs{nSlice} = thisBs;
    traceGs{nSlice} = thisGs;
    traceSNs{nSlice} = thisSNs;
    traceSnScales{nSlice} = thisSnScale;
end