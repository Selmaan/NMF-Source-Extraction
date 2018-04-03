function [Cf,AY] = update_temporal_components_fromTiff(acqObj,nSlice)

if ~exist('nSlice','var') || isempty(nSlice)
    nSlice = 1;
end

%% Load components and precalculate pinv
load(acqObj.roiInfo.slice(nSlice).NMF.filename),
Ab = full([A,b]);
pA = pinv(Ab);

%% 
movSizes = acqObj.correctedMovies.slice(nSlice).channel(1).size;
imSize = movSizes(1,1:2);
ref = reshape(meanRef(acqObj),prod(imSize),1);
AY = nan(size(pA,1),sum(movSizes(:,3)),'single');
Cf = nan(size(pA,1),sum(movSizes(:,3)),'single');
for nMov = 1:size(movSizes,1)
    fprintf('Processing Movie %d of %d ... \n',nMov,size(movSizes,1)),
    frameOffset = sum(movSizes(1:nMov-1,3));
    frameInd = (1:movSizes(nMov,3))+frameOffset;
    tmpMov = readCor(acqObj,nMov,'single',nSlice);
    tmpMov = reshape(tmpMov,prod(imSize),[]);
    if sum(~isfinite(tmpMov(:)))
        replaceFrames = find(sum(~isfinite(tmpMov),1)>0);
        for nFrame = replaceFrames
            replacePix = ~isfinite(tmpMov(:,nFrame));
            tmpMov(replacePix,nFrame) = ref(replacePix);
        end
    end
    AY(:,frameInd) = Ab' * tmpMov;
    Cf(:,frameInd) = pA * tmpMov;
end

%% Coordinate-descent with non-negative HALS
nSources = size(A,2);
nBkgd = size(b,2);
C = Cf(1:nSources,:);
C(C<0) = 0;
Cf(1:nSources,:) = C;
clear C
AA = full([A,b]'*[A,b]);
aa = diag(AA);

maxIter = 20;
repeat = 1;
thisIter = 0;
while repeat
    thisIter = thisIter+1;
    Cf_ = Cf;
    for k=1:(nSources+nBkgd)
        ck = Cf(k,:) + (AY(k,:)-AA(k,:)*Cf)/aa(k);
        if k <= nSources
            Cf(k,:) = max(0,ck);
        else
            Cf(k,:) = ck;
        end
    end
    
    repeat = (thisIter<maxIter) && norm(Cf-Cf_,'fro')/norm(Cf_,'fro') > 1e-4;
end

%% Save Traces (and AY)
saveDir = acqObj.defaultDir;
fName = sprintf('Slice%0.2d_fullResTraces_v180318',nSlice);
save(fullfile(saveDir,fName),'Cf','AY'),
acqObj.roiInfo.slice(nSlice).NMF.traceFn = fullfile(saveDir,fName);
