function [C,f] = coDesc_temporal_components(acqObj_data,nSlice,A,b,C,f)

Cf = [C;f];
clear C f,

%% 
if ~isnumeric(acqObj_data)
    memMap = matfile(acqObj_data.indexedMovie.slice(nSlice).channel(1).memMap);
    imSize = acqObj_data.correctedMovies.slice(nSlice).channel(1).size(1,1:2);
    nFrames = size(memMap,'Y',3);
    ref = reshape(meanRef(acqObj_data),prod(imSize),1);
    useMemMap = true;
else
    imSize = [size(acqObj_data,1),size(acqObj_data,2)];
    nFrames = size(acqObj_data,3);
    ref = reshape(nanmean(acqObj_data,3),prod(imSize),1);
    useMemMap = false;
end
ref(~isfinite(ref)) = nanmean(ref);
step = 1000;
frameBatches = 1:step:nFrames;
frameBatches(2,:) = min(frameBatches(1,:) + step - 1, nFrames);

AY = nan(size(A,2)+size(b,2), nFrames);
Ab = full([A,b]);
for frameBatch = 1:size(frameBatches,2)
    fInd = frameBatches(1,frameBatch):frameBatches(2,frameBatch);
    if useMemMap
        tempY = memMap.Yr(:,fInd);
    else
        tempY = reshape(acqObj_data(:,:,fInd),prod(imSize),length(fInd));
    end
    if sum(~isfinite(tempY(:)))
        replaceFrames = find(sum(~isfinite(tempY),1)>0);
        for nFrame = replaceFrames
            replacePix = ~isfinite(tempY(:,nFrame));
            tempY(replacePix,nFrame) = ref(replacePix);
        end
    end
    
    AY(:,fInd) = Ab' * tempY;
end


%% Do non-negative coordinate descent on source traces
nSources = size(A,2);
nBkgd = size(b,2);
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

C = Cf(1:nSources,:);
f = Cf(nSources+1:nSources+nBkgd,:);
