function result = patchInitNMF(acqObj,nSlice,patches,patchNum,nFactors,tBin)

%% Setup and Data
result = struct();
imSize = acqObj.correctedMovies.slice(nSlice).channel.size(1,1:2);
nFrames = sum(acqObj.correctedMovies.slice(nSlice).channel.size(:,3));

dMap = memmapfile(acqObj.indexedMovie.slice(nSlice).channel(1).fileName,...
    'Format', {'int16', [nFrames, prod(imSize)], 'mov'});
mov = dMap.data.mov;

% Get indices for current patch and convert to column-major format of binary file
patchMask = zeros(imSize);
patchMask(patches{patchNum}(1):patches{patchNum}(2),patches{patchNum}(3):patches{patchNum}(4)) = 1;
[matRow, matCol] = ind2sub(imSize, find(patchMask(:)));
binPixInd = sub2ind([imSize(2), imSize(1)], matCol, matRow);
Y = mov(:,binPixInd);
Y = permute(Y,[2 1]);
clear mov
clear dMap

[P,Y] = preprocess_data(single(Y));

zInds = 1:nFrames-mod(nFrames,tBin);
Ybin = squeeze(mean(reshape(Y(:,zInds),length(binPixInd),tBin,floor(nFrames/tBin)),2));

%% Extract Factors and eliminate redundandant sources
[w,t] = NMF_SNC_Factors(Ybin,nFactors);
b = w(:,nFactors+1);
Ysub = Ybin - (b * t(nFactors+1,:));
for s = 1:nFactors
    sProj(s,:) = w(:,s)'*Ysub;
end
corrGraph = graph(corrcoef(sProj')>0.9);
conComp = conncomp(corrGraph,'OutputForm','cell');
nSources = length(conComp);
A=nan(size(w,1),length(conComp));
for nComp = 1:nSources
    compSources = conComp{nComp};
    A(:,nComp) = mean(w(:,compSources),2);
end
%% Extract high res timeseries and format output
t = pinv([A,b])*Y;

tNorm = sqrt(sum(t.^2,2));
t = bsxfun(@rdivide, t, tNorm);
A = bsxfun(@times, A, tNorm(1:nSources)');
b = b * tNorm(nSources+1);

%% Output in results structure format
result.conComp = conComp;
result.A = A;
result.b = b;
result.C = t(1:nSources,:);
result.f = t(nSources+1,:);
result.P = P;