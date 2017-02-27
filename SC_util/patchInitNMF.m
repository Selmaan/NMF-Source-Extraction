function result = patchInitNMF(acqObj,nSlice,patches,patchNum,nFactors)

%% Setup and Data
result = struct();

memMap = matfile(acqObj.indexedMovie.slice(nSlice).channel(1).memMap);
Y = memMap.Y(patches{patchNum}(1):patches{patchNum}(2),patches{patchNum}(3):patches{patchNum}(4),:);
Y = reshape(Y,size(Y,1)*size(Y,2),size(Y,3));

%% Extract Factors and eliminate redundandant sources
[w,t] = NMF_SNC_Factors(Y,nFactors);
b = w(:,nFactors+1);
Ysub = Y - (b * t(nFactors+1,:));
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
result.P = [];