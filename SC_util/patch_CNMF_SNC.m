function result = patch_CNMF_SNC(data,patches,p,patchNum,nFactors,tBin)

if ~exist('nFactors','var') || isempty(nFactors)
    nFactors = 10;
end

if ~exist('tBin','var') || isempty(tBin)
    tBin = 6;   
end
%% Setup and Parameters
result = struct();
sizY = data.sizY;
rInds = patches{patchNum}(1):patches{patchNum}(2);
cInds = patches{patchNum}(3):patches{patchNum}(4);
Y = nan(length(rInds),length(cInds),sizY(3));
for i = 1:length(cInds)
    cInd = cInds(i);
    thisOffset = (cInd-1)*sizY(1);
    thisInds = rInds+thisOffset;
    Y(:,i,:) = reshape(data.Yr(thisInds,:),length(rInds),1,sizY(3));
end
Y(isnan(Y)) = 0;

[d1,d2,T] = size(Y);
d3 = 1;
d = d1*d2*d3;
[P,Y] = preprocess_data(double(Y),p);

zInds = 1:data.sizY(1,3)-mod(data.sizY(1,3),tBin);
Yvec = squeeze(mean(reshape(Y(:,:,zInds),size(Y,1),size(Y,2),tBin,[]),3));
Yvec = reshape(Yvec,size(Y,1)*size(Y,2),[]);

%% Extract Factors and eliminate redundandant sources
[w,t] = NMF_SNC_Factors(Yvec,nFactors);
b = w(:,nFactors+1);
Ysub = Yvec - (b * t(nFactors+1,:));
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
Y = reshape(Y,d,T);
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