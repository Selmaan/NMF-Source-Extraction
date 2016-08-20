function [w,t] = NMF_SNC_Factors(Yvec,nFactorsInit)

%% get corrmat
corrMat =corrcoef(Yvec');
invC = 1-corrMat;
pilC = median(invC(~isnan(invC(:))));
corrMat = exp(-1/(1*pilC^2) * invC.^2);

%% Initialize Factors

iD = sqrt(diag(1./sum(corrMat)));
A = iD*corrMat*iD;
w = rand(size(corrMat,1),nFactorsInit+1)/10;
wSym = pinv(w)*A;
its = 0;
converged = 0;
convThresh = 1e-6;
maxIts = 2e2;
wBaseline = 25;
while ~converged && its<maxIts
    its = its+1;
    wSym0 = wSym;
    [w,wSym] = updateSimNMF(A,wSym,nFactorsInit,wBaseline);
    wSym(wSym<0) = 0;
    wSym = sqrt(w'.*wSym);
    converged = mean(abs(wSym0(:)-wSym(:)))<convThresh;
end

%% Refine Factors
t = pinv(w)*Yvec;
its = 0;
converged = 0;
convThresh = 1e-4;
maxIts = 1e2;
wBaseline = 10;
while ~converged && its<maxIts
    its = its+1;
    w0 = w;
    [w,t] = updateSimNMF(Yvec,t,nFactorsInit,wBaseline);
    converged = mean(abs(w0(:)-w(:)))<convThresh;
end

%% Sparsify Factors
its = 0;
converged = 0;
convThresh = 1e-5;
maxIts = 1e3;
wBaseline = 80;
while ~converged && its<maxIts
    its = its+1;
    w0 = w;
    [w,t] = updateSimNMF(Yvec,t,nFactorsInit,wBaseline);
    converged = mean(abs(w0(:)-w(:)))<convThresh;
end