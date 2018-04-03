function [w,t,nFactorsInit] = NMF_SNC_Factors(Yvec,nFactorsInit,initImages)

epsilon = 1e-10; % To avoid division by zero.

%% get corrmat
corrMat =corrcoef(Yvec');
corrMat(isnan(corrMat)) = 0; % If two traces have exactly the same values, their corrcoef is nan.
invC = 1-corrMat;
pilC = median(invC(~isnan(invC(:))));
corrMat = exp(-1/(1*pilC^2 + epsilon) * invC.^2);

%% Initialize Factors

iD = sqrt(diag(1./(sum(corrMat) + epsilon)));
A = iD*corrMat*iD;
w = rand(size(corrMat,1),nFactorsInit+1)/10;
wSym = pinv(w)*A;
its = 0;
converged = 0;
convThresh = 1e-6;
maxIts = 2e2;
wBaseline = 10;
while ~converged && its<maxIts
    its = its+1;
    wSym0 = wSym;
    [w,wSym] = updateSimNMF(A,wSym,nFactorsInit,wBaseline);
    wSym(wSym<0) = 0;
    wSym = sqrt(w'.*wSym);
    converged = mean(abs(wSym0(:)-wSym(:)))<convThresh;
end

%% Add factors from images

if ~isempty(initImages)
    useIm = zeros(size(initImages,3),1);
    for nImg = 1:size(initImages,3)
        if bwarea(bwareafilt(medfilt2(initImages(:,:,nImg))>0,1)) > 50
            useIm(nImg) = 1;
        else
            useIm(nImg) = 0;
        end
    end
    
    if sum(useIm)>0
        imgVec = reshape(initImages,size(w,1),size(initImages,3));
        w = [w(:,1:nFactorsInit),imgVec(:,find(useIm)),w(:,nFactorsInit+1)];
        nFactorsInit = nFactorsInit + sum(useIm);
    end
end

% Normalize w
w = bsxfun(@rdivide,w,sqrt(sum(w.^2)));
w(~isfinite(w)) = 0;
%% Refine Factors
t = pinv(w)*Yvec;
its = 0;
converged = 0;
convThresh = 2e-5;
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
convThresh = 2e-6;
maxIts = 250;
wBaseline = 80;
while ~converged && its<maxIts
    its = its+1;
    w0 = w;
    [w,t] = updateSimNMF(Yvec,t,nFactorsInit,wBaseline);
    converged = mean(abs(w0(:)-w(:)))<convThresh;
end