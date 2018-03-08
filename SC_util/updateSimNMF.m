function [w,t] = updateSimNMF(oMov,t,nCells,wBaseline)

t(t<0) = 0;
try
    pT = pinv(t);
catch err
    fprintf('Number of nans in t: %d\n', sum(isnan(t(:))));
    fprintf('Number of infs in t: %d\n', sum(isinf(t(:))));
    disp(t)
    rethrow(err)
end
w = oMov*pT;
for nCell = 1:nCells
    minW = max(0,prctile(w(:,nCell),wBaseline));
    w(:,nCell) = w(:,nCell)-minW;
end
w(w<0) = 0;
% Enforce L2 Norm
w = bsxfun(@rdivide,w,sqrt(sum(w.^2)));
w(~isfinite(w)) = 0;
pW = pinv(w);
t = pW*oMov;