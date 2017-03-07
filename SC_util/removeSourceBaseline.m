function [fSub, f_] = removeSourceBaseline(f,binSize,basePrct)

nF = numel(f);
fOrig = f;
f(nF+1-rem(nF,binSize):end)=[];
f = reshape(f, binSize, []);
baseEstimate = prctile(f, basePrct, 1);

x = 1:nF;
xBin = (binSize:binSize:length(baseEstimate)*binSize) - binSize/2;
warnState = warning('off', 'stats:statrobustfit:IterationLimit');
try
    b = robustfit(xBin',baseEstimate', 'bisquare', 2);
catch
    keyboard,
end
warning(warnState);

f = fOrig;
f_ = [ones(nF, 1), x'] * b;
f_ = f_';
fSub = f - f_;
end