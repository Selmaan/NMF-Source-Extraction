%%
nA = full(sqrt(sum(A.^2)));
nr = length(nA);
A = full(A/spdiags(nA(:),0,nr,nr));
C = spdiags(nA(:),0,nr,nr)*C;
S = spdiags(nA(:),0,nr,nr)*S;

cID = find(sourceProps.idx==sourceProps.cID);
jID = find(sourceProps.idx==sourceProps.jID);
pID = find(sourceProps.idx==sourceProps.pID);
cA = A(:,cID);
jA = A(:,jID);
pA = A(:,pID);

% corJ=A(:,cID)'*A(:,jID);
% corP=A(:,cID)'*A(:,pID);
corA = A(:,cID)'*A;

%%
a=289;
b=646;
aIm = reshape(A(:,cID(a)),512,512);
bIm = reshape(A(:,b),512,512);
figure(11),imshow(cat(3,bIm/max(bIm(:)),aIm/max(aIm(:)),zeros(512)+.1)),
figure(12),imagesc(reshape(A(:,cID(a))+A(:,b),512,512))
figure(13),clf,hold on
plot(C(cID(a),:)),
plot(C(b,:)),axis tight,

corStep = 1e3;
corVals = nan(27,1);
for step=1:27
    ind = (step-1)*1e3+1:step*1e3;
    corVals(step,1) = corr(C(cID(a),ind)',C(b,ind)');
    corVals(step,2) = corr(S(cID(a),ind)',S(b,ind)');
end
%         corVals(isnan(corVals)) = 0;
fprintf('Local Corr: %0.2f \n',nanmedian(corVals(:,1))),
fprintf('Spike Corr: %0.2f \n',nanmedian(corVals(:,2))),

%%
for nROI = 1:length(cID)
% for nROI = 52:length(cID)
    nROI,
    theseOvs = find(corA(nROI,:)>0 & corA(nROI,:)<1-1e-10);
    sourceProps.cOv{nROI} = theseOvs;
    for ov = 1:length(theseOvs)
        thisDist = sqrt(sum((sourceProps.allCentroids(cID(nROI),:)-...
            sourceProps.allCentroids(theseOvs(ov),:)).^2));
        if thisDist<7
            aIm = reshape(A(:,cID(nROI)),512,512);
            bIm = reshape(A(:,theseOvs(ov)),512,512);
            figure(11),imshow(cat(3,bIm*6,aIm*6,zeros(512)+.2)),
            xlim([sourceProps.allCentroids(cID(nROI),1)-35 sourceProps.allCentroids(cID(nROI),1)]),
            ylim([sourceProps.allCentroids(cID(nROI),2)-35 sourceProps.allCentroids(cID(nROI),2)]),
            figure(12),clf,hold on
            plot(C(cID(nROI),:)),
            plot(C(theseOvs(ov),:)),
            axis tight,

            corStep = 1e3;
            corVals = nan(27,2);
            for step=1:27
                ind = (step-1)*1e3+1:step*1e3;
                corVals(step,1) = corr(C(cID(nROI),ind)',C(theseOvs(ov),ind)');
                corVals(step,2) = corr(S(cID(nROI),ind)',S(theseOvs(ov),ind)');
            end
    %         corVals(isnan(corVals)) = 0;
            fprintf('Local Corr: %0.2f \n',nanmedian(corVals(:,1))),
            fprintf('Spike Corr: %0.2f \n',nanmedian(corVals(:,2))),

            rating = input('Input 1 for merge, 0 for exclude, 2 for different cell: ');
            sourceProps.ratOv{nROI}(ov) = rating;
            sourceProps.corOv{nROI}(ov) = corA(nROI,theseOvs(ov));
            sourceProps.temOv{nROI}(ov,:) = nanmedian(corVals);
        else
            sourceProps.ratOv{nROI}(ov) = 2;
            sourceProps.corOv{nROI}(ov) = corA(nROI,theseOvs(ov));
            sourceProps.temOv{nROI}(ov,:) = nan(1,2);
        end
    end
end
        
%%

ovRank = [];
ovDist = [];
ovCor = [];
ovTem = [];
ovSpk = [];
ovIdx = [];

for nROI = 1:length(cID)
    theseOvs = sourceProps.cOv{nROI};
    ovIdx = [ovIdx, sourceProps.idx(theseOvs)'];
    ovCor = [ovCor, sourceProps.corOv{nROI}];
    for ov = 1:length(theseOvs)
        ovRank = [ovRank, sourceProps.ratOv{nROI}(ov)];
        ovDist = [ovDist, sqrt(sum((sourceProps.allCentroids(cID(nROI),:)-...
            sourceProps.allCentroids(theseOvs(ov),:)).^2))];
        
        corStep = 1e3;
        corVals = nan(27,2);
        for step=1:27
            ind = (step-1)*1e3+1:step*1e3;
            corVals(step,1) = corr(C(cID(nROI),ind)',C(theseOvs(ov),ind)');
            corVals(step,2) = corr(S(cID(nROI),ind)',S(theseOvs(ov),ind)');
        end
%         corVals(isnan(corVals)) = 0;
        
        ovTem = [ovTem,nanmedian(corVals(:,1))];
        ovSpk = [ovSpk,nanmedian(corVals(:,2))];
    end
end


%%

figure,hold on
for thisRank = 0:2
    plot(ovDist(ovRank==thisRank),ovTem(ovRank==thisRank),'.','markersize',15),
end        
        
distThresh = 5;
spkThresh = 0.2;
difOv = find(ovDist>distThresh);
exOv = find(ovDist<distThresh & ovSpk<spkThresh);
merOv = find(ovDist<distThresh & ovSpk>spkThresh);