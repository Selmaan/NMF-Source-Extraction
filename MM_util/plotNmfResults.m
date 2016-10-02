%%
cIm = zeros(512,512,3);
cIm(:,:,1) = reshape(sum(A(:,sourceProps.idx==sourceProps.sigRank(1)),2),512,512);
cIm(:,:,2) = reshape(sum(A(:,sourceProps.idx==sourceProps.sigRank(2)),2),512,512);
cIm(:,:,3) = reshape(sum(A(:,sourceProps.idx==sourceProps.sigRank(3)),2),512,512);
cIm(:,:,3) = cIm(:,:,3) + reshape(sum(A(:,sourceProps.idx==sourceProps.sigRank(4)),2),512,512);
cIm(:,:,2) = cIm(:,:,2) + reshape(sum(A(:,sourceProps.idx==sourceProps.sigRank(4)),2),512,512);

figure(2)
imshow(cIm*4, 'parent', gca)

%% Color randomly:

cIm = zeros(512,512,3);
bounds = [];
nc = 20;
clrs = hsv(nc);

% Randomize reproducibly:
ord = sourceProps.allCentroids(:, 1)-floor(sourceProps.allCentroids(:, 1));
[~, ord] = sort(ord);

for i = 1:size(A, 2)
    
    % Exclude bad sources:
    if ismember(sourceProps.idx(i), sourceProps.sigRank([3]))
        continue
    end
    
    ic = mod(ord(i), nc);
    
    if ic==0
        ic = nc;
    end
    aHere = full(reshape(A(:, i), 512, 512));
    imgHere = bsxfun(@times, aHere, ...
        permute(clrs(ic, :), [1, 3, 2]));
    cIm = cIm + imgHere;
    
    bounds = cat(1, bounds, bwboundaries(aHere>0));
end

figure(1)
clf
imshow(cIm*4, 'parent', gca)

hold on
for i = 1:numel(bounds)
    plot(bounds{i}(:, 2), bounds{i}(:, 1), 'color', [0.3 0.3 0.3], ...
        'linewidth', 1)
end

%% Plot traces
sel = sourceProps.idx==sourceProps.sigRank(2);
figure(99)
clf
hold on
offset = 0;
for i = 1:numel(sel)
    
    if ~sel(i)
        continue
    end
    
    plot(C(i, :) + offset)
    
    offset = offset + max(C(i, :));
end
    

