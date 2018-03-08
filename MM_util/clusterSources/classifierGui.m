classdef classifierGui < handle
    
    properties
        sources
        labels
        classNames
        centroids
        currentClass
        hFig
        hImg
        hAx
        overlayImg % An optional image to overlay the sources, e.g. for red label in inhibitory neurons.
    end
    
    methods
        function o = classifierGui(sources, labels, centroids)
            o.sources = sources;
            o.labels = labels;
            o.centroids = centroids-12.5; % Centroids in selmaans cluster function are shifted.
            o.currentClass = 1;
            o.drawFig;
            
            % These names are the classification scheme that I (Matthias)
            % am using. To add more classes, append them at the end.
            o.classNames = {'doughnutSoma', 'outOfPlaneSoma', ...
                'smallRoundProcess', 'complexProcess', ...
                'messyAndMultiples', 'cutEdge'};
            
            if any(isnan(o.labels))
                nan2newLabel(o, 'originalNans')
            end
        end
        
        function drawFig(o)
            o.hFig = figure(234234);
            clf
            o.hAx = gca;

            if isempty(o.overlayImg)
                o.hImg = imshow(mat2gray(full(...
                    reshape(sum(o.sources(:,o.labels==o.currentClass), 2), 512, 512))), ...
                    'InitialMagnification', 'fit');
            else
                cIm = zeros(512,512,3);
                cIm(:,:,1) = o.overlayImg;
                cIm(:,:,2) = mat2gray(full(...
                    reshape(sum(o.sources(:,o.labels==o.currentClass), 2), 512, 512)));
                o.hImg = imshow(cIm, 'InitialMagnification', 'fit');
            end
            
%             cIm = zeros(512,512,3);
%             nc = 100;
%             clrs = hsv(nc);
%             sel = o.labels==o.currentClass;
% 
%             % Randomize reproducibly:
%             ord = o.centroids(:, 1)-floor(o.centroids(:, 1));
%             [~, ord] = sort(ord);
% 
%             for i = find(sel(:))'
%     
%                 ic = mod(ord(i), nc);
%     
%                 if ic==0
%                     ic = nc;
%                 end
%                 aHere = mat2gray(full(reshape(o.sources(:,i), 512, 512)));
%                 cIm = cIm + bsxfun(@times, aHere, permute(clrs(ic, :), [1, 3, 2]));
%             end
%             o.hImg = imshow(cIm, 'initialMagnification', 'fit');
            
            o.hImg.ButtonDownFcn = @o.cbRemoveSource;
%             hold on
%             isDisplayed = o.labels == o.currentClass;
%             plot(o.centroids(isDisplayed, 1)-12.5, o.centroids(isDisplayed, 2)-12.5, '.')
            if numel(o.classNames) >= o.currentClass
                title(o.classNames{o.currentClass});
            end
        end
        
        function setCurrentClass(o, c)
            o.currentClass = c;
            o.drawFig;
        end
        
        function cbRemoveSource(o, ~, ~)
            isDisplayed = o.labels == o.currentClass;
            coords = abs(bsxfun(@minus, o.centroids, o.hAx.CurrentPoint(1, [2, 1])));
            dist = hypot(coords(:, 1), coords(:, 2));
            dist(~isDisplayed) = inf;
            [~, i] = min(dist);
            o.labels(i) = nan;
            o.drawFig;
        end
        
        function nan2newLabel(o, newLabelName)
            o.labels(isnan(o.labels)) = max(o.labels)+1;
            o.classNames{end+1} = newLabelName;
        end
        
        function nan2label(o, label)
            o.labels(isnan(o.labels)) = label;
        end
        
        function saveData(o, fileName)
            sources = o.sources;
            labels = o.labels;
            classNames = o.classNames;
            save(fileName, 'sources', 'labels', 'classNames', '-v7.3')
        end
    end
end
