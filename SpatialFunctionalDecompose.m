function [OUT, h] = SpatialFunctionalDecompose(IN_Score, IN_Cluster, IN_X, IN_Y, NumConnected)
    if nargin < 5
        NumConnected = 8;
    end
    MaxIter = 15;
    DisplayFig = 3;
    OrgNumCluster = length(IN_Cluster);
    for c = 1:OrgNumCluster
        Iter = 1;
        while Iter <= MaxIter
            if Iter == 1
                Idx = LouvainJaccard(IN_Score(IN_Cluster{c}(:), :), NumConnected);                
            else
                CurScore = IN_Score;
                CurCluster = unique(Idx);
                CurNumCluster = length(CurCluster);
                for i = 1:CurNumCluster
                    CurCluObj = bwconncomp(SpatialIdx == CurCluster(i));
                    for j = 1:length(CurCluObj.PixelIdxList)
                        CurScore(CurCluObj.PixelIdxList{j}, :) = repmat(mean(CurScore(CurCluObj.PixelIdxList{j}, :), 1),...
                            numel(CurCluObj.PixelIdxList{j}), 1);
                    end
                end
                Idx = LouvainJaccard(CurScore(IN_Cluster{c}(:), :), NumConnected);
            end
            SpatialIdx = zeros(size(IN_Cluster{c}(:)));
            SpatialIdx(IN_Cluster{c}(:)) = Idx;
            SpatialIdx = reshape(SpatialIdx, IN_X, IN_Y);
            if Iter >= MaxIter-DisplayFig
                h{c, Iter-(MaxIter-DisplayFig)+1} = figure('Visible','off'); imagesc(SpatialIdx);colorbar;
            end
            Iter = Iter + 1;
            fprintf('Orginal Cluster: %d, Iter: %d \n', c, Iter');
        end
        for i = 1:CurNumCluster
            if i == 1 & c == 1
                OUT{1} = SpatialIdx == i;
            else
                OUT{end+1} = SpatialIdx == i;
            end
        end        
    end
end