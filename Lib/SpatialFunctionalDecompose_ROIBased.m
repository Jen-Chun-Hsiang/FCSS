function [ROI_Cluster, Converge] = SpatialFunctionalDecompose_ROIBased(IN_Score, ROI_Cluster, NumConnected)
    if nargin < 5
        NumConnected = 8;
    end
    MaxIter = 30;
    DisplayFig = 3;
    Iter = 1;
    Converge = nan(MaxIter, 1);
    ReachConvergent = 1;
    while Iter <= MaxIter && ReachConvergent <= 3
        % 1). Generate the ROI nearby matrix (all the ROIs and wherther they are
        % next to each other
        NumCluster = length(ROI_Cluster);
        fprintf('Original Number of ROIs: %d \n', NumCluster);
        ROINearbyMat = nan(NumCluster, NumCluster);
        for i = 1:NumCluster
            for j = 1:NumCluster
                if j >= i
                    ROINearbyMat(i, j) = IsROINextToEach(ROI_Cluster{i}, ROI_Cluster{j});
                else
                    ROINearbyMat(i, j) = ROINearbyMat(j, i);
                end
            end
        end
        
        ROI_Score = nan(NumCluster, size(IN_Score, 2));
        for i = 1:NumCluster
            ROI_Score(i, :) = median(IN_Score(ROI_Cluster{i}(:), :), 1);
        end
        
        Idx = LouvainJaccard(ROI_Score, NumConnected);
        
        NewCluster = unique(Idx);
        NumNewCluster = length(NewCluster);
        RmvClusterIdx = [];
        AddROIIter = 1;
        for i = 1:NumNewCluster
            CurROIIdx = find(Idx == NewCluster(i));
            NumROIinCluster = numel(CurROIIdx);
            
%             % Visual inspection for the nearby ROIs - before
%             SpatialIdx = zeros(size(ROI_Cluster{1}));
%             for c = 1:NumROIinCluster
%                 SpatialIdx(ROI_Cluster{CurROIIdx(c)}(:)) = c;
%             end
%             figure; imagesc(SpatialIdx);colorbar;
%             title('Before');
            %
            if NumROIinCluster > 1
                ROICombine = nan(NumROIinCluster, 1);
                AssignedGroup = 1;
                for j = 1:NumROIinCluster
                    for k = 1:NumROIinCluster
                        if k > j && isnan(ROICombine(j))
                            if ROINearbyMat(CurROIIdx(j), CurROIIdx(k))
                                ROICombine(j) = AssignedGroup;
                                ROICombine(k) = AssignedGroup;
                                AssignedGroup = AssignedGroup + 1;
                            end
                        end
                        if k == NumROIinCluster && isnan(ROICombine(j))
                             ROICombine(j) = AssignedGroup;
                             AssignedGroup = AssignedGroup + 1;
                        end
                    end
                end
                ROICombineNum = unique(ROICombine);
                for j = 1:length(ROICombineNum)
                    CombineIdx = find(ROICombine == ROICombineNum(j));
                    if length(CombineIdx) > 1
                        NewROICluster = zeros(size(ROI_Cluster{1}));
                        for k = 1:length(CombineIdx)
                            NewROICluster = NewROICluster | ROI_Cluster{CurROIIdx(CombineIdx(k))};
                        end
                        RmvClusterIdx = [RmvClusterIdx; CurROIIdx(CombineIdx)];
                        AddROICluster{AddROIIter} = NewROICluster;
                        AddROIIter = AddROIIter + 1;
                        fprintf('Merging... \n');                        
                    end
                end
            end            
        end
        if ~ isempty(RmvClusterIdx)
            ROI_Cluster(RmvClusterIdx) = [];
            ROI_Cluster = [ROI_Cluster AddROICluster];
            clear AddROICluster
            Converge(Iter) = 0;
        else
            Converge(Iter) = 1;
            disp('Segementation is converged. Stop the iteration.');             
            ReachConvergent = ReachConvergent + 1;
        end
        
        
%         % Visual inspection for the nearby ROIs - after
%         SpatialIdx = zeros(size(ROI_Cluster{1}));
%         for i = 1:length(ROI_Cluster)
%             SpatialIdx(ROI_Cluster{i}(:)) = i;
%         end
%         figure; imagesc(SpatialIdx);colorbar;
%         title('After');
%         %
        Iter = Iter + 1;
        fprintf('Iter: %d \n', Iter');        
    end
 