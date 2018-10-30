function Community = LouvainJaccard(Dataset, k)
    
    [NumNode, NumDim] = size(Dataset);
    if nargin == 1
        k = floor(log(NumNode)*3);
    end
    
    
    kNNMat = knnsearch(Dataset, Dataset, 'k', k);
    
    AdjMat = nan(NumNode, NumNode);
    JaccardInd = @(x, y) numel(intersect(x, y))/numel(union(x, y));
    
    for i = 1:NumNode
        for j = 1:NumNode
            if isnan(AdjMat(j, i))
                Neibor_i = kNNMat(i, :);
                Neibor_j = kNNMat(j, :);
                AdjMat(i, j) = JaccardInd(Neibor_i, Neibor_j);
            else
                AdjMat(i, j) = AdjMat(j, i);
            end
        end
    end
    
    [Community, ~] = fast_mo(AdjMat);
%     Coms = unique(Community);
%     NumCom = length(Community);

end