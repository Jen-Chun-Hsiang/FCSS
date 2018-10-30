% USAGE INSTRUCTION
% 1). Load the registered, denoised time-lapse images (in a 3D array) in a variable named 'I' 
%   - In this demo, I load it from 'Demo.mat'
% 2). After intensity function clustering (~ 25 mins), you need to pick those clusters
% that are not the background nor the noise-contaminated pixels (a dialog box will pop up)
%   - In this demo, I picked the first 3 clusters, and skip the rest.
% 3). The final result (total time ~ 45 mins) is that the pixel-wised ROIs will be stored in a variable called
% 'SlcROI', which is a cell structure containing binary spatial info for each selected ROI
%
% MATLAB REQUIREMENT
% Please make sure you have installed 'statistics and machine learning
% toolbox', 'image processing toolbox'


%% Define Parameters
Fz = 9.48; % Sampling rate of 2-P imaging (frames/second)
MinDimByExpla = 0.8;
LJClusterK = 8;
IsProcessVisualized = 1;
CluMinSize = 4;

%% Load data and functions
% load a variable 'I' as the 3D array of the registered time-lapse 2-p imaging 
% I: [Height, Width, Time]
load('Demo.mat'); 
addpath(genpath('Lib'));
addpath(genpath('3rdPartyLib'));
%% Pixel-wise Functional Clustering
ImgClm = I;
ImgClm = reshape(ImgClm(:), size(ImgClm, 1)*size(ImgClm, 2), []); % Convert 3D to 2D (pixel, timepoints)
NumTimePoints = size(ImgClm, 2);

FilNam = 'Demo_ROIs_Cluster';
if exist(['./ROIs/' FilNam '.mat'], 'file')
    disp('Loading pixel-wise functional ROI clustering dataset...');
    load(['./ROIs/' FilNam '.mat']);
else
    % Perform Dimensional Reduction - PCA    
    ReducedImgClm = TemporalDimensionReduction(ImgClm, MinDimByExpla);
    
    % Clustering
    Idx = LouvainJaccard(ReducedImgClm, LJClusterK);
    
    % Sort clusters by the average intensity
    Clusters = unique(Idx);
    NumCluster = length(Clusters);
    ClusterMean = nan(1, NumCluster);  
    ClusterMeanTraces = nan(NumCluster, NumTimePoints);  
    for i = 1:NumCluster
        CurClusterMean = mean(ImgClm(Idx == Clusters(i), :), 1);
        ClusterMeanTraces(i, :) = CurClusterMean;        
        ClusterMean(i) = mean(CurClusterMean);
    end
    [~, SortedIdx] = sort(ClusterMean, 'descend');
    
    if IsProcessVisualized
        Colors = parula(NumCluster);  
        figure;
        subplot(1, 2, 1);
        SpatialIdx = reshape(Idx, size(I, 1), size(I, 2));
        imagesc(SpatialIdx);colormap(Colors);colorbar;
        
        subplot(1, 2, 2);
        T = (0:NumTimePoints-1)/Fz;
        ClusterMeanTraces = ClusterMeanTraces./range(ClusterMeanTraces, 2);
        for i = 1:NumCluster
            plot(T', ClusterMeanTraces(i, :)', 'Color', Colors(i,:)); hold on
        end
    end
    
    save(['./ROIs/' FilNam '.mat'],...
        'Idx', 'SortedIdx', 'SpatialIdx', 'T', 'Clusters', 'NumCluster', 'Colors', 'ClusterMean', 'ClusterMeanTraces');
end

%% Manually select the responsive clusters and perform further functional segmentation
FilNam = 'Demo_ROIs_SemiCluster';

clear SelecClusterMask
% Original (without rerun)
if exist(['./ROIs/' FilNam '.mat'], 'file')
    disp('Loading visually-selected pixel-wise functional ROI clustering dataset...');
    load(['./ROIs/' FilNam '.mat']);
else
    if ~exist('SelecClusterMask', 'var')
        % Display The Functional Contour Image
        figure;
        AvgIgC = mean(I, 3);
        imagesc(AvgIgC);colorbar;title('Average raw data for each pixel');
        
        SelecCluster = nan(1, NumCluster);
        for i = 1:NumCluster
            h1 = figure;
            plot(T', ClusterMeanTraces(SortedIdx(i),:)', 'Color', Colors(SortedIdx(i),:));
            h2 = figure;
            SpatialIdx = reshape(Idx ==Clusters(SortedIdx(i)), size(I, 1), size(I, 2));
            imagesc(SpatialIdx);colorbar;
            keyboard;
            % (Instruction) for the data from Demo.mat, I take the first three, and skip
            % the rest.
            Answer = keepCluster;
            if Answer == 2
                SelecCluster(i:end) = 0;
                break;
            else
                SelecCluster(i) = Answer;
            end
            close(h1);
            close(h2);
        end
        SelecCluster = Clusters(SortedIdx(SelecCluster==1));
    end
    clear SelecClusterMask
    for i = 1:length(SelecCluster)
        SelecClusterMask{i} = reshape(Idx ==SelecCluster(i), size(I, 1), size(I, 2));
    end
    
    % Normalize intensity for each pixel
    ImgClm = (ImgClm-quantile(ImgClm, 0.25, 2))./std(ImgClm, [], 2);
        
    % Functionally decompose the selected intensity cluster
    ReducedImgClm = TemporalDimensionReduction(ImgClm, MinDimByExpla);
    [SelecClusterMask, FigureH] = SpatialFunctionalDecompose(ReducedImgClm, SelecClusterMask, size(I, 1), size(I, 2), LJClusterK);
    
    if IsProcessVisualized
        NumCluster = length(SelecClusterMask);
        Colors = parula(NumCluster);  
        figure;
        
        SpatialIdx = nan(size(SelecClusterMask{1}));
        for i = 1:NumCluster
           SpatialIdx(SelecClusterMask{i}) = i;
        end
        imagesc(SpatialIdx);colormap(Colors);colorbar;
    end
    save(['./ROIs/' FilNam '.mat'],...
        'SelecCluster', 'SelecClusterMask', 'SpatialIdx');
end

%% Spatially separate the selected clusters for further segmentation
FilNam = 'Demo_ROIs_AutoCluster';
if exist(['./ROIs/' FilNam '.mat'], 'file')
    disp('Loading ROI dataset...');
    load(['./ROIs/' FilNam '.mat']);
else
    StdPix = std(I, [], 3);
    RmvPix = zeros(size(SelecClusterMask{1}));
    for i = 1:length(SelecClusterMask)
        CurCluObj = bwconncomp(~RmvPix & SelecClusterMask{i});
        if i == 1
            CluObj = CurCluObj;
        else
            CluObj.NumObjects = CluObj.NumObjects + CurCluObj.NumObjects;
            CluObj.PixelIdxList = [CluObj.PixelIdxList, CurCluObj.PixelIdxList];
        end
    end
    CluCelInd = cellfun(@numel, CluObj.PixelIdxList);
    CluCelInd = CluCelInd >= CluMinSize;
    CluCel = CluObj.PixelIdxList(CluCelInd);
    SelClu = zeros(size(RmvPix));
    MultiSelClu = zeros(size(RmvPix));
    for i = 1:length(CluCel)
        SelClu(CluCel{i}) = 1;
        MultiSelClu(CluCel{i}) = i;
        CurSlcROI = zeros(size(RmvPix));
        CurSlcROI(CluCel{i}) = 1;
        SlcROI{i} = CurSlcROI == 1;
    end
    
    if IsProcessVisualized
        figure; 
        subplot(1, 2, 1);
        imagesc(MultiSelClu);colorbar;
        title('Segmentation of ROIs');
        subplot(1, 2, 2);
        imagesc(StdPix.*SelClu);colorbar;
        title('STD of ROIs');
    end
    
    save(['./ROIs/' FilNam '.mat'],...
        'CluCel', 'SlcROI', 'SelClu', 'RmvPix', 'CluObj', 'CluMinSize');
end

%% Merge neighboring components exhibiting similar function and iterate until unchanged
FilNam = 'Demo_ROIs_ConvergentCluster';
if exist(['./ROIs/' FilNam '.mat'], 'file')
    disp('Loading Convergent ROI dataset...');
    load(['./ROIs/' FilNam '.mat']);
else
    ProcessingdRecord.NumROI_BeforeConverge = length(SlcROI);
    [SlcROIMerge, ProcessingdRecord.IterConverge] = SpatialFunctionalDecompose_ROIBased(ReducedImgClm, SlcROI, LJClusterK);
    ProcessingdRecord.NumROI_AfterConverge = length(SlcROIMerge);    
    clear CluCel
    
    if IsProcessVisualized
        figure; 
        subplot(1, 2, 1);
        SpatialIdx = zeros(size(SlcROI{1}));
        for i = 1:length(SlcROI)
            SpatialIdx(SlcROI{i}(:)) = i;
        end
        imagesc(SpatialIdx); colorbar;title('Before convergent merging');
        subplot(1, 2, 2);
        SpatialIdx = zeros(size(SlcROIMerge{1}));
        for i = 1:length(SlcROIMerge)
            SpatialIdx(SlcROIMerge{i}(:)) = i;
            CluCel{i} = find(SlcROIMerge{i} == 1);
        end
        imagesc(SpatialIdx); colorbar;title('After convergent merging');
    end
    SlcROI = SlcROIMerge;
    save(['./ROIs/' FilNam '.mat'],...
        'CluCel', 'SlcROI', 'ProcessingdRecord');
end