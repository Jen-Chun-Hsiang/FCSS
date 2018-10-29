function OUT = IsROINextToEach(ROIa, ROIb)
    Cluster_a = bwconncomp(ROIa);
    Cluster_b = bwconncomp(ROIb);
    if Cluster_a.NumObjects ~= 1 || Cluster_b.NumObjects ~= 1
        error('The number of input ROI should only have one pair at a time!');
    end
    ROI_merge = ROIa | ROIb;
    Cluster_merge = bwconncomp(ROI_merge);
    if Cluster_merge.NumObjects == 1
        OUT = 1;
    else
        OUT = 0;
    end
end