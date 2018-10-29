function ReducedData = TemporalDimensionReduction(Data, Explanation)
    % Data: [Pixels, Timepoints]
    % Explanation: The percentage of the variation explained by reduced
    % dimension (0 ~ 1)
    
    [Coeff, Score, Latent, ~, Explained] = pca(Data);
    CumExplainedVar = cumsum(Explained);
    DimRdc = find(CumExplainedVar>Explanation*100);
    DimRdc = DimRdc(1);
    figure; plot(Explained, 'b'); hold on;
    plot(CumExplainedVar, 'g');hold on;
    plot(DimRdc, CumExplainedVar(DimRdc), 'rx');
    legend({'Explained Dimemsions', 'Cumulative Explanation', 'Selected Point'});
    
    ReducedData = Score(:, 1:DimRdc);
end