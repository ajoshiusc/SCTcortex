function [ testStats fit ] = hypothesis_test_regression_all (sctGram_all, sct_all, cogScore)

% (1)
% at each voxel, perform regression (Nadaraya-Watson)
% optimize bandwidth parameter (LOO cross validation) to optimize expected prediction error
% compute prediction error, for each subject

numOfSubjects = max (size (sct_all, 2), size (sctGram_all, 2));
numOfVertices = max (size (sct_all, 3), size (sctGram_all, 3));
subsetSize = 1024;
numOfSubsets = ceil (numOfVertices / subsetSize);

fit.valuesPredicted_all  = nan (numOfSubjects, numOfVertices, 'single');
fit.predictionErrors_all = nan (numOfSubjects, numOfVertices, 'single');
fit.predictionRMSE_all   = nan (numOfVertices, 1, 'single');
fit.conc_all             = nan (numOfVertices, 1, 'single');

for subsetNum = 1 : numOfSubsets
    % subsetNum
    % select a subset of the vertices
    indicesVertex = 1+(subsetNum-1)*subsetSize : min(numOfVertices,subsetNum*subsetSize);
    % regression
    
    %[ valuesPredicted_part, predictionErrors_part, predictionRMSE_part, conc_part ] ...
    %    = kernel_regression_hyperSphere_all (sctGram_all, indicesVertex, cogScore);
    
    [ valuesPredicted_part, predictionErrors_part, predictionRMSE_part ] ...
        = linear_regression_all (sct_all, indicesVertex, cogScore);
    
    % save
    fit.valuesPredicted_all  (:,indicesVertex) = valuesPredicted_part;
    fit.predictionErrors_all (:,indicesVertex) = predictionErrors_part;
    fit.predictionRMSE_all     (indicesVertex) = predictionRMSE_part;
    % fit.conc_all               (indicesVertex) = conc_part;
end

% (2)
% at each voxel, perform regression (constant function)
% compute prediction error, for each subject
valuesPredictedConst = mean (cogScore);
predictionErrors_const = cogScore - valuesPredictedConst;
predictionRMSE_const = sqrt (mean (predictionErrors_const.^2));

% (3)
% at each voxel, perform 2-sample F test between the previously computed 2 groups of perdiction error values
sample1 = predictionErrors_const;
for jj = 1 : numOfVertices
    if mod(jj,1e3) == 0
        % jj
    end
    %
    sample2 = fit.predictionErrors_all (:,jj);
    % F test for equality of variances
    [h,p,ci,stats] = vartest2 (sample1, sample2);
    testStats.Ftest_all   (jj) = stats.fstat;
    testStats.p_Ftest_all (jj) = p;
    
%     % Bartlett's test for equality of variances
%     [p,stats] = vartestn ([sample1', sample2], 'Display','off');
%     testStats.Bartlett_all   (jj) = stats.chisqstat;
%     testStats.p_Bartlett_all (jj) = p;
%     % Brown-Forsythe's test for equality of variances
%     [p,stats] = vartestn ([sample1', sample2], 'TestType', 'BrownForsythe', 'Display','off');
%     testStats.BrownForsythe_all   (jj) = stats.fstat;
%     testStats.p_BrownForsythe_all (jj) = p;
%     % Levene's test (absolute) for equality of variances
%     [p,stats] = vartestn ([sample1', sample2], 'TestType', 'LeveneAbsolute', 'Display','off');
%     testStats.LeveneAbsolute_all   (jj) = stats.fstat;
%     testStats.p_LeveneAbsolute_all (jj) = p;
%     % Ansari-Bradley test for equality of variances
%     [h,p,stats] = ansaribradley (sample1, sample2);
%     testStats.AnsariBradley_all   (jj) = stats.W;
%     testStats.p_AnsariBradley_all (jj) = p;
end

return