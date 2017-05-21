function [ valuesPredicted_all, predictionErrors_all, predictionRMSE_all ] ...
    = linear_regression_all (features, indicesVertex, values)

% independent variable = features (feature-dimension X number of subjects X number of vertices)
% dependent variable = values (1 X number of subjects)

features = features (:,:,indicesVertex);

numOfHistBins = size (features, 1);
numOfSubjects = size (features, 2);
numOfVertices = size (features, 3);

% dependent variable
Y = values';
Y_all = zeros (numOfSubjects * numOfVertices, 1, 'double');
for ii = 1 : numOfVertices
   Y_all ((ii-1) * numOfSubjects + 1 : ii * numOfSubjects) = double (Y);
end
% independent variables (add a constant vector to allow a shift)
allOnes = ones (numOfSubjects, 1, 'double');
features_cell = cell (numOfVertices,1);
for ii = 1 : numOfVertices
  features_cell{ii} = sparse ([ double(features(:,:,ii)') allOnes ]);
end
X_all = blkdiag (features_cell{:});
% prevent too many warnings
%warning off stats:regress:RankDefDesignMat
%[B_all,BINT_all,R_all,RINT_all,STATS_all] = regress (Y_all, X_all);
warning off MATLAB:rankDeficientMatrix
B_all = X_all \ Y_all;
% B = regress (Y, X) gives coefficients B such that Y = B * X
%
valuesPredicted_all = reshape (full (single (X_all * B_all)), numOfSubjects, numOfVertices);
clear X_all B_all Y_all

% for jj = indicesVertex
%     idx = jj - indicesVertex (1) + 1;
%     % independent variables (add a constant vector to allow a shift)
%     X = [ features(:,:,idx)' allOnes ];
%     % prevent too many warnings
%     warning off stats:regress:RankDefDesignMat
%     [B,BINT,R,RINT,STATS] = regress (Y, X);
%     % B = regress (Y, X) gives coefficients B such that Y = B * X
%     %
%     valuesPredicted_all (:,idx) = X * B;
% end

% prediction errors for the optimal coefficients B
predictionErrors_all = bsxfun (@minus, valuesPredicted_all, values');
predictionRMSE_all = sqrt (sum (predictionErrors_all.^2) / numOfSubjects);

return