function [ pValUnpermuted, testStatsPermutedExtreme ] = permutation_test_p_values (testStatsUnpermuted, testStatsPermuted_all, outlierRemovalPrctile)

% testStatsUnpermuted   : numOfVertices X 1
% testStatsPermuted_all : numOfVertices X numOfPermutations
% pValUnpermuted : numOfVertices X 1
% t2PermutedExtreme : numOfPermutations X 1

numOfPermutations = size (testStatsPermuted_all, 2);
for permNum = 1 : numOfPermutations 
    % store the extreme (most positive) testStat value over all nodes, BUT exclude outliers
    testStatsPermutedExtreme (permNum) = prctile (testStatsPermuted_all (:,permNum), outlierRemovalPrctile);
end

% use histogram of extreme T2 values to compute p values at all nodes
for nodeNum = 1 : size (testStatsUnpermuted, 1)
    pValUnpermuted (nodeNum) = sum (testStatsPermutedExtreme > testStatsUnpermuted (nodeNum)) / numOfPermutations;
end

return