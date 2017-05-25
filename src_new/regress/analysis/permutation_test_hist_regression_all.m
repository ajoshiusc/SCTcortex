function [ testStatsUnpermuted_all, ...
    testStatsPermuted_Ftest_all, ...
    testStatsPermuted_Bartlett_all, ...
    testStatsPermuted_BrownForsythe_all, ...
    testStatsPermuted_LeveneAbsolute_all, ...
    testStatsPermuted_AnsariBradley_all ] ...
    = permutation_test_hist_regression_all (histsSqrtCrossExp_all, hists, cogScore, numOfPermutations)

% dimensions of histsCross_all
% number of subjects X number of subjects X number of mesh vertices

% dimensions of cogScore
% 1 X number of subjects

numOfVertices = max ( size (histsSqrtCrossExp_all,3), size (hists,3) );

testStatsPermuted_Ftest_all          = nan ([ numOfVertices 1+numOfPermutations ], 'single');
testStatsPermuted_Bartlett_all       = nan ([ numOfVertices 1+numOfPermutations ], 'single');
testStatsPermuted_BrownForsythe_all  = nan ([ numOfVertices 1+numOfPermutations ], 'single');
testStatsPermuted_LeveneAbsolute_all = nan ([ numOfVertices 1+numOfPermutations ], 'single');
testStatsPermuted_AnsariBradley_all  = nan ([ numOfVertices 1+numOfPermutations ], 'single');

parfor permNum = 1 : numOfPermutations + 1
    permNum
    if permNum == 1
        % NO permutation
        cogScorePerm = cogScore;
    else
        % permute cogScore
        rng (permNum);
        cogScorePerm = cogScore (randperm (length (cogScore)));
    end
    % perform hypothesis test at each node, using permuted labels
    [ testStats_all ] = hypothesis_test_hist_regression_all (histsSqrtCrossExp_all, hists, cogScorePerm);
    % save
    testStatsPermuted_Ftest_all          (:,permNum) = testStats_all.Ftest_all';
%     testStatsPermuted_Bartlett_all       (:,permNum) = testStats_all.Bartlett_all';
%     testStatsPermuted_BrownForsythe_all  (:,permNum) = testStats_all.BrownForsythe_all';
%     testStatsPermuted_LeveneAbsolute_all (:,permNum) = testStats_all.LeveneAbsolute_all';
%     testStatsPermuted_AnsariBradley_all  (:,permNum) = testStats_all.AnsariBradley_all';
end
%
%save ('sct_hist_simulation_testStats_right_all.mat', ...
%    'testStatsPermuted_Ftest_all', ...
%    'testStatsPermuted_Bartlett_all', ...
%    'testStatsPermuted_BrownForsythe_all', ...
%    'testStatsPermuted_LeveneAbsolute_all', ...
%    'testStatsPermuted_AnsariBradley_all' );%

testStatsUnpermuted_all.Ftest_all          = testStatsPermuted_Ftest_all          (:,1);
testStatsUnpermuted_all.Bartlett_all       = testStatsPermuted_Bartlett_all       (:,1);
testStatsUnpermuted_all.BrownForsythe_all  = testStatsPermuted_BrownForsythe_all  (:,1);
testStatsUnpermuted_all.LeveneAbsolute_all = testStatsPermuted_LeveneAbsolute_all (:,1);
testStatsUnpermuted_all.AnsariBradley_all  = testStatsPermuted_AnsariBradley_all  (:,1);

return