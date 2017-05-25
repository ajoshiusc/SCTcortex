function [ t2unpermuted t2permuted_all ] = permutation_test (SCT1, SCT2, numOfPermutations)

% dimensions of SCT1:
% number of features X number of subjects1 X number of mesh nodes

% dimensions of SCT2:
% number of features X number of subjects2 X number of mesh nodes

% merge SCT1 and SCT 2 into SCT
% dimensions of SCT:
% number of features X (number of subjects1 + number of subjects 2) X number of mesh nodes

SCT = cat (2, SCT1, SCT2);
groupLabels = [ ones(size(SCT1,2), 1,'uint8') ; 2*ones(size(SCT2,2), 1,'uint8') ];
clear SCT1 SCT2;
t2permuted_all = zeros ([ size(SCT,3) 1+numOfPermutations ], 'single');

% save the t2 values at all nodes for the unpermuted data
permNum = 1
[~,t2unpermuted] = hotelling_t2_test_all (SCT(:,groupLabels==1,:), SCT(:,groupLabels==2,:));
t2permuted_all (:, 1) = t2unpermuted;

parfor permNum = 2 : numOfPermutations + 1
    permNum
    % permute labels
    groupLabelsPerm = groupLabels (randperm (length (groupLabels)));
    % perform Hotelling's T2 test at each node, using permuted labels
    [~,t2_all] = hotelling_t2_test_all (SCT(:,groupLabelsPerm==1,:), SCT(:,groupLabelsPerm==2,:));
    t2permuted_all (:, permNum) = t2_all;
end

return