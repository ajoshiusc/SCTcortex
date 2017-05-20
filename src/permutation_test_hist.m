function [ t2unpermuted, t2permuted_all ] = permutation_test_hist (SCT_hist_grp1, SCT_hist_grp2, numOfPermutations)

% dimensions of SCT_1_hist:
% number of histogram bins X number of subjects1 X number of mesh nodes

% dimensions of SCT_2_hist:
% number of histogram bins X number of subjects2 X number of mesh nodes

% merge SCT_1_hist and SCT_2_hist into SCT_hist
% dimensions of SCT_hist:
% number of histogram bins X (number of subjects1 + number of subjects 2) X number of mesh nodes

SCT_hist = cat (2, SCT_hist_grp1, SCT_hist_grp2);
groupLabels = [ ones(size(SCT_hist_grp1,2), 1,'uint8') ; 2*ones(size(SCT_hist_grp2,2), 1,'uint8') ];
clear SCT_1_hist SCT_2_hist;
t2permuted_all = zeros ([ size(SCT_hist,3) 1+numOfPermutations ], 'single');

parfor permNum = 1 : numOfPermutations + 1
    permNum
    %
    if permNum == 1
        % NO permutation
        groupLabelsPerm = groupLabels;
        permutation = [1:length(groupLabels)];
    else
        % permute labels
        rng (permNum);
        permutation = randperm (length (groupLabels));
        groupLabelsPerm = groupLabels (permutation);
    end
    % perform Hotelling's T2 test at each node, using permuted labels
    % t2_all = hotelling_t2_test_hist_all (SCT_hist(:,groupLabelsPerm==1,:), SCT_hist(:,groupLabelsPerm==2,:));
    
    t2_all = hotelling_t2_test_hist_kernel_all (SCT_hist(:,groupLabelsPerm==1,:), SCT_hist(:,groupLabelsPerm==2,:));
    
    % save
    t2permuted_all (:, permNum) = t2_all;
end

t2unpermuted = t2permuted_all (:,1);

return