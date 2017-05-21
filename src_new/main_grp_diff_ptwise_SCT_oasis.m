clc;clear all;close all;

% This script performs pointwise group differences in cortical thicknesses
% for the left hemisphere
% modify as appropriate for left hemisphere

%group1 and group2 directory names
addpath(genpath('C:\Users\ajoshi\Documents\coding_ground\svreg-matlab\dev'));
addpath(genpath('C:\Users\ajoshi\Documents\coding_ground\svreg-matlab\src'));
% 
% dirname1='C:\Users\ajoshi\Downloads\oasis_data';
% dirname2='C:\Users\ajoshi\Downloads\Beijing';
% NLevels=1;
% aa1=dir(dirname1);
% aa1=aa1(3:end);
% aa2=dir(dirname2);
% aa2=aa2(3:end);

jjj=1;
npts=1000;
tar_orig=readdfs('C:\Users\ajoshi\Documents\coding_ground\svreg-matlab/BrainSuiteAtlas1\mri.left.mid.cortex.dfs');
tar=reducepatch(tar_orig,10000);
[~,ia,indx_reduced]=intersect(tar.vertices,tar_orig.vertices,'rows','stable')


[N,T,R]=xlsread('E:\sipi_data\oasis_full_data\data\oasis_cross-sectional.csv');
fnames= R(2:end,1);
temp  = R(2:end,8);
for jj=1:length(fnames)
    CDRs(jj) = temp{jj};
end

% get size of each histogram to preallocate
mat_fname=['E:\sipi_data\oasis_full_data\data\',fnames{1},'\PROCESSED\MPRAGE\T88_111\',fnames{1},'_aaj.SCT.mat'];
tmp = load (mat_fname);
numOfScales = size (tmp.SCTatlasleft, 2)
numOfMeshNodes = length(tar.vertices);%size (tmp.SCTatlas_left, 2)

% load SCT histograms for group 1 (controls, label 0 in OASIS)
numOfSubjectsGrp1 = sum (CDRs == 0);
fnamesGrp1 = fnames (CDRs == 0);
SCTgrp1 = zeros (numOfScales, numOfMeshNodes, numOfSubjectsGrp1, 'single');
count = 0;
% generate SCT histograms for group 1, save them
% parfor jj=1:numOfSubjectsGrp1
%      subbasename=['E:\sipi_data\oasis_full_data\data\',fnamesGrp1{jj},'\PROCESSED\MPRAGE\T88_111\',fnamesGrp1{jj},'_aaj'];
%      sub_SCT_hist(subbasename);
%      fprintf('%d/%d SCT done\n',jj,numOfSubjectsGrp1);
% end
 % generate SCT histograms for group 2, save them
numOfSubjectsGrp2 = sum (CDRs == 1);
fnamesGrp2 = fnames (CDRs == 1);
count = 0;
SCTgrp2 = zeros (numOfScales, numOfMeshNodes, numOfSubjectsGrp2, 'single');

%  parfor jj=1:numOfSubjectsGrp2
%      subbasename=['E:\sipi_data\oasis_full_data\data\',fnamesGrp2{jj},'\PROCESSED\MPRAGE\T88_111\',fnamesGrp2{jj},'_aaj'];
%      sub_SCT_hist(subbasename);
%      fprintf('%d/%d SCT done\n',jj,numOfSubjectsGrp2);
% end

for jj = 1 : numOfSubjectsGrp1
    mat_fname=['E:\sipi_data\oasis_full_data\data\',fnamesGrp1{jj},'\PROCESSED\MPRAGE\T88_111\',fnamesGrp1{jj},'_aaj.SCT.mat'];
    if ~exist(mat_fname,'file')
        continue;
    end
    mat_fname
    count = count + 1
    tmp = load (mat_fname, 'SCTatlasleft');
    SCTgrp1 (:,:,count) = single (tmp.SCTatlasleft (indx_reduced,:))';
end
SCTgrp1 = SCTgrp1 (:,:,1:count);
SCTgrp1 = permute (SCTgrp1, [1 3 2]);

% load SCT histograms for group 2 (moderate AD, label 1 in OASIS)
numOfSubjectsGrp2 = sum (CDRs == 1);
fnamesGrp2 = fnames (CDRs == 1);
SCTgrp2 = zeros (numOfScales, numOfMeshNodes, numOfSubjectsGrp2, 'single');
count = 0;
for jj = 1 : numOfSubjectsGrp2
    mat_fname=['E:\sipi_data\oasis_full_data\data\',fnamesGrp2{jj},'\PROCESSED\MPRAGE\T88_111\',fnamesGrp2{jj},'_aaj.SCT.mat'];
    if ~exist(mat_fname,'file')
        continue;
    end
    mat_fname
    count = count + 1
    tmp = load (mat_fname, 'SCTatlasleft');
    SCTgrp2 (:,:,count) = single (tmp.SCTatlasleft (indx_reduced,:))';
end
SCTgrp2 = SCTgrp2 (:,:,1:count);
SCTgrp2 = permute (SCTgrp2, [1 3 2]);

save ('data_oasis_sct_left.mat', 'SCTgrp1', 'SCTgrp2', 'indx_reduced','tar');

return

load ('data_oasis_sct_right.mat');
%SCTgrp2=SCTgrp2(:,1:18,:);
tar_orig=readdfs('C:\Users\ajoshi\Documents\coding_ground\svreg-matlab/BrainSuiteAtlas1\mri.right.mid.cortex.dfs');
%tar=reducepatch(tar_orig,10000);
tarsm=smooth_cortex_fast(tar,.1,200);
cmap=cbrewer('seq','Oranges',100);

% SCT: call permutation testing for histograms
numOfPermutations = 10
[ t2Unpermuted t2permuted_all ] = permutation_test(SCTgrp1, SCTgrp2(:,1:10,:), numOfPermutations);
patch('vertices',tarsm.vertices,'faces',tar.faces,'facevertexcdata',t2Unpermuted,'edgecolor','none','facecolor','interp');axis equal;axis off;camlight;material dull;lighting phong; colormap jet; colorbar;
save ([ 'results_oasisright_10.mat']);
outlierRemovalPrctile = 97
[ pValUnpermuted t2PermutedExtreme ] = permutation_test_p_values (t2Unpermuted, t2permuted_all, outlierRemovalPrctile);
%patch('vertices',tarsm.vertices,'faces',tar.faces,'facevertexcdata',pValUnpermuted','edgecolor','none','facecolor','interp');axis equal off tight;material dull;lighting phong; colormap jet;view(90,0);camlight;
show_surf_val(tarsm,1-pValUnpermuted,[0,1],cmap)
export_fig -a1 -m2 -transparent Pics/results_oasis_SCT_right_hemi_10.png, close

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% S histogram only
Sgrp1 = SCTgrp1(1:5,:,:);%squeeze (sum (sum (reshape (SCTgrp1, 5, 5, 5, size(SCTgrp1,2), size(SCTgrp1,3)), 2), 3));
Sgrp2 = SCTgrp2(1:5,:,:);%squeeze (sum (sum (reshape (SCTgrp2, 5, 5, 5, size(SCTgrp2,2), size(SCTgrp2,3)), 2), 3));
[ S_t2Unpermuted S_t2permuted_all ] = permutation_test(Sgrp1, Sgrp2, numOfPermutations);
outlierRemovalPrctile = 97
[ S_pValUnpermuted S_t2PermutedExtreme ] = permutation_test_p_values (S_t2Unpermuted, S_t2permuted_all, outlierRemovalPrctile);
show_surf_val(tarsm,1-S_pValUnpermuted,[0,1],cmap)
export_fig -a1 -m2 -transparent Pics/results_oasis_S_left_hemi.png, close

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% C histogram only
Cgrp1 = SCTgrp1(6:10,:,:);%squeeze (sum (sum (reshape (SCTgrp1, 5, 5, 5, size(SCTgrp1,2), size(SCTgrp1,3)), 1), 3));
Cgrp2 = SCTgrp2(6:10,:,:);%squeeze (sum (sum (reshape (SCTgrp2, 5, 5, 5, size(SCTgrp2,2), size(SCTgrp2,3)), 1), 3));
[ C_t2Unpermuted C_t2permuted_all ] = permutation_test(Cgrp1, Cgrp2, numOfPermutations);
outlierRemovalPrctile = 97
[ C_pValUnpermuted C_t2PermutedExtreme ] = permutation_test_p_values (C_t2Unpermuted, C_t2permuted_all, outlierRemovalPrctile);
show_surf_val(tarsm,1-C_pValUnpermuted,[0,1],cmap)
export_fig -a1 -m2 -transparent Pics/results_oasis_C_left_hemi.png, close

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% T histogram only
Tgrp1 = SCTgrp1(11:15,:,:);%squeeze (sum (sum (reshape (SCTgrp1, 5, 5, 5, size(SCTgrp1,2), size(SCTgrp1,3)), 1), 2));
Tgrp2 = SCTgrp2(11:15,:,:);%;%squeeze (sum (sum (reshape (SCTgrp2, 5, 5, 5, size(SCTgrp2,2), size(SCTgrp2,3)), 1), 2));
[ T_t2Unpermuted T_t2permuted_all ] = permutation_test(Tgrp1, Tgrp2, numOfPermutations);
outlierRemovalPrctile = 97
[ T_pValUnpermuted T_t2PermutedExtreme ] = permutation_test_p_values (T_t2Unpermuted, T_t2permuted_all, outlierRemovalPrctile);
show_surf_val(tarsm,1-T_pValUnpermuted,[0,1],cmap)
export_fig -a1 -m2 -transparent Pics/results_oasis_T_left_hemi.png, close

return
