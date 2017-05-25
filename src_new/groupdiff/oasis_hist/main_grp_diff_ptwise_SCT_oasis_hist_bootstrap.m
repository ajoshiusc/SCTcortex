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
[~,ia,indx_reduced]=intersect(tar.vertices,tar_orig.vertices,'rows','stable');


[N,T,R]=xlsread('E:\sipi_data\oasis_full_data\data\oasis_cross-sectional.csv');
fnames= R(2:end,1);
temp  = R(2:end,8);
for jj=1:length(fnames)
  CDRs(jj) = temp{jj};
end

% get size of each histogram to preallocate
mat_fname=['E:\sipi_data\oasis_full_data\data\',fnames{1},'\PROCESSED\MPRAGE\T88_111\',fnames{1},'_aaj.SCT_hist.mat'];
tmp = load (mat_fname);
numOfHistBins = size (tmp.SCT_hist_atlas_left, 1);
numOfMeshNodes = length(tar.vertices);%size (tmp.SCT_hist_atlas_left, 2)

% load SCT histograms for group 1 (controls, label 0 in OASIS)
numOfSubjectsGrp1 = sum (CDRs == 0);
fnamesGrp1 = fnames (CDRs == 0);
SCT_hist_grp1 = zeros (numOfHistBins, numOfMeshNodes, numOfSubjectsGrp1, 'single');
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
SCT_hist_grp2 = zeros (numOfHistBins, numOfMeshNodes, numOfSubjectsGrp2, 'single');

%  parfor jj=1:numOfSubjectsGrp2
%      subbasename=['E:\sipi_data\oasis_full_data\data\',fnamesGrp2{jj},'\PROCESSED\MPRAGE\T88_111\',fnamesGrp2{jj},'_aaj'];
%      sub_SCT_hist(subbasename);
%      fprintf('%d/%d SCT done\n',jj,numOfSubjectsGrp2);
% end

for jj = 1 : numOfSubjectsGrp1
  mat_fname=['E:\sipi_data\oasis_full_data\data\',fnamesGrp1{jj},'\PROCESSED\MPRAGE\T88_111\',fnamesGrp1{jj},'_aaj.SCT_hist.mat'];
  if ~exist(mat_fname,'file')
    continue;
  end
  mat_fname
  count = count + 1
  tmp = load (mat_fname, 'SCT_hist_atlas_left');
  SCT_hist_grp1 (:,:,count) = single (tmp.SCT_hist_atlas_left (:,indx_reduced));
end
SCT_hist_grp1 = SCT_hist_grp1 (:,:,1:count);
SCT_hist_grp1 = permute (SCT_hist_grp1, [1 3 2]);

% load SCT histograms for group 2 (moderate AD, label 1 in OASIS)
numOfSubjectsGrp2 = sum (CDRs == 1);
fnamesGrp2 = fnames (CDRs == 1);
SCT_hist_grp2 = zeros (numOfHistBins, numOfMeshNodes, numOfSubjectsGrp2, 'single');
count = 0;
for jj = 1 : numOfSubjectsGrp2
  mat_fname=['E:\sipi_data\oasis_full_data\data\',fnamesGrp2{jj},'\PROCESSED\MPRAGE\T88_111\',fnamesGrp2{jj},'_aaj.SCT_hist.mat'];
  if ~exist(mat_fname,'file')
    continue;
  end
  mat_fname
  count = count + 1
  tmp = load (mat_fname, 'SCT_hist_atlas_left');
  SCT_hist_grp2 (:,:,count) = single (tmp.SCT_hist_atlas_left (:,indx_reduced));
end
SCT_hist_grp2 = SCT_hist_grp2 (:,:,1:count);
SCT_hist_grp2 = permute (SCT_hist_grp2, [1 3 2]);

save ('data_oasis_sct_hist2_left.mat', 'SCT_hist_grp1', 'SCT_hist_grp2', 'indx_reduced','tar');

return
%%
%load ('data_oasis_sct_hist2_left.mat');
load ('data_oasis_sct_hist2_right.mat');
%SCT_hist_grp2=SCT_hist_grp2(:,19:end,:);

%tar_orig=readdfs('C:\Users\ajoshi\Documents\coding_ground\svreg-matlab/BrainSuiteAtlas1\mri.left.mid.cortex.dfs');
%tar=reducepatch(tar_orig,10000);
tarsm=smooth_cortex_fast(tar,.1,200);
cmap=cbrewer('seq','Oranges',100);

% SCT: call permutation testing for histograms
numOfPermutations = 100
for nb=1:1 % 20
  nb
  if nb==1
    [ t2Unpermuted t2permuted_all ] = permutation_test_hist (SCT_hist_grp1, SCT_hist_grp2, numOfPermutations);
  else
    [ t2Unpermuted t2permuted_all ] = permutation_test_hist (datasample(SCT_hist_grp1,size(SCT_hist_grp1,2),2), datasample(SCT_hist_grp2,size(SCT_hist_grp2,2),2), numOfPermutations);        
  end
  outlierRemovalPrctile = 97
  [ pValUnpermuted t2PermutedExtreme ] = permutation_test_p_values (t2Unpermuted, t2permuted_all, outlierRemovalPrctile);
  %
  pValUnpermuted_boot(:,nb)=pValUnpermuted';
  %
  t2Unpermuted_boot (:,nb) = t2Unpermuted;
  t2permuted_all_boot (:,:,nb) = t2permuted_all;
  %
  fprintf('nb=%d/20',nb);    
  'save start ...'
  save results_oasis_hist_kernel_conc1p5_dim20_perm100nb1_R.mat
  'save done'
end

show_surf_val(tarsm,(0.1-pValUnpermuted_boot(:,1)').*(pValUnpermuted_boot(:,1)'<=0.1),[0,0.1],cmap); axis equal;view(90,0); camlight;

show_surf_val(tarsm,var(pValUnpermuted_boot(:,[2:end]),[],2)',[0,0.3],cmap); axis equal;view(90,0); camlight;

figure (1); export_fig -a1 -m2 -transparent Pics/Hist_RiemannianKPCAdim20_SCT_R.png, pause (1)
figure (1); export_fig -a1 -m2 -transparent Pics/Var_Hist_RiemannianKPCAdim20_SCT_R.png

return
%%
%patch('vertices',tarsm.vertices,'faces',tar.faces,'facevertexcdata',pValUnpermuted','edgecolor','none','facecolor','interp');axis equal off tight;material dull;lighting phong; colormap jet;view(90,0);camlight;
show_surf_val(tarsm,(var(pValUnpermuted,1)),[0,.3],cmap);
export_fig -a1 -m2 -transparent Pics/bootstrap_SCT_right.png, close
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% S histogram only
S_hist_grp1 = squeeze (sum (sum (reshape (SCT_hist_grp1, 5, 5, 5, size(SCT_hist_grp1,2), size(SCT_hist_grp1,3)), 2), 3));
S_hist_grp2 = squeeze (sum (sum (reshape (SCT_hist_grp2, 5, 5, 5, size(SCT_hist_grp2,2), size(SCT_hist_grp2,3)), 2), 3));
[ S_t2Unpermuted S_t2permuted_all ] = permutation_test_hist (S_hist_grp1, S_hist_grp2, numOfPermutations);
outlierRemovalPrctile = 97
[ S_pValUnpermuted S_t2PermutedExtreme ] = permutation_test_p_values (S_t2Unpermuted, S_t2permuted_all, outlierRemovalPrctile);
show_surf_val(tarsm,1-S_pValUnpermuted,[0,1],cmap)
export_fig -a1 -m2 -transparent Pics/results_oasis_S_left.png, close

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% C histogram only
C_hist_grp1 = squeeze (sum (sum (reshape (SCT_hist_grp1, 5, 5, 5, size(SCT_hist_grp1,2), size(SCT_hist_grp1,3)), 1), 3));
C_hist_grp2 = squeeze (sum (sum (reshape (SCT_hist_grp2, 5, 5, 5, size(SCT_hist_grp2,2), size(SCT_hist_grp2,3)), 1), 3));
[ C_t2Unpermuted C_t2permuted_all ] = permutation_test_hist (C_hist_grp1, C_hist_grp2, numOfPermutations);
outlierRemovalPrctile = 97
[ C_pValUnpermuted C_t2PermutedExtreme ] = permutation_test_p_values (C_t2Unpermuted, C_t2permuted_all, outlierRemovalPrctile);
show_surf_val(tarsm,1-C_pValUnpermuted,[0,1],cmap)
export_fig -a1 -m2 -transparent Pics/results_oasis_C_left.png, close

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% T histogram only
T_hist_grp1 = squeeze (sum (sum (reshape (SCT_hist_grp1, 5, 5, 5, size(SCT_hist_grp1,2), size(SCT_hist_grp1,3)), 1), 2));
T_hist_grp2 = squeeze (sum (sum (reshape (SCT_hist_grp2, 5, 5, 5, size(SCT_hist_grp2,2), size(SCT_hist_grp2,3)), 1), 2));
[ T_t2Unpermuted T_t2permuted_all ] = permutation_test_hist (T_hist_grp1, T_hist_grp2, numOfPermutations);
outlierRemovalPrctile = 97
[ T_pValUnpermuted T_t2PermutedExtreme ] = permutation_test_p_values (T_t2Unpermuted, T_t2permuted_all, outlierRemovalPrctile);
show_surf_val(tarsm,1-T_pValUnpermuted,[0,1],cmap)
export_fig -a1 -m2 -transparent Pics/results_oasis_T_left.png, close

%[fp,tp]=roi_pval(pValSCT,(tar_orig.labels(indx_reduced)==226),1000);
%show_surf_val(tarsm,(tar_orig.labels(indx_reduced)==226),[0,1],cmap)
return
