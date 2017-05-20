clc;clear all;close all;

% This script performs pointwise group differences in cortical thicknesses
% for the right hemisphere
% modify as appropriate for right hemisphere

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
tar_orig=readdfs('C:\Users\ajoshi\Documents\coding_ground\svreg-matlab/BrainSuiteAtlas1\mri.right.mid.cortex.dfs');
tar=reducepatch(tar_orig,10000);
[~,ia,indx_reduced]=intersect(tar.vertices,tar_orig.vertices,'rows','stable');

% % generate SCT histograms for group 1, save them
% parfor jj=1:length(aa1)
%     subbasename=sprintf('%s/%s/%s_mpr_n4_anon_sbj_111.RAS',dirname1,aa1(jj).name,aa1(jj).name);
%     sub_SCT_hist(subbasename,NLevels);
%     fprintf('%d/%d SCT done\n',jj,length(aa1));
% end
% % generate SCT histograms for group 2, save them
% parfor jj=1:length(aa2)
%     subbasename=sprintf('%s/%s/anat/mprage_anonymized',dirname2,aa2(jj).name);
%     sub_SCT_hist(subbasename,NLevels);
%     fprintf('%d/%d SCT done\n',jj,length(aa1));
% end

[N,T,R]=xlsread('E:\sipi_data\oasis_full_data\data\oasis_cross-sectional.csv');
fnames= R(2:end,1);
temp  = R(2:end,8);
for jj=1:length(fnames)
    CDRs(jj) = temp{jj};
end

% get size of each histogram to preallocate
mat_fname=['E:\sipi_data\oasis_full_data\data\',fnames{1},'\PROCESSED\MPRAGE\T88_111\',fnames{1},'_aaj.SCT_hist.mat'];
tmp = load (mat_fname);
numOfHistBins = size (tmp.SCT_hist_atlas_right, 1)
numOfMeshNodes = length(tar.vertices);%size (tmp.SCT_hist_atlas_right, 2)

% load SCT histograms for group 1 (controls, label 0 in OASIS)
numOfSubjectsGrp1 = sum (CDRs == 0);
fnamesGrp1 = fnames (CDRs == 0);
SCT_hist_grp1 = zeros (numOfHistBins, numOfMeshNodes, 50, 'single');
SCT_hist_grp2 = zeros (numOfHistBins, numOfMeshNodes, numOfSubjectsGrp1-50, 'single');

count1 = 0; count2=0;
for jj = 1 : numOfSubjectsGrp1

    if jj<=50
        mat_fname=['E:\sipi_data\oasis_full_data\data\',fnamesGrp1{jj},'\PROCESSED\MPRAGE\T88_111\',fnamesGrp1{jj},'_aaj.SCT_hist_simulation_0.1.mat'];
    else
        mat_fname=['E:\sipi_data\oasis_full_data\data\',fnamesGrp1{jj},'\PROCESSED\MPRAGE\T88_111\',fnamesGrp1{jj},'_aaj.SCT_hist.mat'];
%        delete(['E:\sipi_data\oasis_full_data\data\',fnamesGrp1{jj},'\PROCESSED\MPRAGE\T88_111\',fnamesGrp1{jj},'_aaj.SCT_hist_simulation.mat']);
    end
    if ~exist(mat_fname,'file')
        continue;
    end
    mat_fname

    tmp = load (mat_fname, 'SCT_hist_atlas_right');
    if jj<=50
        count1 = count1 + 1;
        SCT_hist_grp1 (:,:,count1) = single (tmp.SCT_hist_atlas_right (:,indx_reduced));
    else
        count2 = count2 + 1;        
        SCT_hist_grp2 (:,:,count2) = single (tmp.SCT_hist_atlas_right (:,indx_reduced));
    end
end
SCT_hist_grp1 = SCT_hist_grp1 (:,:,1:count1);
SCT_hist_grp1 = permute (SCT_hist_grp1, [1 3 2]);
SCT_hist_grp2 = SCT_hist_grp2 (:,:,1:count2);
SCT_hist_grp2 = permute (SCT_hist_grp2, [1 3 2]);

save ('data_oasis_sct_hist_simulation_0.1_right.mat', 'SCT_hist_grp1', 'SCT_hist_grp2', 'indx_reduced','tar');

mean_hist_diff = squeeze (sqrt (sum ( (mean(SCT_hist_grp1,2) - mean(SCT_hist_grp2,2)).^2 )));
tarsm=smooth_cortex_fast(tar,.5,200);
patch('vertices',tarsm.vertices,'faces',tar.faces,'facevertexcdata',mean_hist_diff,'edgecolor','none','facecolor','interp');axis equal;axis off;camlight;material dull;lighting phong;colormap jet;view(90,0);
export_fig -a1 -m2 Pics/oasis_sim_right.png; close

return

clear all
close all
load ('data_oasis_sct_hist_simulation_0.1_right.mat');

tar_orig=readdfs('C:\Users\ajoshi\Documents\coding_ground\svreg-matlab/BrainSuiteAtlas1\mri.right.mid.cortex.dfs');
grnd_truth=double(tar_orig.labels==226);grnd_truth=grnd_truth(indx_reduced);
tar=reducepatch(tar_orig,10000);
tarsm=smooth_cortex_fast(tar,.1,200);

numOfPermutations = 10
cmap=cbrewer('seq','Oranges',100);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SCT histogram
% compute permutation distribution of test statistic
% [ t2unpermuted t2permuted_all ] = permutation_test_hist (SCT_hist_grp1, SCT_hist_grp2, numOfPermutations);
[ t2unpermuted t2permuted_all ] = permutation_test_hist(SCT_hist_grp1, SCT_hist_grp2, numOfPermutations);
% see test statistiCS
%show_surf_val(tarsm,t2unpermuted); pause, close
%show_surf_val(tarsm,t2permuted_all(:,end)); pause, close
% compute p values
outlierRemovalPrctile = 97
[ pValUnpermuted t2PermutedExtreme ] = permutation_test_p_values (t2unpermuted, t2permuted_all, outlierRemovalPrctile);
% see p values
show_surf_val(tarsm,(1-pValUnpermuted),[0,1],cmap);
% show_surf_val(tarsm,pValUnpermuted,[0,1]);
%show_surf_val(tarsm,pValUnpermuted,[0,1]);
%patch('vertices',tarsm.vertices,'faces',tar.faces,'facevertexcdata',pValUnpermuted','edgecolor','none','facecolor','interp');axis equal off tight;camlight;material dull;lighting phong; colormap jet;camlight;view(90,0);
export_fig -a1 -m2 -transparent Pics/results_oasis_sim_SCT_right_Riemannian.png, close
[false_pos_rate,true_pos_rate]=roc_pval(pValUnpermuted,grnd_truth,1000);
figure;plot(false_pos_rate,true_pos_rate)
export_fig -a1 -m2 -transparent Pics/ROC_sim_SCT_hist_right.png, close

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% S histogram only
% show_surf_val(tarsm,S_pValUnpermuted,[0,1]);
S_hist_grp1 = squeeze (sum (sum (reshape (SCT_hist_grp1, 5, 5, 5, size(SCT_hist_grp1,2), size(SCT_hist_grp1,3)), 2), 3));
S_hist_grp2 = squeeze (sum (sum (reshape (SCT_hist_grp2, 5, 5, 5, size(SCT_hist_grp2,2), size(SCT_hist_grp2,3)), 2), 3));

% S_hist_grp1_noisy = S_hist_grp1 + 0.0 * rand (size (S_hist_grp1));
% S_hist_grp1_noisy = reshape (S_hist_grp1_noisy, size(S_hist_grp1,1), size(S_hist_grp1,2)*size(S_hist_grp1,3));
% S_hist_grp1_noisy = bsxfun (@times, S_hist_grp1_noisy, 1./squeeze(sum (S_hist_grp1_noisy,1)));
% S_hist_grp1_noisy = reshape (S_hist_grp1_noisy, size(S_hist_grp1,1), size(S_hist_grp1,2), size(S_hist_grp1,3));
% 
% S_hist_grp2_noisy = S_hist_grp2 + 0.0 * rand (size (S_hist_grp2));
% S_hist_grp2_noisy = reshape (S_hist_grp2_noisy, size(S_hist_grp2,1), size(S_hist_grp2,2)*size(S_hist_grp2,3));
% S_hist_grp2_noisy = bsxfun (@times, S_hist_grp2_noisy, 1./squeeze(sum (S_hist_grp2_noisy,1)));
% S_hist_grp2_noisy = reshape (S_hist_grp2_noisy, size(S_hist_grp2,1), size(S_hist_grp2,2), size(S_hist_grp2,3));

%numOfSubSelected = 100
% [ S_t2unpermuted S_t2permuted_all ] = permutation_test_hist (S_hist_grp1_noisy(:,1:min(end,numOfSubSelected),:), S_hist_grp2_noisy(:,1:min(end,numOfSubSelected),:), numOfPermutations);
[ S_t2unpermuted S_t2permuted_all ] = permutation_test_hist(S_hist_grp1, S_hist_grp2, numOfPermutations);
%[ S_t2unpermuted S_t2permuted_all ] = permutation_test (S_hist_grp1, S_hist_grp2, numOfPermutations);
outlierRemovalPrctile = 97
[ S_pValUnpermuted S_t2PermutedExtreme ] = permutation_test_p_values (S_t2unpermuted, S_t2permuted_all, outlierRemovalPrctile);
show_surf_val(tarsm,1-S_pValUnpermuted,[0,1],cmap); 
export_fig -a1 -m2 -transparent Pics/results_oasis_sim_S_right_Riemannian.png, close
[false_pos_rate,true_pos_rate]=roc_pval(S_pValUnpermuted,grnd_truth,1000);
figure;plot(false_pos_rate,true_pos_rate)
export_fig -a1 -m2 -transparent Pics/ROC_sim_S_hist_right.png, close

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% C histogram only
C_hist_grp1 = squeeze (sum (sum (reshape (SCT_hist_grp1, 5, 5, 5, size(SCT_hist_grp1,2), size(SCT_hist_grp1,3)), 1), 3));
C_hist_grp2 = squeeze (sum (sum (reshape (SCT_hist_grp2, 5, 5, 5, size(SCT_hist_grp2,2), size(SCT_hist_grp2,3)), 1), 3));
% [ C_t2unpermuted C_t2permuted_all ] = permutation_test_hist (C_hist_grp1, C_hist_grp2, numOfPermutations);
[ C_t2unpermuted C_t2permuted_all ] = permutation_test_hist (C_hist_grp1, C_hist_grp2, numOfPermutations);
outlierRemovalPrctile = 97
[ C_pValUnpermuted C_t2PermutedExtreme ] = permutation_test_p_values (C_t2unpermuted, C_t2permuted_all, outlierRemovalPrctile);
show_surf_val(tarsm,1-C_pValUnpermuted,[0,1],cmap); 
% show_surf_val(tarsm,C_pValUnpermuted,[0,1]);
export_fig -a1 -m2 -transparent Pics/results_oasis_sim_C_right_Riemannian.png, close
pValC=C_pValUnpermuted;
[false_pos_rate,true_pos_rate]=roc_pval(C_pValUnpermuted,grnd_truth,1000);
figure;plot(false_pos_rate,true_pos_rate)
export_fig -a1 -m2 -transparent Pics/ROC_sim_C_hist_right.png, close

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% T histogram only
T_hist_grp1 = squeeze (sum (sum (reshape (SCT_hist_grp1, 5, 5, 5, size(SCT_hist_grp1,2), size(SCT_hist_grp1,3)), 1), 2));
T_hist_grp2 = squeeze (sum (sum (reshape (SCT_hist_grp2, 5, 5, 5, size(SCT_hist_grp2,2), size(SCT_hist_grp2,3)), 1), 2));
[ T_t2unpermuted T_t2permuted_all ] = permutation_test_hist (T_hist_grp1, T_hist_grp2, numOfPermutations);
%[ T_t2unpermuted T_t2permuted_all ] = permutation_test (T_hist_grp1, T_hist_grp2, numOfPermutations);
outlierRemovalPrctile = 97
[ T_pValUnpermuted T_t2PermutedExtreme ] = permutation_test_p_values (T_t2unpermuted, T_t2permuted_all, outlierRemovalPrctile);
show_surf_val(tarsm,1-T_pValUnpermuted,[0,1],cmap); 
% show_surf_val(tarsm,T_pValUnpermuted,[0,1]);
export_fig -a1 -m2-transparent  Pics/results_oasis_sim_T_right_Riemannian.png, close
pValT=T_pValUnpermuted;
[false_pos_rate,true_pos_rate]=roc_pval(T_pValUnpermuted,grnd_truth,1000);
figure;plot(false_pos_rate,true_pos_rate)
export_fig -a1 -m2 -transparent Pics/ROC_sim_T_hist_right.png, close

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % SC histogram only
% SC_hist_grp1 = squeeze (sum (reshape (SCT_hist_grp1, 5, 5, 5, size(SCT_hist_grp1,2), size(SCT_hist_grp1,3)), 3));
% SC_hist_grp1 = reshape (SC_hist_grp1, 25, size(SCT_hist_grp1,2), size(SCT_hist_grp1,3));
% SC_hist_grp2 = squeeze (sum (reshape (SCT_hist_grp2, 5, 5, 5, size(SCT_hist_grp2,2), size(SCT_hist_grp2,3)), 3));
% SC_hist_grp2 = reshape (SC_hist_grp2, 25, size(SCT_hist_grp2,2), size(SCT_hist_grp2,3));
% % [ SC_t2unpermuted SC_t2permuted_all ] = permutation_test_hist (SC_hist_grp1, SC_hist_grp2, numOfPermutations);
% [ SC_t2unpermuted SC_t2permuted_all ] = permutation_test_hist(SC_hist_grp1, SC_hist_grp2, numOfPermutations);
% outlierRemovalPrctile = 97
% [ SC_pValUnpermuted SC_t2PermutedExtreme ] = permutation_test_p_values (SC_t2unpermuted, SC_t2permuted_all, outlierRemovalPrctile);
% show_surf_val(tarsm,1-SC_pValUnpermuted,[0,1],cmap); 
% % show_surf_val(tarsm,SC_pValUnpermuted,[0,1]);
% export_fig -a1 -m2 -transparent Pics/results_oasis_sim_SC_right_Riemannian.png
% pValSC=SC_pValUnpermuted;
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % ST histogram only
% ST_hist_grp1 = squeeze (sum (reshape (SCT_hist_grp1, 5, 5, 5, size(SCT_hist_grp1,2), size(SCT_hist_grp1,3)), 2));
% ST_hist_grp1 = reshape (ST_hist_grp1, 25, size(SCT_hist_grp1,2), size(SCT_hist_grp1,3));
% ST_hist_grp2 = squeeze (sum (reshape (SCT_hist_grp2, 5, 5, 5, size(SCT_hist_grp2,2), size(SCT_hist_grp2,3)), 2));
% ST_hist_grp2 = reshape (ST_hist_grp2, 25, size(SCT_hist_grp2,2), size(SCT_hist_grp2,3));
% % [ ST_t2unpermuted ST_t2permuted_all ] = permutation_test_hist (ST_hist_grp1, ST_hist_grp2, numOfPermutations);
% [ ST_t2unpermuted ST_t2permuted_all ] = permutation_test_hist(ST_hist_grp1, ST_hist_grp2, numOfPermutations);
% outlierRemovalPrctile = 97
% [ ST_pValUnpermuted ST_t2PermutedExtreme ] = permutation_test_p_values (ST_t2unpermuted, ST_t2permuted_all, outlierRemovalPrctile);
% show_surf_val(tarsm,1-ST_pValUnpermuted,[0,1],cmap); 
% % show_surf_val(tarsm,ST_pValUnpermuted,[0,1]);
% export_fig -a1 -m2-transparent Pics/results_oasis_sim_ST_right_Riemannian.png
% pValST=ST_pValUnpermuted;
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % CT histogram only
% CT_hist_grp1 = squeeze (sum (reshape (SCT_hist_grp1, 5, 5, 5, size(SCT_hist_grp1,2), size(SCT_hist_grp1,3)), 1));
% CT_hist_grp1 = reshape (CT_hist_grp1, 25, size(SCT_hist_grp1,2), size(SCT_hist_grp1,3));
% CT_hist_grp2 = squeeze (sum (reshape (SCT_hist_grp2, 5, 5, 5, size(SCT_hist_grp2,2), size(SCT_hist_grp2,3)), 1));
% CT_hist_grp2 = reshape (CT_hist_grp2, 25, size(SCT_hist_grp2,2), size(SCT_hist_grp2,3));
% % [ CT_t2unpermuted CT_t2permuted_all ] = permutation_test_hist (CT_hist_grp1, CT_hist_grp2, numOfPermutations);
% [ CT_t2unpermuted CT_t2permuted_all ] = permutation_test_hist(CT_hist_grp1, CT_hist_grp2, numOfPermutations);
% outlierRemovalPrctile = 97
% [ CT_pValUnpermuted CT_t2PermutedExtreme ] = permutation_test_p_values (CT_t2unpermuted, CT_t2permuted_all, outlierRemovalPrctile);
% show_surf_val(tarsm,1-CT_pValUnpermuted,[0,1],cmap); 
% % show_surf_val(tarsm,CT_pValUnpermuted,[0,1]);
% export_fig -a1 -m2 -transparent Pics/results_oasis_sim_CT_right_Riemannian.png, close
% pValCT=CT_pValUnpermuted;
% 
% save ( 'results_simulation_0.1_right_Riemannian.mat');
% [fp,tp]=roi_pval(pValSCT,(tar_orig.labels(indx_reduced)==226),10000);hold on;plot(fp,tp,'bo-');
% show_surf_val(tarsm,1-pValSCT,[0,1],cmap)
% save('iter1000_SCT_hist_sim.mat')

return