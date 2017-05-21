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
%     sub_SCT(subbasename,NLevels);
%     fprintf('%d/%d SCT done\n',jj,length(aa1));
% end
% % generate SCT histograms for group 2, save them
% parfor jj=1:length(aa2)
%     subbasename=sprintf('%s/%s/anat/mprage_anonymized',dirname2,aa2(jj).name);
%     sub_SCT(subbasename,NLevels);
%     fprintf('%d/%d SCT done\n',jj,length(aa1));
% end

[N,T,R]=xlsread('E:\sipi_data\oasis_full_data\data\oasis_cross-sectional.csv');
fnames= R(2:end,1);
temp  = R(2:end,8);
for jj=1:length(fnames)
    CDRs(jj) = temp{jj};
end

% get size of each histogram to preallocate
mat_fname=['E:\sipi_data\oasis_full_data\data\',fnames{1},'\PROCESSED\MPRAGE\T88_111\',fnames{1},'_aaj.SCT.mat'];
tmp = load (mat_fname);
numLevels = size (tmp.SCTatlasright, 2);
numOfMeshNodes = length(tar.vertices);%size (tmp.SCTatlas_right, 2)

% load SCT histograms for group 1 (controls, label 0 in OASIS)
numOfSubjectsGrp1 = sum (CDRs == 0);
fnamesGrp1 = fnames (CDRs == 0);
SCTgrp1 = zeros (numLevels, numOfMeshNodes, 50, 'single');
SCTgrp2 = zeros (numLevels, numOfMeshNodes, numOfSubjectsGrp1-50, 'single');

count1 = 0; count2=0;
for jj = 1 : numOfSubjectsGrp1

    if jj<=50
        mat_fname=['E:\sipi_data\oasis_full_data\data\',fnamesGrp1{jj},'\PROCESSED\MPRAGE\T88_111\',fnamesGrp1{jj},'_aaj.SCT_simulation_0.1.mat'];
    else
        mat_fname=['E:\sipi_data\oasis_full_data\data\',fnamesGrp1{jj},'\PROCESSED\MPRAGE\T88_111\',fnamesGrp1{jj},'_aaj.SCT.mat'];
%        delete(['E:\sipi_data\oasis_full_data\data\',fnamesGrp1{jj},'\PROCESSED\MPRAGE\T88_111\',fnamesGrp1{jj},'_aaj.SCTsimulation.mat']);
    end
    if ~exist(mat_fname,'file')
        continue;
    end
    mat_fname

    tmp = load (mat_fname, 'SCTatlasright');
    if jj<=50
        count1 = count1 + 1;
        SCTgrp1 (:,:,count1) = single (tmp.SCTatlasright (indx_reduced,:))';
    else
        count2 = count2 + 1;        
        SCTgrp2 (:,:,count2) = single (tmp.SCTatlasright (indx_reduced,:))';
    end
end
SCTgrp1 = SCTgrp1 (:,:,1:count1);
SCTgrp1 = permute (SCTgrp1, [1 3 2]);
SCTgrp2 = SCTgrp2 (:,:,1:count2);
SCTgrp2 = permute (SCTgrp2, [1 3 2]);

save ('data_oasis_sct_simulation_0.1_right.mat', 'SCTgrp1', 'SCTgrp2', 'indx_reduced','tar');

meandiff = squeeze (sqrt (sum ( (mean(SCTgrp1(11:15,:,:),2) - mean(SCTgrp2(11:15,:,:),2)).^2 )));
tarsm=smooth_cortex_fast(tar,.5,200);
patch('vertices',tarsm.vertices,'faces',tar.faces,'facevertexcdata',meandiff,'edgecolor','none','facecolor','interp');axis equal;axis off;camlight;material dull;lighting phong;colormap jet;view(90,0);
export_fig -a1 -m2 Pics/oasis_sim_right_t.png; close

return

clear all;
close all;
load ('data_oasis_sct_simulation_0.1_right.mat');

tar_orig=readdfs('C:\Users\ajoshi\Documents\coding_ground\svreg-matlab/BrainSuiteAtlas1\mri.right.mid.cortex.dfs');
grnd_truth=double(tar_orig.labels==226);grnd_truth=grnd_truth(indx_reduced);
tar=reducepatch(tar_orig,10000);
tarsm=smooth_cortex_fast(tar,.1,200);

numOfPermutations = 10
cmap=cbrewer('seq','Oranges',100);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SCT histogram
% compute permutation distribution of test statistic
% [ t2unpermuted t2permuted_all ] = permutation_test (SCTgrp1, SCTgrp2, numOfPermutations);
[ t2unpermuted t2permuted_all ] = permutation_test(SCTgrp1(6:15,:,:), SCTgrp2(6:15,:,:), numOfPermutations);
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
[false_pos_rate,true_pos_rate]=roc_pval(pValUnpermuted,grnd_truth,1000);
export_fig -a1 -m2 -transparent Pics/results_oasis_sim_SCT_right.png, close
figure;plot(false_pos_rate,true_pos_rate)
export_fig -a1 -m2 -transparent Pics/ROC_sim_SCT_right.png, close

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% S histogram only
% show_surf_val(tarsm,S_pValUnpermuted,[0,1]);
Sgrp1 = SCTgrp1(1:5,:,:);% squeeze (sum (sum (reshape (SCTgrp1, 5, 5, 5, size(SCTgrp1,2), size(SCTgrp1,3)), 2), 3));
Sgrp2 = SCTgrp2(1:5,:,:);%squeeze (sum (sum (reshape (SCTgrp2, 5, 5, 5, size(SCTgrp2,2), size(SCTgrp2,3)), 2), 3));

% Sgrp1_noisy = Sgrp1 + 0.0 * rand (size (Sgrp1));
% Sgrp1_noisy = reshape (Sgrp1_noisy, size(Sgrp1,1), size(Sgrp1,2)*size(Sgrp1,3));
% Sgrp1_noisy = bsxfun (@times, Sgrp1_noisy, 1./squeeze(sum (Sgrp1_noisy,1)));
% Sgrp1_noisy = reshape (Sgrp1_noisy, size(Sgrp1,1), size(Sgrp1,2), size(Sgrp1,3));
% 
% Sgrp2_noisy = Sgrp2 + 0.0 * rand (size (Sgrp2));
% Sgrp2_noisy = reshape (Sgrp2_noisy, size(Sgrp2,1), size(Sgrp2,2)*size(Sgrp2,3));
% Sgrp2_noisy = bsxfun (@times, Sgrp2_noisy, 1./squeeze(sum (Sgrp2_noisy,1)));
% Sgrp2_noisy = reshape (Sgrp2_noisy, size(Sgrp2,1), size(Sgrp2,2), size(Sgrp2,3));

%numOfSubSelected = 100
% [ S_t2unpermuted S_t2permuted_all ] = permutation_test (Sgrp1_noisy(:,1:min(end,numOfSubSelected),:), Sgrp2_noisy(:,1:min(end,numOfSubSelected),:), numOfPermutations);
[ S_t2unpermuted S_t2permuted_all ] = permutation_test(Sgrp1, Sgrp2, numOfPermutations);
%[ S_t2unpermuted S_t2permuted_all ] = permutation_test (Sgrp1, Sgrp2, numOfPermutations);
outlierRemovalPrctile = 97
[ S_pValUnpermuted S_t2PermutedExtreme ] = permutation_test_p_values (S_t2unpermuted, S_t2permuted_all, outlierRemovalPrctile);
show_surf_val(tarsm,S_pValUnpermuted,[0,1],cmap); 
export_fig -a1 -m2 -transparent Pics/results_oasis_sim_S_right.png, close
[false_pos_rate,true_pos_rate]=roc_pval(S_pValUnpermuted,grnd_truth,1000);
figure;plot(false_pos_rate,true_pos_rate)
export_fig -a1 -m2 -transparent Pics/ROC_sim_S_right.png, close

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% C histogram only
Cgrp1 = SCTgrp1(6:10,:,:);%squeeze (sum (sum (reshape (SCTgrp1, 5, 5, 5, size(SCTgrp1,2), size(SCTgrp1,3)), 1), 3));
Cgrp2 = SCTgrp2(6:10,:,:);%squeeze (sum (sum (reshape (SCTgrp2, 5, 5, 5, size(SCTgrp2,2), size(SCTgrp2,3)), 1), 3));
% [ C_t2unpermuted C_t2permuted_all ] = permutation_test (Cgrp1, Cgrp2, numOfPermutations);
[ C_t2unpermuted C_t2permuted_all ] = permutation_test (Cgrp1, Cgrp2, numOfPermutations);
outlierRemovalPrctile = 97
[ C_pValUnpermuted C_t2PermutedExtreme ] = permutation_test_p_values (C_t2unpermuted, C_t2permuted_all, outlierRemovalPrctile);
show_surf_val(tarsm,1-C_pValUnpermuted,[0,1],cmap); 
% show_surf_val(tarsm,C_pValUnpermuted,[0,1]);
export_fig -a1 -m2 -transparent Pics/results_oasis_sim_C_right.png, close
[false_pos_rate,true_pos_rate]=roc_pval(C_pValUnpermuted,grnd_truth,1000);
figure;plot(false_pos_rate,true_pos_rate)
export_fig -a1 -m2 -transparent Pics/ROC_sim_C_right.png, close

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% T histogram only
Tgrp1 = SCTgrp1(11:15,:,:);%squeeze (sum (sum (reshape (SCTgrp1, 5, 5, 5, size(SCTgrp1,2), size(SCTgrp1,3)), 1), 2));
Tgrp2 = SCTgrp2(11:15,:,:);%squeeze (sum (sum (reshape (SCTgrp2, 5, 5, 5, size(SCTgrp2,2), size(SCTgrp2,3)), 1), 2));
[ T_t2unpermuted T_t2permuted_all ] = permutation_test (Tgrp1, Tgrp2, numOfPermutations);
%[ T_t2unpermuted T_t2permuted_all ] = permutation_test (Tgrp1, Tgrp2, numOfPermutations);
outlierRemovalPrctile = 97
[ T_pValUnpermuted T_t2PermutedExtreme ] = permutation_test_p_values (T_t2unpermuted, T_t2permuted_all, outlierRemovalPrctile);
show_surf_val(tarsm,1-T_pValUnpermuted,[0,1],cmap); 
% show_surf_val(tarsm,T_pValUnpermuted,[0,1]);
export_fig -a1 -m2-transparent  Pics/results_oasis_sim_T_right.png, close
[false_pos_rate,true_pos_rate]=roc_pval(T_pValUnpermuted,grnd_truth,1000);
figure;plot(false_pos_rate,true_pos_rate)
export_fig -a1 -m2 -transparent Pics/ROC_sim_T_right.png, close

% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % SC histogram only
% SCgrp1 = squeeze (sum (reshape (SCTgrp1, 5, 5, 5, size(SCTgrp1,2), size(SCTgrp1,3)), 3));
% SCgrp1 = reshape (SCgrp1, 25, size(SCTgrp1,2), size(SCTgrp1,3));
% SCgrp2 = squeeze (sum (reshape (SCTgrp2, 5, 5, 5, size(SCTgrp2,2), size(SCTgrp2,3)), 3));
% SCgrp2 = reshape (SCgrp2, 25, size(SCTgrp2,2), size(SCTgrp2,3));
% % [ SC_t2unpermuted SC_t2permuted_all ] = permutation_test (SCgrp1, SCgrp2, numOfPermutations);
% [ SC_t2unpermuted SC_t2permuted_all ] = permutation_test(SCgrp1, SCgrp2, numOfPermutations);
% outlierRemovalPrctile = 97
% [ SC_pValUnpermuted SC_t2PermutedExtreme ] = permutation_test_p_values (SC_t2unpermuted, SC_t2permuted_all, outlierRemovalPrctile);
% show_surf_val(tarsm,1-SC_pValUnpermuted,[0,1],cmap); 
% % show_surf_val(tarsm,SC_pValUnpermuted,[0,1]);
% export_fig -a1 -m2 -transparent Pics/results_oasis_sim_SC_right.png
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % ST histogram only
% STgrp1 = squeeze (sum (reshape (SCTgrp1, 5, 5, 5, size(SCTgrp1,2), size(SCTgrp1,3)), 2));
% STgrp1 = reshape (STgrp1, 25, size(SCTgrp1,2), size(SCTgrp1,3));
% STgrp2 = squeeze (sum (reshape (SCTgrp2, 5, 5, 5, size(SCTgrp2,2), size(SCTgrp2,3)), 2));
% STgrp2 = reshape (STgrp2, 25, size(SCTgrp2,2), size(SCTgrp2,3));
% % [ ST_t2unpermuted ST_t2permuted_all ] = permutation_test (STgrp1, STgrp2, numOfPermutations);
% [ ST_t2unpermuted ST_t2permuted_all ] = permutation_test(STgrp1, STgrp2, numOfPermutations);
% outlierRemovalPrctile = 97
% [ ST_pValUnpermuted ST_t2PermutedExtreme ] = permutation_test_p_values (ST_t2unpermuted, ST_t2permuted_all, outlierRemovalPrctile);
% show_surf_val(tarsm,1-ST_pValUnpermuted,[0,1],cmap); 
% % show_surf_val(tarsm,ST_pValUnpermuted,[0,1]);
% export_fig -a1 -m2-transparent Pics/results_oasis_sim_ST_right.png
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % CT histogram only
% CTgrp1 = squeeze (sum (reshape (SCTgrp1, 5, 5, 5, size(SCTgrp1,2), size(SCTgrp1,3)), 1));
% CTgrp1 = reshape (CTgrp1, 25, size(SCTgrp1,2), size(SCTgrp1,3));
% CTgrp2 = squeeze (sum (reshape (SCTgrp2, 5, 5, 5, size(SCTgrp2,2), size(SCTgrp2,3)), 1));
% CTgrp2 = reshape (CTgrp2, 25, size(SCTgrp2,2), size(SCTgrp2,3));
% % [ CT_t2unpermuted CT_t2permuted_all ] = permutation_test (CTgrp1, CTgrp2, numOfPermutations);
% [ CT_t2unpermuted CT_t2permuted_all ] = permutation_test(CTgrp1, CTgrp2, numOfPermutations);
% outlierRemovalPrctile = 97
% [ CT_pValUnpermuted CT_t2PermutedExtreme ] = permutation_test_p_values (CT_t2unpermuted, CT_t2permuted_all, outlierRemovalPrctile);
% show_surf_val(tarsm,1-CT_pValUnpermuted,[0,1],cmap); 
% % show_surf_val(tarsm,CT_pValUnpermuted,[0,1]);
% export_fig -a1 -m2 -transparent Pics/results_oasis_sim_CT_right.png, close
% 
% save ( 'results_simulation_0.1_right.mat');
% 
% return