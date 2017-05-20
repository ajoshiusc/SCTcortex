
opengl software

clc;clear ;close all;restoredefaultpath;
addpath(genpath('/home/ajoshi/git_sandbox/svreg-matlab/3rdParty'))
addpath(genpath('/home/ajoshi/git_sandbox/svreg-matlab/MEX_Files'))
addpath(genpath('/home/ajoshi/git_sandbox/svreg-matlab/src'))
addpath(genpath('/home/ajoshi/git_sandbox/svreg-matlab/dev'));

%addpath(genpath('/home/ajoshi/git_sandbox/svreg-matlab\dev\cortex_analysis'));

lst1=dir('/home/ajoshi/HCP900/data/*');

% build histograms at each voxel
count=0;
 parfor jj=1:length(lst1)
     fname=['/home/ajoshi/HCP900/data/',lst1(jj).name,'/colin_brain/',lst1(jj).name];
     try
     if exist([fname,'.SCT_hist.mat'],'file')
         sub_SCT(fname);
     end
     catch
         fprintf(fname);
     end
 end
%%
% read all hists
subno=1;
load('/home/ajoshi/HCP900/xls_hcp900.mat');
%[num,txt,raw]=xlsread('/home/ajoshi/HCP900/hcp_unrestricted_aajoshi_4_17_2016_14_56_18.csv');

npts=1000;
tar_orig=readdfs('/home/ajoshi/git_sandbox/svreg-matlab/BrainSuiteAtlas1/mri.left.mid.cortex.dfs');
tar=reducepatch(tar_orig,10000);
[~,ia,indx_reduced]=intersect(tar.vertices,tar_orig.vertices,'rows','stable');

for jj=1:length(num(:,1))
    subname=['/home/ajoshi/HCP900/data/',num2str(num(jj,1)),'/colin_brain/',num2str(num(jj,1))];
    if isnan(num(jj,182)) % 
        continue
    end
    try
        tmp=load([subname,'.SCT.mat']);
        %SCT_atlas_left_all (:,:,subno)=single(tmp.SCTatlasleft (indx_reduced,:)');
        SCT_atlas_right_all(:,:,subno)=single(tmp.SCTatlasright(indx_reduced,:)');
        cogScore (subno)=num(jj,182);
        subno
        subno=subno+1;
    catch
        % subname
    end
end
SCT_atlas_left_all = permute (SCT_atlas_left_all, [1 3 2]);
%SCT_hist_atlas_right_all = permute (SCT_atlas_right_all, [1 3 2]);
% now dimensions of SCT_hist_atlas_left_all, SCT_hist_atlas_right_all are
% (histLength x numOfsubjects x numOfVoxels)
% save sct_measures_left_all
save sct_measures_right_all

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
opengl software
%
addpath(genpath('/home/ajoshi/git_sandbox/svreg-matlab/src'))
addpath(genpath('/home/ajoshi/git_sandbox/svreg-matlab/dev'))
%
load ('sct_measures_left_all.mat', 'SCT_atlas_left_all', 'cogScore', 'tar');
% load ('sct_measures_right_all.mat', 'SCT_atlas_right_all', 'cogScore', 'tar');
%
SCT_all = SCT_atlas_left_all;
% get similarity matrix
numOfSubjects = size (SCT_all, 2);
numOfVertices = size (SCT_all, 3);
SCTgram_all = nan (numOfSubjects, numOfSubjects, numOfVertices, 'single');
for jj = 1 : numOfVertices
    SCTgram_all (:,:,jj) = exp (- 0.5 * pdist2 (SCT_all(:,:,jj)', SCT_all(:,:,jj)').^2);
end
%SCT_atlas_all = [];
SCTgram_all = [];
clear SCT_atlas_left_all SCT_atlas_right_all
whos
%save sct_histsSqrtCrossExp_left_all % do NOT save this. Take way too much time.

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%load sct_histsSqrtCrossExp_left_all
% permutation testing - compute null distribution of extremal test statistic
% matlabpool close force local
% matlabpool 6
numOfPermutations = 20
for nb = 1 : 20
    '--------------------------------------------------------', nb
    % create a bootstrap sample
    if nb == 1
        bootstrapSampleIndices = [1:numOfSubjects];
    else
        rng (100 + nb);
        [ bootstrapSampleDummy bootstrapSampleIndices ]  = datasample ([1:numOfSubjects], numOfSubjects);
        bootstrapSampleIndices = sort (bootstrapSampleIndices); % sort to reduce data movement in next step
    end
    %
    [ testStatsUnpermuted_all, ...
        testStatsPermuted_Ftest_all, ...
        testStatsPermuted_Bartlett_all, ...
        testStatsPermuted_BrownForsythe_all, ...
        testStatsPermuted_LeveneAbsolute_all, ...
        testStatsPermuted_AnsariBradley_all ] ...
        = permutation_test_hist_regression_all ...
        (SCTgram_all (bootstrapSampleIndices, bootstrapSampleIndices, :), ...
        SCT_all, ... % (:, bootstrapSampleIndices, :), ...
        cogScore (bootstrapSampleIndices), ...
        numOfPermutations);
    %
    testStatsUnpermuted_all_boot (:,nb)                = testStatsUnpermuted_all;
    testStatsPermuted_Ftest_all_boot (:,:,nb)          = testStatsPermuted_Ftest_all;
    testStatsPermuted_Bartlett_all_boot (:,:,nb)       = testStatsPermuted_Bartlett_all;
    testStatsPermuted_BrownForsythe_all_boot (:,:,nb)  = testStatsPermuted_BrownForsythe_all;
    testStatsPermuted_LeveneAbsolute_all_boot (:,:,nb) = testStatsPermuted_LeveneAbsolute_all;
    testStatsPermuted_AnsariBradley_all_boot (:,:,nb)  = testStatsPermuted_AnsariBradley_all;
    %
    save ('sct_testStats_left_all_nonlinear.mat', ...
        'testStatsUnpermuted_all_boot', ...
        'testStatsPermuted_Ftest_all_boot', ...
        'testStatsPermuted_Bartlett_all_boot', ...
        'testStatsPermuted_BrownForsythe_all_boot', ...
        'testStatsPermuted_LeveneAbsolute_all_boot', ...
        'testStatsPermuted_AnsariBradley_all_boot' );
end

%% permutation testing -compute p values
outlierRemovalPrctile = 90
%
for nb = 1 : size (testStatsPermuted_Ftest_all_boot, 3)
    [ pValUnpermuted_Ftest_all, testStatsPermutedExtreme_Ftest_all ] ...
        = permutation_test_p_values (testStatsPermuted_Ftest_all_boot(:,1,nb), testStatsPermuted_Ftest_all_boot(:,:,nb), outlierRemovalPrctile);
    pValUnpermuted_Ftest_all_boot (:,nb) = pValUnpermuted_Ftest_all;
end

% visualize unpermuted p values
%right_hemi_surf=readdfs('/home/ajoshi/git_sandbox/svreg-matlab/BrainSuiteAtlas1/mri.right.mid.cortex.dfs');
right_hemi_surf=readdfs('/home/ajoshi/git_sandbox/svreg-matlab/BrainSuiteAtlas1/mri.left.mid.cortex.dfs');
right_hemi_surf=smooth_cortex_fast(tar,.5,300);
tarsm=smooth_cortex_fast(tar,.1,200);
cmap=cbrewer('seq','Oranges',100);
%
pp=(pValUnpermuted_Ftest_all_boot(:,1)<0.1).*(0.1-pValUnpermuted_Ftest_all_boot(:,1));
pv=smooth_surf_function(right_hemi_surf,pp,1,1);

%
show_surf_val(tarsm, pv, [0,0.05], cmap); view(-90,0); camlight; colorbar off;
pause, export_fig -a1 -m2 -transparent Pics/HCPregress_Nonlinear_SCT_L.png, close
show_surf_val(tarsm, var(pValUnpermuted_Ftest_all_boot,[],2), [0,0.3], cmap); view(-90,0); camlight; colorbar off;
pause, export_fig -a1 -m2 -transparent Pics/HCPregress_Var_Nonlinear_SCT_L.png, close

%% permutation testing -compute p values
outlierRemovalPrctile = 90
%
[ pValUnpermuted_Ftest_all, testStatsPermutedExtreme_Ftest_all ] ...
    = permutation_test_p_values (testStatsPermuted_Ftest_all(:,1), testStatsPermuted_Ftest_all (:,1:11), outlierRemovalPrctile);
hist (pValUnpermuted_Ftest_all); pause, close
hist (testStatsPermutedExtreme_Ftest_all, 30); pause, close
hist (testStatsPermuted_Ftest_all (:,1), 30); pause, close

%
[ pValUnpermuted_Bartlett_all, testStatsPermutedExtreme_Bartlett_all ] ...
    = permutation_test_p_values (testStatsPermuted_Bartlett_all(:,1), testStatsPermuted_Bartlett_all (:,1:350), outlierRemovalPrctile);
hist (pValUnpermuted_Bartlett_all); pause, close
%
[ pValUnpermuted_BrownForsythe_all, testStatsPermutedExtreme_BrownForsythe_all ] ...
    = permutation_test_p_values (testStatsPermuted_BrownForsythe_all(:,1), testStatsPermuted_BrownForsythe_all (:,1:350), outlierRemovalPrctile);
hist (pValUnpermuted_BrownForsythe_all); pause, close
%
[ pValUnpermuted_LeveneAbsolute_all, testStatsPermutedExtreme_LeveneAbsolute_all ] ...
    = permutation_test_p_values (testStatsPermuted_LeveneAbsolute_all(:,1), testStatsPermuted_LeveneAbsolute_all (:,1:350), outlierRemovalPrctile);
hist (pValUnpermuted_LeveneAbsolute_all); pause, close
%
[ pValUnpermuted_AnsariBradley_all, testStatsPermutedExtreme_AnsariBradley_all ] ...
    = permutation_test_p_values (testStatsPermuted_AnsariBradley_all(:,1), testStatsPermuted_AnsariBradley_all (:,1:350), outlierRemovalPrctile);
hist (pValUnpermuted_AnsariBradley_all); pause, close



% visualize of std dev in prediction error for all vertices
left_hemi_surf=readdfs('/home/ajoshi/git_sandbox/svreg-matlab/BrainSuiteAtlas1/mri.left.mid.cortex.dfs');
left_hemi_surf=smooth_cortex_fast(tar,.5,300);
figure;
%patch('faces',left_hemi_surf.faces,'vertices',left_hemi_surf.vertices,'facevertexcdata',mean(fit.predictionErrors_all)','edgecolor','none','facecolor','interp');
patch('faces',left_hemi_surf.faces,'vertices',left_hemi_surf.vertices,'facevertexcdata',std(fit.predictionErrors_all)','edgecolor','none','facecolor','interp');
axis equal;axis off; view (-90,0); camlight; colorbar; material dull;
pause, close

% visualize RMSE
left_hemi_surf=readdfs('/home/ajoshi/git_sandbox/svreg-matlab/BrainSuiteAtlas1/mri.left.mid.cortex.dfs');
left_hemi_surf=smooth_cortex_fast(tar,.5,300);
figure;
patch('faces',left_hemi_surf.faces,'vertices',left_hemi_surf.vertices,'facevertexcdata',fit.predictionRMSE_all,'edgecolor','none','facecolor','interp');
axis equal;axis off; view (-90,0); camlight; colorbar; material dull;
pause, close

% visualize conc
left_hemi_surf=readdfs('/home/ajoshi/git_sandbox/svreg-matlab/BrainSuiteAtlas1/mri.left.mid.cortex.dfs');
left_hemi_surf=smooth_cortex_fast(tar,.5,300);
figure;
patch('faces',left_hemi_surf.faces,'vertices',left_hemi_surf.vertices,'facevertexcdata',fit.conc_all,'edgecolor','none','facecolor','interp');
axis equal;axis off; view (-90,0); camlight; colorbar; material dull;
pause, close

% visualize p values
left_hemi_surf=smooth_cortex_fast(tar,.5,300);
figure;
%patch('faces',left_hemi_surf.faces,'vertices',left_hemi_surf.vertices,'facevertexcdata',testStats.p_Ftest_all','edgecolor','none','facecolor','interp');
%patch('faces',left_hemi_surf.faces,'vertices',left_hemi_surf.vertices,'facevertexcdata',testStats.p_Bartlett_all','edgecolor','none','facecolor','interp');
%patch('faces',left_hemi_surf.faces,'vertices',left_hemi_surf.vertices,'facevertexcdata',testStats.p_BrownForsythe_all','edgecolor','none','facecolor','interp');
%patch('faces',left_hemi_surf.faces,'vertices',left_hemi_surf.vertices,'facevertexcdata',testStats.p_LeveneAbsolute_all','edgecolor','none','facecolor','interp');
pfdr=FDR(p_AnsariBradley_all',0.05);pfdr=0.5
patch('faces',left_hemi_surf.faces,'vertices',left_hemi_surf.vertices,'facevertexcdata',(pfdr-testStats.p_AnsariBradley_all').*double(testStats.p_AnsariBradley_all'< pfdr),'edgecolor','none','facecolor','interp');
axis equal;axis off;camlight; colorbar; view (-90,0); camlight; material dull;%caxis([0,pfdr]);j1=flipud(jet);
pause, close
