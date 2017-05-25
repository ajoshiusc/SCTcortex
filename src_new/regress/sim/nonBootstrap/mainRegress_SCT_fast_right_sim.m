
opengl software

clc;clear ;close all;restoredefaultpath;
addpath(genpath('/home/ajoshi/git_sandbox/svreg-matlab/3rdParty'))
addpath(genpath('/home/ajoshi/git_sandbox/svreg-matlab/MEX_Files'))
addpath(genpath('/home/ajoshi/git_sandbox/svreg-matlab/src'))
addpath(genpath('/home/ajoshi/git_sandbox/svreg-matlab/dev'));

%addpath(genpath('/home/ajoshi/git_sandbox/svreg-matlab\dev\cortex_analysis'));

lst1=dir('/home/ajoshi/HCP900/data/*');

% build histograms at each voxel
% for jj=1:length(lst1)
%     fname=['E:\sipi_data\HCP900\data\',lst1(jj).name,'/colin_brain/',lst1(jj).name];
%     try
%         sub_SCT_hist(fname);
%     catch
%         fprintf(fname);
%     end
% end
%%
% Compute and save simulated ROI data

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% simulate changes in cortex proportional to IQ, ONLY for control subjects

load('/home/ajoshi/HCP900/xls_hcp900.mat');
% matlabpool 1
%
rng (0);
permutation = randperm (length(num(:,182)));
%
cogScores = num(:,182);
cogScores_perm = cogScores (permutation);
%
idNumbers = num (:,1);
idNumbers_perm = idNumbers (permutation);
%
save ('cogScores_perm.mat');
%
subno = 0
parfor jj=1:length(num(:,1))
    subname=['/home/ajoshi/HCP900/data/',num2str(num(jj,1)),'/colin_brain/',num2str(num(jj,1))];
    if isnan(cogScores_perm (jj)) % 
        continue
    end
    try
        %SCT_hist_atlas_right_all(:,:,subno)=single(tmp.SCT_hist_atlas_right(:,indx_reduced));
        cogScore = cogScores_perm (jj);
        modify=0.3*(100-cogScore)/50;
        modify=round((modify*100))/100;

        if ~exist([subname,'.SCT_hist_simulation_regress_',num2str(modify),'.mat'],'file')  
            subname
            subno = subno + 1
            sub_SCT_hist_simulation (subname, 16, modify);
            % sub_SCT_simulation (subname, 5, modify);
            close all;drawnow;
        end
    catch
        % subname
    end
end

%%
% subno=1;
% load('/home/ajoshi/HCP900/xls_hcp900.mat');
% 
% npts=1000;
% tar_orig=readdfs('/home/ajoshi/git_sandbox/svreg-matlab/BrainSuiteAtlas1/mri.right.mid.cortex.dfs');
% tar=reducepatch(tar_orig,10000);
% [~,ia,indx_reduced]=intersect(tar.vertices,tar_orig.vertices,'rows','stable');
% %[num,txt,raw]=xlsread('/home/ajoshi/HCP900/hcp_unrestricted_aajoshi_4_17_2016_14_56_18.csv');
% 
% for jj=1:length(num(:,1))
%     subname=['/home/ajoshi/HCP900/data/',num2str(num(jj,1)),'/colin_brain/',num2str(num(jj,1))];
%     if isnan(num(jj,182)) % 
%         continue;
%     end
%     try
%         sub_SCT_hist_simulation(subname,16,modify);
%     end
% end

%%
load('/home/ajoshi/HCP900/xls_hcp900.mat'); % reads 'num'

%npts = 1000;
tar_orig = readdfs('/home/ajoshi/git_sandbox/svreg-matlab/BrainSuiteAtlas1/mri.right.mid.cortex.dfs');
tar = reducepatch (tar_orig, 10000);
[~,ia,indx_reduced] = intersect (tar.vertices,tar_orig.vertices,'rows','stable');
%
SCT_atlas_right_all = single([]);
%
%cogScore = single ([]);
load ('cogScore_perm.mat',' cogScores_perm');
%
subno = 0
for jj=1:length(num(:,1))
    subname=['/home/ajoshi/HCP900/data/',num2str(num(jj,1)),'/colin_brain/',num2str(num(jj,1))];
    if isnan (cogScores_perm (jj)) % 
        continue;
    end
    try
        modify = 0.3 * (100 - cogScores_perm (jj)) / 50;
        modify = round ((modify*100))/100;
        
        if exist       ([subname,'.SCT_simulation_regress_',num2str(modify),'.mat'],'file')
            % get cortical features
            tmp = load ([subname,'.SCT_simulation_regress_' num2str(modify) '.mat']);
            subno = subno + 1
            SCT_atlas_right_all (:,:,subno) = single (tmp.SCTatlasright (indx_reduced,:)');
            % get IQ score
            cogScore (subno) = single (cogScores_perm (jj));
        end
    catch
        % subname
    end
end
SCT_atlas_right_all = permute (SCT_atlas_right_all, [1 3 2]);
% now dimensions of SCT_atlas_right_all :: (numOfFeatures x numOfSubects x numOfVoxels)
save sct_simulation_regress_measures_right_all

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% opengl software
%
addpath(genpath('/home/ajoshi/git_sandbox/svreg-matlab/src'))
addpath(genpath('/home/ajoshi/git_sandbox/svreg-matlab/dev'))
%
load ('sct_simulation_measures_right_all.mat', 'SCT_atlas_right_all', 'cogScore', 'tar');
%
sct_all = SCT_atlas_right_all;
%
numOfSubjects = size (sct_all, 2);
numOfVertices = size (sct_all, 3);
sctGram_all = nan (numOfSubjects, numOfSubjects, numOfVertices, 'single');
for jj = 1 : numOfVertices
    sctGram_all (:,:,jj) = exp (- 0.5 * ( pdist2 (sct_all (:,:,jj)', sct_all (:,:,jj)') ).^2); % isotropic Gaussian kernel with v = 1 (this is tuned automatically later)
end
whos

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% permutation testing - compute null distribution of extremal test statistic
% parpool close force local
% parpool 8
numOfPermutations = 50
[ testStatsUnpermuted_all, ...
    testStatsPermuted_Ftest_all, ...
    testStatsPermuted_Bartlett_all, ...
    testStatsPermuted_BrownForsythe_all, ...
    testStatsPermuted_LeveneAbsolute_all, ...
    testStatsPermuted_AnsariBradley_all ] ...
    = permutation_test_regression_all (sctGram_all, sct_all, cogScore, numOfPermutations);
save ('sct_simulation_testStats_right_all_nonlinear.mat', ...
    'testStatsUnpermuted_all', ...
    'testStatsPermuted_Ftest_all', ...
    'testStatsPermuted_Bartlett_all', ...
    'testStatsPermuted_BrownForsythe_all', ...
    'testStatsPermuted_LeveneAbsolute_all', ...
    'testStatsPermuted_AnsariBradley_all' );

%% visualize unpermuted p values
% permutation testing -compute p values
outlierRemovalPrctile = 90
%
[ pValUnpermuted_Ftest_all, testStatsPermutedExtreme_Ftest_all ] ...
    = permutation_test_p_values (testStatsPermuted_Ftest_all(:,1), testStatsPermuted_Ftest_all, outlierRemovalPrctile);

right_hemi_surf=readdfs('/home/ajoshi/coding_ground/svreg-matlab/BrainSuiteAtlas1/mri.right.mid.cortex.dfs');
right_hemi_surf=smooth_cortex_fast(tar,.5,300);
pp=(pValUnpermuted_Ftest_all'<0.02).*(0.02-pValUnpermuted_Ftest_all');
pv=smooth_surf_function(right_hemi_surf,pp,1,1);
figure;
patch('faces',right_hemi_surf.faces,'vertices',right_hemi_surf.vertices,'facevertexcdata',pv,'edgecolor','none','facecolor','interp');
axis equal;axis off; view (90,0); camlight; colorbar; material dull;
figure;
patch('faces',right_hemi_surf.faces,'vertices',right_hemi_surf.vertices,'facevertexcdata',testStatsPermuted_AnsariBradley_all(:,2),'edgecolor','none','facecolor','interp');
axis equal;axis off; view (-90,0); camlight; colorbar; material dull;

pause, close

% permutation testing -compute p values
outlierRemovalPrctile = 90
%
[ pValUnpermuted_Ftest_all, testStatsPermutedExtreme_Ftest_all ] ...
    = permutation_test_p_values (testStatsPermuted_Ftest_all(:,1), testStatsPermuted_Ftest_all, outlierRemovalPrctile);
hist (pValUnpermuted_Ftest_all); pause, close
hist (testStatsPermutedExtreme_Ftest_all, 100); pause, close
hist (testStatsPermuted_Ftest_all (:,1), 1000); pause, close
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

% visualize std dev in prediction error for all vertices
right_hemi_surf=readdfs('/home/ajoshi/git_sandbox/svreg-matlab/BrainSuiteAtlas1/mri.right.mid.cortex.dfs');
right_hemi_surf=smooth_cortex_fast(tar,.5,300);
figure;
%patch('faces',right_hemi_surf.faces,'vertices',right_hemi_surf.vertices,'facevertexcdata',mean(fit.predictionErrors_all)','edgecolor','none','facecolor','interp');
patch('faces',right_hemi_surf.faces,'vertices',right_hemi_surf.vertices,'facevertexcdata',std(fit.predictionErrors_all)','edgecolor','none','facecolor','interp');
axis equal;axis off; view (-90,0); camlight; colorbar; material dull;
pause, close

% visualize RMSE
right_hemi_surf=readdfs('/home/ajoshi/git_sandbox/svreg-matlab/BrainSuiteAtlas1/mri.right.mid.cortex.dfs');
right_hemi_surf=smooth_cortex_fast(tar,.5,300);
figure;
patch('faces',right_hemi_surf.faces,'vertices',right_hemi_surf.vertices,'facevertexcdata',fit.predictionRMSE_all,'edgecolor','none','facecolor','interp');
axis equal;axis off; view (-90,0); camlight; colorbar; material dull;
pause, close

% visualize conc
right_hemi_surf=readdfs('/home/ajoshi/git_sandbox/svreg-matlab/BrainSuiteAtlas1/mri.right.mid.cortex.dfs');
right_hemi_surf=smooth_cortex_fast(tar,.5,300);
figure;
patch('faces',right_hemi_surf.faces,'vertices',right_hemi_surf.vertices,'facevertexcdata',fit.conc_all,'edgecolor','none','facecolor','interp');
axis equal;axis off; view (-90,0); camlight; colorbar; material dull;
pause, close

% visualize p values
right_hemi_surf=smooth_cortex_fast(tar,.5,300);
figure;
%patch('faces',right_hemi_surf.faces,'vertices',right_hemi_surf.vertices,'facevertexcdata',testStats.p_Ftest_all','edgecolor','none','facecolor','interp');
%patch('faces',right_hemi_surf.faces,'vertices',right_hemi_surf.vertices,'facevertexcdata',testStats.p_Bartlett_all','edgecolor','none','facecolor','interp');
%patch('faces',right_hemi_surf.faces,'vertices',right_hemi_surf.vertices,'facevertexcdata',testStats.p_BrownForsythe_all','edgecolor','none','facecolor','interp');
%patch('faces',right_hemi_surf.faces,'vertices',right_hemi_surf.vertices,'facevertexcdata',testStats.p_LeveneAbsolute_all','edgecolor','none','facecolor','interp');
pfdr=0.05 % FDR(pValUnpermuted_AnsariBradley_all',0.05);
%patch('faces',right_hemi_surf.faces,'vertices',right_hemi_surf.vertices,'facevertexcdata',(pfdr-pValUnpermuted_Ftest_all').*double(pValUnpermuted_Ftest_all'< pfdr),'edgecolor','none','facecolor','interp');
patch('faces',right_hemi_surf.faces,'vertices',right_hemi_surf.vertices,'facevertexcdata',(pfdr-pValUnpermuted_LeveneAbsolute_all').*double(pValUnpermuted_LeveneAbsolute_all'< pfdr),'edgecolor','none','facecolor','interp');
%patch('faces',right_hemi_surf.faces,'vertices',right_hemi_surf.vertices,'facevertexcdata',(pfdr-pValUnpermuted_AnsariBradley_all').*double(pValUnpermuted_AnsariBradley_all'< pfdr),'edgecolor','none','facecolor','interp');
%patch('faces',right_hemi_surf.faces,'vertices',right_hemi_surf.vertices,'facevertexcdata',(pfdr-pValUnpermuted_Bartlett_all').*double(pValUnpermuted_Bartlett_all'< pfdr),'edgecolor','none','facecolor','interp');
%patch('faces',right_hemi_surf.faces,'vertices',right_hemi_surf.vertices,'facevertexcdata',(pfdr-pValUnpermuted_BrownForsythe_all').*double(pValUnpermuted_BrownForsythe_all'< pfdr),'edgecolor','none','facecolor','interp');
axis equal;axis off; camlight; colorbar; view (-90,0); camlight; material dull;%caxis([0,pfdr]);j1=flipud(jet);
pause, close
