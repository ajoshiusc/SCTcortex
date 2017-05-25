opengl software

clc;clear ;close all;restoredefaultpath;
addpath(genpath('/home/ajoshi/coding_ground/svreg-matlab/3rdParty'))
addpath(genpath('/home/ajoshi/coding_ground/svreg-matlab/MEX_Files'))
addpath(genpath('/home/ajoshi/coding_ground/svreg-matlab/src'))
addpath(genpath('/home/ajoshi/coding_ground/svreg-matlab/dev'));

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
subno = 0
% matlabpool 1
for jj=1:length(num(:,1))
    subname=['/home/ajoshi/HCP900/data/',num2str(num(jj,1)),'/colin_brain/',num2str(num(jj,1))];
    if isnan(num(jj,182)) % 
        continue;
    end
    try
        %SCT_hist_atlas_right_all(:,:,subno)=single(tmp.SCT_hist_atlas_right(:,indx_reduced));
        modify=0.3*(100-num(jj,182))/50;
        modify=round((modify*100))/100;

        if exist([subname,'.SCT_hist_simulation_',num2str(modify),'.mat'],'file')        
            % sub_SCT_hist_simulation(subname, 16, modify);
            sub_SCT_simulation (subname, 5, modify);
        end
        subno=subno+1
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

subno = 0
for jj=1:length(num(:,1))
    subname=['/home/ajoshi/HCP900/data/',num2str(num(jj,1)),'/colin_brain/',num2str(num(jj,1))];
    if isnan(num(jj,182)) % 
        continue;
    end
    try
        modify=0.3*(100-num(jj,182))/50;
        modify=round((modify*100))/100;
        
        if exist       ([subname,'.SCT_hist_simulation_',num2str(modify),'.mat'],'file')
            % get cortical features
            tmp = load ([subname,'.SCT_hist_simulation_' num2str(modify) '.mat']);
            subno = subno + 1
            SCT_hist_atlas_right_all (:,:,subno) = single (tmp.SCT_hist_atlas_right (:,indx_reduced));
            % get IQ score
            cogScore (subno) = num(jj,182);
        end
    catch
        % subname
    end
end
SCT_hist_atlas_right_all = permute (SCT_hist_atlas_right_all, [1 3 2]);
%SCT_hist_atlas_right_all = permute (SCT_hist_atlas_right_all, [1 3 2]);
% now dimensions of SCT_hist_atlas_right_all, SCT_hist_atlas_right_all are
% (numOfsubjects x histLength x numOfVoxels)
save sct_hist_simulation_measures_right_all


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%opengl software
%
addpath(genpath('/home/ajoshi/git_sandbox/svreg-matlab/src'))
addpath(genpath('/home/ajoshi/git_sandbox/svreg-matlab/dev'))
%
load ('sct_hist_simulation_measures_right_all.mat', 'SCT_hist_atlas_right_all', 'cogScore', 'tar');
% get hists
hists = SCT_hist_atlas_right_all;
% map hists to hypersphere
histsSqrt = sqrt (hists);
% get dot products of all pairs of sqrt(hists)
numOfSubjects = size (SCT_hist_atlas_right_all, 2);
numOfVertices = size (SCT_hist_atlas_right_all, 3);
%
histsSqrtCrossExp_all = nan (numOfSubjects, numOfSubjects, numOfVertices, 'single');
for jj = 1 : numOfVertices
    histsSqrtCrossExp_all (:,:,jj) = exp (histsSqrt (:,:,jj)' * histsSqrt (:,:,jj));
end
clear histsSqrt SCT_hist_simulation_atlas_right_all SCT_hist_atlas_right_all
%save histsSqrtCrossExp_all % do NOT save this. Take way too much time.
whos

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%load histsSqrtCrossExp_all

% permutation testing - compute null distribution of extremal test statistic
% matlabpool close force local
% matlabpool 4
numOfPermutations = 50
[ testStatsUnpermuted_all, ...
    testStatsPermuted_Ftest_all, ...
    testStatsPermuted_Bartlett_all, ...
    testStatsPermuted_BrownForsythe_all, ...
    testStatsPermuted_LeveneAbsolute_all, ...
    testStatsPermuted_AnsariBradley_all ] ...
    = permutation_test_hist_regression_all (histsSqrtCrossExp_all, hists, cogScore, numOfPermutations);
save ('sct_hist_simulation_testStats_right_nonlinear.mat', ...
    'testStatsUnpermuted_all', ...
    'testStatsPermuted_Ftest_all', ...
    'testStatsPermuted_Bartlett_all', ...
    'testStatsPermuted_BrownForsythe_all', ...
    'testStatsPermuted_LeveneAbsolute_all', ...
    'testStatsPermuted_AnsariBradley_all' );

%% visualize unpermuted p values
right_hemi_surf=readdfs('/home/ajoshi/coding_ground/svreg-matlab/BrainSuiteAtlas1/mri.right.mid.cortex.dfs');
right_hemi_surf=smooth_cortex_fast(tar,.5,300);
pp=(pValUnpermuted_Ftest_all'<0.02).*(0.02-pValUnpermuted_Ftest_all');
pv=smooth_surf_function(right_hemi_surf,pp,1,1);
figure;
patch('faces',right_hemi_surf.faces,'vertices',right_hemi_surf.vertices,'facevertexcdata',pv,'edgecolor','none','facecolor','interp');
axis equal;axis off; view (90,0); camlight; colorbar; material dull; caxis([0,0.02])
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
    = permutation_test_p_values (testStatsPermuted_Bartlett_all(:,1), testStatsPermuted_Bartlett_all, outlierRemovalPrctile);
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