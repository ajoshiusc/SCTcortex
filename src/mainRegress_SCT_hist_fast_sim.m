
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
for jj=1:length(num(:,1))
    subname=['/home/ajoshi/HCP900/data/',num2str(num(jj,1)),'/colin_brain/',num2str(num(jj,1))];
    if isnan(cogScores_perm (jj)) % 
        continue
    end
    try
        %SCT_hist_atlas_right_all(:,:,subno)=single(tmp.SCT_hist_atlas_right(:,indx_reduced));
        cogScore = cogScores_perm (jj);
        modify=0.3*(100-cogScore)/50;
        modify=round((modify*100))/100;

        if exist([subname,'.SCT_hist_simulation_regress_',num2str(modify),'.mat'],'file')  
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
clear all
tar_orig = readdfs('/home/ajoshi/git_sandbox/svreg-matlab/BrainSuiteAtlas1/mri.right.mid.cortex.dfs');
tar = reducepatch (tar_orig, 10000);
[~,ia,indx_reduced] = intersect (tar.vertices,tar_orig.vertices,'rows','stable');
%
SCT_hist_atlas_right_all = single([]);
%
load ('cogScores_perm.mat', 'cogScores_perm', 'idNumbers');
%
subno = 0
for jj = 1 : length (cogScores_perm)
    subname=['/home/ajoshi/HCP900/data/',num2str(idNumbers(jj)),'/colin_brain/',num2str(idNumbers(jj))];
%    subname=['/home/ajoshi/HCP900/data/',num2str(num(jj,1)),'/colin_brain/',num2str(num(jj,1))];

    if isnan (cogScores_perm (jj)) %
        continue;
    end
    %    try
    modify = 0.3 * (100 - cogScores_perm (jj)) / 50;
    modify = round ((modify*100))/100;
%        cogScore = cogScores_perm (jj);
%        modify=0.3*(100-cogScore)/50;
%        modify=round((modify*100))/100;
    
%    exist([subname,'.SCT_hist_simulation_regress_',num2str(modify),'.mat'],'file')  
    filename = [subname,'.SCT_hist_simulation_regress_',num2str(modify),'.mat'];
    if exist       (filename,'file')
        % get cortical features
%        tmp = load (filename);
        subno = subno + 1
%        SCT_hist_atlas_right_all (:,:,subno) = single (tmp.SCT_hist_atlas_right (:,indx_reduced));
        % get IQ score
        cogScore (subno) = single (cogScores_perm (jj));
        %
        whos SCT_hist_atlas_right_all cogScore %, pause
    end
    %     catch
    %         % subname
    %     end
end
SCT_hist_atlas_right_all = permute (SCT_hist_atlas_right_all, [1 3 2]);
% now dimensions of SCT_atlas_right_all :: (numOfFeatures x numOfSubects x numOfVoxels)
save sct_hist_simulation_regress_measures_right_all

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% opengl software
%
clear all
addpath(genpath('/home/ajoshi/git_sandbox/svreg-matlab/src'))
addpath(genpath('/home/ajoshi/git_sandbox/svreg-matlab/dev'))
%
load ('sct_hist_simulation_regress_measures_right_all.mat', 'SCT_hist_atlas_right_all', 'cogScore', 'tar');
%
% get dot products of all pairs of sqrt(hists)
numOfSubjects = size (SCT_hist_atlas_right_all, 2);
numOfVertices = size (SCT_hist_atlas_right_all, 3);
histsSqrtCrossExp_all = nan (numOfSubjects, numOfSubjects, numOfVertices, 'single');
for jj = 1 : numOfVertices
    histsSqrtCrossExp_all (:,:,jj) = exp (SCT_hist_atlas_right_all (:,:,jj)' * SCT_hist_atlas_right_all (:,:,jj));
end
%SCT_hist_atlas_right_all = [];
%histsSqrtCrossExp_all = [];
whos

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% permutation testing - compute null distribution of extremal test statistic
% parpool close force local
% parpool 8
numOfPermutations = 100
for nb = 1:20
    '--------------------------------------------------------', nb
    % create a bootstrap sample
    if nb == 1
        bootstrapSampleIndices = [1:numOfSubjects];
    else
        rng (100 + nb);
        [ bootstrapSampleDummy bootstrapSampleIndices ]  = datasample ([1:numOfSubjects], numOfSubjects);
        bootstrapSampleIndices = sort (bootstrapSampleIndices); % sort to reduce data movement in next step
    end
    % resample the data
    [ testStatsUnpermuted_all, ...
        testStatsPermuted_Ftest_all, ...
        testStatsPermuted_Bartlett_all, ...
        testStatsPermuted_BrownForsythe_all, ...
        testStatsPermuted_LeveneAbsolute_all, ...
        testStatsPermuted_AnsariBradley_all ] ...
        = permutation_test_hist_regression_all ...
        (histsSqrtCrossExp_all, ... % (bootstrapSampleIndices, bootstrapSampleIndices, :), ...
        SCT_hist_atlas_right_all (:, bootstrapSampleIndices, :), ...
        cogScore (bootstrapSampleIndices), ...
        numOfPermutations);
    %
    testStatsUnpermuted_all_boot (:,nb)              = testStatsUnpermuted_all;
    testStatsPermuted_Ftest_all_boot (:,:,nb)          = testStatsPermuted_Ftest_all;
    testStatsPermuted_Bartlett_all_boot (:,:,nb)       = testStatsPermuted_Bartlett_all;
    testStatsPermuted_BrownForsythe_all_boot (:,:,nb)  = testStatsPermuted_BrownForsythe_all;
    testStatsPermuted_LeveneAbsolute_all_boot (:,:,nb) = testStatsPermuted_LeveneAbsolute_all;
    testStatsPermuted_AnsariBradley_all_boot (:,:,nb)  = testStatsPermuted_AnsariBradley_all;
    %
    save ('sct_hist_simulation_testStats_right_all_linear_new.mat', ...
        'testStatsUnpermuted_all_boot', ...
        'testStatsPermuted_Ftest_all_boot', ...
        'testStatsPermuted_Bartlett_all_boot', ...
        'testStatsPermuted_BrownForsythe_all_boot', ...
        'testStatsPermuted_LeveneAbsolute_all_boot', ...
        'testStatsPermuted_AnsariBradley_all_boot' );
end

%% permutation testing -compute p values
outlierRemovalPrctile = 97
for nb = 1 : size (testStatsUnpermuted_all_boot, 2)
    %
    [ pValUnpermuted_Ftest_all, testStatsPermutedExtreme_Ftest_all ] ...
        = permutation_test_p_values (testStatsPermuted_Ftest_all_boot (:,1,nb), testStatsPermuted_Ftest_all_boot (:,:,nb), outlierRemovalPrctile);
    pValUnpermuted_Ftest_all_boot (:,nb) = pValUnpermuted_Ftest_all;
    %
    %[ pValUnpermuted_Bartlett_all, testStatsPermutedExtreme_Bartlett_all ] ...
    %    = permutation_test_p_values (testStatsPermuted_Bartlett_all(:,1), testStatsPermuted_Bartlett_all (:,:,1), outlierRemovalPrctile);
    %
    %[ pValUnpermuted_BrownForsythe_all, testStatsPermutedExtreme_BrownForsythe_all ] ...
    %    = permutation_test_p_values (testStatsPermuted_BrownForsythe_all(:,1), testStatsPermuted_BrownForsythe_all (:,:,1), outlierRemovalPrctile);
    %
    %[ pValUnpermuted_LeveneAbsolute_all, testStatsPermutedExtreme_LeveneAbsolute_all ] ...
    %    = permutation_test_p_values (testStatsPermuted_LeveneAbsolute_all(:,1), testStatsPermuted_LeveneAbsolute_all (:,:,1), outlierRemovalPrctile);
    %
    %[ pValUnpermuted_AnsariBradley_all, testStatsPermutedExtreme_AnsariBradley_all ] ...
    %    = permutation_test_p_values (testStatsPermuted_AnsariBradley_all(:,1), testStatsPermuted_AnsariBradley_all (:,:,1), outlierRemovalPrctile);
end

% visualize p values
tar_orig = readdfs('/home/ajoshi/git_sandbox/svreg-matlab/BrainSuiteAtlas1/mri.right.mid.cortex.dfs');
tar = reducepatch (tar_orig, 10000);
right_hemi_surf=readdfs('/home/ajoshi/git_sandbox/svreg-matlab/BrainSuiteAtlas1/mri.right.mid.cortex.dfs');
right_hemi_surf=smooth_cortex_fast(tar,.1,200);
%
thresh = 0.1
%
pp (:,1) = (pValUnpermuted_Ftest_all_boot(:,1)'<thresh).*(thresh-pValUnpermuted_Ftest_all_boot(:,1)');
%pp (:,2) = (pValUnpermuted_Bartlett_all'      <thresh).*(thresh-pValUnpermuted_Bartlett_all');
%pp (:,3) = (pValUnpermuted_BrownForsythe_all' <thresh).*(thresh-pValUnpermuted_BrownForsythe_all');
%pp (:,4) = (pValUnpermuted_LeveneAbsolute_all'<thresh).*(thresh-pValUnpermuted_LeveneAbsolute_all');
%pp (:,5) = (pValUnpermuted_AnsariBradley_all' <thresh).*(thresh-pValUnpermuted_AnsariBradley_all');
cmap=cbrewer('seq','Oranges',100);
for testNum = 1 : 1 % 5
    testNum
    %
    pv (:,testNum) = smooth_surf_function(right_hemi_surf, pp (:,testNum),1,1);
    %
    figure;
    patch('faces',right_hemi_surf.faces,'vertices',right_hemi_surf.vertices,'facevertexcdata',pp (:,testNum),'edgecolor','none','facecolor','interp');
    axis equal;axis off; view (90,0); camlight; colorbar; material dull;colormap(cmap);lighting phong;
    % pause, close
    % export_fig -a1 -m2 -transparent Pics/Var_Hist_RiemannianKernel_SCT.png, close
end
figure;
patch('faces',right_hemi_surf.faces,'vertices',right_hemi_surf.vertices,'facevertexcdata',var(pValUnpermuted_Ftest_all_boot,[],2),'edgecolor','none','facecolor','interp');
axis equal;axis off; view (90,0); camlight; colorbar; material dull;colormap(cmap);lighting phong;caxis([0 0.3]);


%% visualize std dev in prediction error for all vertices
right_hemi_surf=readdfs('/home/ajoshi/git_sandbox/svreg-matlab/BrainSuiteAtlas1/mri.right.mid.cortex.dfs');
right_hemi_surf=smooth_cortex_fast(tar,.5,300);
figure;
%patch('faces',right_hemi_surf.faces,'vertices',right_hemi_surf.vertices,'facevertexcdata',mean(fit.predictionErrors_all)','edgecolor','none','facecolor','interp');
patch('faces',right_hemi_surf.faces,'vertices',right_hemi_surf.vertices,'facevertexcdata',std(fit.predictionErrors_all)','edgecolor','none','facecolor','interp');
axis equal;axis off; view (-90,0); camlight; colorbar; material dull;
pause, close

%% visualize RMSE
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
