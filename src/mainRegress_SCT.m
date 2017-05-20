
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

% read all hists
subno=1;
load('/home/ajoshi/HCP900/xls_hcp900.mat');
%[num,txt,raw]=xlsread('/home/ajoshi/HCP900/hcp_unrestricted_aajoshi_4_17_2016_14_56_18.csv');

npts=1000;
tar_orig=readdfs('/home/ajoshi/git_sandbox/svreg-matlab/BrainSuiteAtlas1/mri.left.mid.cortex.dfs');
tar=reducepatch(tar_orig,100000);
[~,ia,indx_reduced]=intersect(tar.vertices,tar_orig.vertices,'rows','stable');

for jj=1:length(num(:,1))
    subname=['/home/ajoshi/HCP900/data/',num2str(num(jj,1)),'/colin_brain/',num2str(num(jj,1))];
    if isnan(num(jj,182))
        continue;
    end
    try
        tmp=load([subname,'.SCT_hist.mat']);
        SCT_hist_atlas_left_all(:,:,subno)=single(tmp.SCT_hist_atlas_left(:,indx_reduced));
        %SCT_hist_atlas_right_all(:,:,subno)=single(tmp.SCT_hist_atlas_right(:,indx_reduced));
        cogScore (subno)=num(jj,182);
        subno=subno+1
    catch
        % subname
    end
end
SCT_hist_atlas_left_all = permute (SCT_hist_atlas_left_all, [1 3 2]);
%SCT_hist_atlas_right_all = permute (SCT_hist_atlas_right_all, [1 3 2]);
% now dimensions of SCT_hist_atlas_left_all, SCT_hist_atlas_right_all are
% (numOfsubjects x histLength x numOfVoxels)
load sct_hist_measures_left_all_hires

% (1)
% at each voxel, perform regression (Nadaraya-Watson)
% optimize bandwidth parameter (LOO cross validation) to optimize expected prediction error
% compute prediction error, for each subject
for jj=1:size(SCT_hist_atlas_left_all,3)
    if mod(jj,10) == 0
        jj
    end
    hists = SCT_hist_atlas_left_all(:,:,jj);
    values = cogScore;
    [ valuesPredicted predictionErrors predictionRMSE conc ] = kernel_regression_hyperSphere (sqrt(hists), values);
    %
    valuesPredicted_all (:,jj) = valuesPredicted;
    conc_all (jj) = conc;
    predictionErrors_all (:,jj) = predictionErrors;
    predictionRMSE_all (jj) = predictionRMSE;
    %
    % plot (values, valuesPredicted, 'ro'); grid on; axis equal tight; pause, close
end
'done all regressions'
save sct_hist_regression_left_all_hires

% visualize of std dev in prediction error for all vertices
left_hemi_surf=readdfs('/home/ajoshi/git_sandbox/svreg-matlab/BrainSuiteAtlas1/mri.left.mid.cortex.dfs');
left_hemi_surf=smooth_cortex_fast(tar,.5,300);
figure;
%patch('faces',left_hemi_surf.faces,'vertices',left_hemi_surf.vertices,'facevertexcdata',mean(predictionErrors_all)','edgecolor','none','facecolor','interp');
patch('faces',left_hemi_surf.faces,'vertices',left_hemi_surf.vertices,'facevertexcdata',std(predictionErrors_all)','edgecolor','none','facecolor','interp');
axis equal;axis off; view (-90,0); camlight; colorbar; material dull;
%pause, close

% visualize RMSE
left_hemi_surf=readdfs('/home/ajoshi/git_sandbox/svreg-matlab/BrainSuiteAtlas1/mri.left.mid.cortex.dfs');
left_hemi_surf=smooth_cortex_fast(tar,.5,300);
figure;
patch('faces',left_hemi_surf.faces,'vertices',left_hemi_surf.vertices,'facevertexcdata',predictionRMSE_all','edgecolor','none','facecolor','interp');
axis equal;axis off; view (-90,0); camlight; colorbar; material dull;
%pause, close

% visualize conc
left_hemi_surf=readdfs('/home/ajoshi/git_sandbox/svreg-matlab/BrainSuiteAtlas1/mri.left.mid.cortex.dfs');
left_hemi_surf=smooth_cortex_fast(tar,.5,300);
figure;
patch('faces',left_hemi_surf.faces,'vertices',left_hemi_surf.vertices,'facevertexcdata',conc_all','edgecolor','none','facecolor','interp');
axis equal;axis off; view (-90,0); camlight; colorbar; material dull;
%pause, close

% (2)
% at each voxel, perform regression (constant function)
% compute prediction error, for each subject
valuesPredictedConst = mean (cogScore);
predictionErrors_const = cogScore - valuesPredictedConst;
predictionRMSE_const = sqrt(mean(predictionErrors_const.^2))

% (3)
% at each voxel, perform 2-sample F test between the previously computed 2 groups of perdiction error values
sample1 = predictionErrors_const;
for jj=1:size(SCT_hist_atlas_left_all,3)
    if mod(jj,1e3) == 0
        jj
    end
    %
    sample2 = predictionErrors_all (:,jj);
    % F test for equality of variances
    [h,p,ci,stats] = vartest2 (sample1, sample2);
    stats_Ftest_all (jj) = stats.fstat;
    p_Ftest_all (jj) = p;
    % Bartlett's test for equality of variances
    [p,stats] = vartestn ([sample1', sample2], 'Display','off');
    stats_Bartlett_all (jj) = stats.chisqstat;
    p_Bartlett_all (jj) = p;
    % Brown-Forsythe's test for equality of variances
    [p,stats] = vartestn ([sample1', sample2], 'TestType', 'BrownForsythe', 'Display','off');
    stats_BrownForsythe_all (jj) = stats.fstat;
    p_BrownForsythe_all (jj) = p;
    % Levene's test (absolute) for equality of variances
    [p,stats] = vartestn ([sample1', sample2], 'TestType', 'LeveneAbsolute', 'Display','off');
    stats_LeveneAbsolute_all (jj) = stats.fstat;
    p_LeveneAbsolute_all (jj) = p;
    % Ansari-Bradley test for equality of variances
    [h,p,stats] = ansaribradley (sample1, sample2);
    stats_AnsariBradley_all (jj) = stats.W;
    p_AnsariBradley_all (jj) = p;
end
save sct_hist_regression_hyptest_left_hires

% visualize p values
left_hemi_surf=smooth_cortex_fast(tar,.5,300);
figure;
%patch('faces',left_hemi_surf.faces,'vertices',left_hemi_surf.vertices,'facevertexcdata',p_Ftest_all','edgecolor','none','facecolor','interp');
%patch('faces',left_hemi_surf.faces,'vertices',left_hemi_surf.vertices,'facevertexcdata',p_Bartlett_all','edgecolor','none','facecolor','interp');
%patch('faces',left_hemi_surf.faces,'vertices',left_hemi_surf.vertices,'facevertexcdata',p_BrownForsythe_all','edgecolor','none','facecolor','interp');
%patch('faces',left_hemi_surf.faces,'vertices',left_hemi_surf.vertices,'facevertexcdata',p_LeveneAbsolute_all','edgecolor','none','facecolor','interp');
pfdr=FDR(p_AnsariBradley_all,0.05);
patch('faces',left_hemi_surf.faces,'vertices',left_hemi_surf.vertices,'facevertexcdata',(pfdr-p_AnsariBradley_all').*double(p_AnsariBradley_all'<pfdr),'edgecolor','none','facecolor','interp');
axis equal;axis off;camlight;camlight; colorbar; material dull
%pause, close
