clc;clear all;close all;

% This script performs pointwise group differences in cortical thicknesses
% for the right hemisphere
% modify as appropriate for left hemisphere

%group1 and group2 directory names
addpath(genpath('C:\Users\ajoshi\Documents\coding_ground\svreg-matlab\dev'));
addpath(genpath('C:\Users\ajoshi\Documents\coding_ground\svreg-matlab/3rdParty'))
addpath(genpath('C:\Users\ajoshi\Documents\coding_ground\svreg-matlab/MEX_Files'))
addpath(genpath('C:\Users\ajoshi\Documents\coding_ground\svreg-matlab/src'))

% dirname1='C:\Users\ajoshi\Downloads\oasis_data';

% % generate SCT histograms for group 1, save them
% parfor jj=1:length(aa1)
%     subbasename=sprintf('%s/%s/%s_mpr_n4_anon_sbj_111.RAS',dirname1,aa1(jj).name,aa1(jj).name);
%     sub_SCT_hist(subbasename,NLevels);
%     fprintf('%d/%d SCT done\n',jj,length(aa1));
% end
% % % generate SCT histograms for group 2, save them
% for jj=1:length(aa2)
%      subbasename=sprintf('%s/%s/anat/mprage_anonymized',dirname2,aa2(jj).name);
%      sub_SCT_hist_simulation(subbasename,NLevels);
%      fprintf('%d/%d SCT done\n',jj,length(aa1));
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
numOfHistBins = size (tmp.SCT_hist_atlas_left, 1);

numOfSubjectsGrp1 = sum (CDRs == 0);
fnamesGrp1 = fnames (CDRs == 0);
count = 0;
for modify = 0.1 % [0.01,0.05,0.1,0.25]
    parfor jj = 1 : numOfSubjectsGrp1
        mat_fname=['E:\sipi_data\oasis_full_data\data\',fnamesGrp1{jj},'\PROCESSED\MPRAGE\T88_111\',fnamesGrp1{jj},'_aaj'];
        if 1%exist([mat_fname,'.SCT_hist.mat'],'file')
            if jj<=50
                sub_SCT_hist_simulation(mat_fname,16,modify);
%            else
%                sub_SCT_hist_simulation(mat_fname,16,0);
                %copyfile([mat_fname,'.SCT_hist.mat'],[mat_fname,'.SCT_hist_simulation.mat'])
            end
        end
        fprintf('%s : %d/%d done\n',mat_fname,jj,numOfSubjectsGrp1);
    end
end
