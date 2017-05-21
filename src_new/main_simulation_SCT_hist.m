
clc
clear all
close all

%group1 and group2 directory names
addpath(genpath('C:\Users\ajoshi\Documents\coding_ground\svreg-matlab\dev'));
addpath(genpath('C:\Users\ajoshi\Documents\coding_ground\svreg-matlab/3rdParty'))
addpath(genpath('C:\Users\ajoshi\Documents\coding_ground\svreg-matlab/MEX_Files'))
addpath(genpath('C:\Users\ajoshi\Documents\coding_ground\svreg-matlab/src'))

% xlsread works only on windows
%[N,T,R]=xlsread('E:\sipi_data\oasis_full_data\data\oasis_cross-sectional.csv');
%fnames= R(2:end,1);
%temp  = R(2:end,8);
%for jj=1:length(fnames)
%    CDRs(jj) = temp{jj};
%end

%
fid = fopen ('../sct_data/oasis_cross-sectional.csv');
out = textscan (fid, '%s%s%s%s%s%s%s%s%s%s%s%s', 'delimiter', ',');
fclose (fid);
%
for jj = 1 : length (out{1}) - 1
  % filename
  fnames{jj} = out{1}{1+jj};
  % CDR dementia rating
  if isempty (out{8}{1+jj})
    CDRs(jj) = nan;
  else
    CDRs(jj) = str2num (out{8}{1+jj});
  end
end

% get size of each histogram to preallocate
mat_fname = ['../sct_data/sct_data/',fnames{1},'_aaj.SCT.mat'];
tmp = load (mat_fname);
numOfHistBins = size (tmp.SCT_hist_atlas_left, 1);

numOfSubjectsGrp1 = sum (CDRs == 0);
fnamesGrp1 = fnames (CDRs == 0);
count = 0;
for modify = 0.1 % [0.01,0.05,0.1,0.25]
  for jj = 1 : numOfSubjectsGrp1
    mat_fname = ['../sct_data/sct_data/',fnamesGrp1{jj},'_aaj'];
    if 1 % exist([mat_fname,'.SCT_hist.mat'],'file')
      if jj <= 50
        sub_SCT_hist_simulation (mat_fname, 16, modify);
      end
    end
    fprintf('%s : %d/%d done\n',mat_fname,jj,numOfSubjectsGrp1);
  end
end
