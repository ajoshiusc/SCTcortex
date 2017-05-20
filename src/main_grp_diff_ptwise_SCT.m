clc;clear all;close all;

% This script performs pointwise group differences in cortical thicknesses
% for the right hemisphere
% modify as appropriate for left hemisphere

%group1 and group2 directory names
addpath(genpath('C:\Users\ajoshi\Documents\coding_ground\svreg-matlab\dev'));
addpath(genpath('C:\Users\ajoshi\Documents\coding_ground\svreg-matlab\src'));

dirname1='C:\Users\ajoshi\Downloads\oasis_data';
dirname2='C:\Users\ajoshi\Downloads\Beijing';
NLevels=1;
aa1=dir(dirname1);
aa1=aa1(3:end);
aa2=dir(dirname2);
aa2=aa2(3:end);

jjj=1;
npts=1000;
tar=readdfs('C:\Users\ajoshi\Documents\coding_ground\svreg-matlab/BrainSuiteAtlas1\mri.right.mid.cortex.dfs');

parfor jj=1:length(aa1)
    subbasename=sprintf('%s/%s/%s_mpr_n4_anon_sbj_111.RAS',dirname1,aa1(jj).name,aa1(jj).name);
    sub_SCT_hist(subbasename,NLevels);
    fprintf('%d/%d SCT done\n',jj,length(aa1));
end
%
parfor jj=1:length(aa2)
    subbasename=sprintf('%s/%s/anat/mprage_anonymized',dirname2,aa2(jj).name);
    sub_SCT_hist(subbasename,NLevels);
    fprintf('%d/%d SCT done\n',jj,length(aa1));
end

for ptset=1:npts:length(tar.vertices)
    ptst1=ptset:min(ptset+npts-1,length(tar.vertices));clear Thickness1 Thickness2
    jj1=1;
    clear SCTgrp1 SCTgrp2;
    for jj=1:length(aa1)
        subname=sprintf('%s/%s/%s_mpr_n4_anon_sbj_111.RAS',dirname1,aa1(jj).name,aa1(jj).name);
        load([subname,'.SCT.mat']);
        SCTgrp1(:,:,jj1)=SCTatlasright(ptst1,:);        jj1=jj1+1;
    end
    jj1=1;
    
    for jj=1:length(aa2)
        subname=sprintf('%s/%s/anat/mprage_anonymized',dirname2,aa2(jj).name);
        load([subname,'.SCT.mat']);
        SCTgrp2(:,:,jj1)=SCTatlasright(ptst1,:);        jj1=jj1+1;
    end
    %   Thickness1(:,13)=[];
    [p(ptst1),~]=hotelling_t2_test(permute(SCTgrp1,[2,3,1]),permute(SCTgrp2,[2,3,1]));
    
    disp(sprintf('%d/%d',ptset, length(tar.vertices)));
    
end

pfdr=FDR(p,0.05);
padj=min(1,p*0.05/pfdr);
h=figure;
patch('vertices',tar.vertices,'faces',tar.faces,'facevertexcdata',1.0*(0.0+padj').*(padj'<10.05),'edgecolor','none','facecolor','interp');axis equal;axis off;camlight;material dull;lighting phong;
colormap jet;colorbar;
tar.vcolor=0*tar.vcolor;tar.vcolor(:,3)=1;
%tar.vcolor(p<0.05,3)=0;tar.vcolor(p<0.05,1)=1;
writedfs('for_bci_linked_dist_right.dfs',tar);
saveas(h,'for_bci_linked_dist_right.fig');
