function sub_SCT(subname,NLevels)
if ~exist('NLevels','var')
    NLevels=5;
end

pth=fileparts(subname);
[SCTleft,SCTatlasleft]= cortical_SCT_measures([subname,'.left.mid.cortex.svreg.dfs'],[pth,'/atlas.left.mid.cortex.svreg.dfs'],NLevels)
%[SCTleft,SCTatlasleft]=cortical_SCT_measures([subname,'.left.mid.cortex.svreg.dfs'],[pth,'/atlas.left.mid.cortex.svreg.dfs'],NLevels);
[SCTright,SCTatlasright]=cortical_SCT_measures([subname,'.right.mid.cortex.svreg.dfs'],[pth,'/atlas.right.mid.cortex.svreg.dfs'],NLevels);

save([subname,'.SCT.mat'],'SCTleft','SCTright','SCTatlasleft','SCTatlasright');

