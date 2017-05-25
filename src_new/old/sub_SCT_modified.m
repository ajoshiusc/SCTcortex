function sub_SCT(subname)
pth=fileparts(subname);

[SCTleft,SCTatlasleft]=cortical_SCT_measures_modified([subname,'.left.mid.cortex.svreg.dfs'],[pth,'/atlas.left.mid.cortex.svreg.dfs']);
[SCTright,SCTatlasright]=cortical_SCT_measures_modified([subname,'.right.mid.cortex.svreg.dfs'],[pth,'/atlas.right.mid.cortex.svreg.dfs']);

save([subname,'.SCT.mat'],'SCTleft','SCTright','SCTatlasleft','SCTatlasright');

