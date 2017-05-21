function sub_SCT_hist (subname,NRingsNbr)

if ~exist('NRingsNbr','var')
    NRingsNbr=16;
end

pth=fileparts(subname);
[ SCT_hist_left , SCT_hist_atlas_left  ] = cortical_SCT_measures_hist ([subname,'.left.mid.cortex.svreg.dfs'], [pth,'/atlas.left.mid.cortex.svreg.dfs'], NRingsNbr);
[ SCT_hist_right, SCT_hist_atlas_right ] = cortical_SCT_measures_hist ([subname,'.right.mid.cortex.svreg.dfs'],[pth,'/atlas.right.mid.cortex.svreg.dfs'],NRingsNbr);

save ([subname,'.SCT_hist.mat'],'SCT_hist_left','SCT_hist_right','SCT_hist_atlas_left','SCT_hist_atlas_right');

return