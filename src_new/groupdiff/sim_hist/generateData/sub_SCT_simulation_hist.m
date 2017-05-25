function sub_SCT_simulation_hist (mat_fname,subname,NRingsNbr,modify)
  
  if ~exist('NRingsNbr','var')
    NRingsNbr = 16;
  end
  NRingsNbr
  
  r = readdfs (['../../dfs/' subname,'.right.mid.cortex.svreg.dfs']);
  attrib = 1.0 - modify * double (r.labels==226);
  
  load ([ mat_fname '.SCT_hist.mat']);
  
  [ SCT_hist_right, SCT_hist_atlas_right ] ...
      = cortical_SCT_measures_simulation_hist (['../../dfs/' subname,'.right.mid.cortex.svreg.dfs'], ['../../dfs/atlas.right.mid.cortex.svreg.dfs'], NRingsNbr, attrib);
  
  save ([subname,'.SCT_hist_simulation_',num2str(modify),'.mat'], 'SCT_hist_right', 'SCT_hist_atlas_right');
  
  return
  