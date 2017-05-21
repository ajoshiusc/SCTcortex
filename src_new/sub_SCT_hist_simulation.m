function sub_SCT_hist_simulation (subname,NRingsNbr,modify)
  
  if ~exist('NRingsNbr','var')
    NRingsNbr = 16;
  end
  NRingsNbr
  
  r = readdfs ([subname,'.right.mid.cortex.svreg.dfs']);
  attrib = 1.0 - modify * double (r.labels==226);
  
  pth = fileparts(subname);
  
  load ([subname,'.SCT_hist.mat']);
  
  [ SCT_hist_right, SCT_hist_atlas_right ] ...
      = cortical_SCT_measures_hist_simulation ([subname,'.right.mid.cortex.svreg.dfs'], [pth,'/atlas.right.mid.cortex.svreg.dfs'], NRingsNbr, attrib);
  
  save ([subname,'.SCT_hist_simulation_regress_',num2str(modify),'.mat'], 'SCT_hist_right', 'SCT_hist_atlas_right');
  
  return
  