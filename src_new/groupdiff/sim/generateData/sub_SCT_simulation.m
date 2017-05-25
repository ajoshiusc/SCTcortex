function sub_SCT_simulation (mat_fname,subname,NLevels,modify)
  
  if ~exist('NLevels','var')
    NLevels = 5;
  end
  NLevels
  
  r = readdfs(['../../dfs/' subname,'.right.mid.cortex.svreg.dfs']);
  attrib = 1.0 - modify * double (r.labels == 226);
  
  load ([ mat_fname '.SCT.mat']);
  
  [ SCT_right, SCT_atlas_right ] ...
      = cortical_SCT_measures_simulation (['../../dfs/' subname,'.right.mid.cortex.svreg.dfs'], ['../../dfs/atlas.right.mid.cortex.svreg.dfs'], NLevels, attrib);
  
  save ([subname,'.SCT_simulation_',num2str(modify),'.mat'], 'SCT_right', 'SCT_atlas_right');
  
  return
  