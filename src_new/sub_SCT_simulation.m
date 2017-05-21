function sub_SCT_simulation (subname,NLevels,modify)
  
  if ~exist('NLevels','var')
    NLevels = 5;
  end
  NLevels
  
  r = readdfs([subname,'.right.mid.cortex.svreg.dfs']);
  attrib = 1.0 - modify * double (r.labels == 226);
  
  pth = fileparts(subname);
  
  %load ([subname,'.SCT.mat']);
  
  [ SCTright, SCTatlasright ] ...
      = cortical_SCT_measures_simulation ([subname,'.right.mid.cortex.svreg.dfs'], [pth,'/atlas.right.mid.cortex.svreg.dfs'], NLevels, attrib);
  
  save ([subname,'.SCT_simulation_',num2str(modify),'.mat'], 'SCTright', 'SCTatlasright');
  
  return
  