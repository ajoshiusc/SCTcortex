function expT1_P1 = expMapHypersphere (p1, t1)
  
% t1 is a column vector: (dimension,1)
% p1 is a column vector: (dimension, 1)
  
t1Norm = norm (t1);
expT1_P1 = cos (t1Norm) * p1 + sin (t1Norm) * t1 / max (1e-10, t1Norm);

return
  