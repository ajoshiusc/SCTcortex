function expT1_P1_all = expMapHypersphere_all (p1_all, t1_all)
  
% t1 is a column vector: (dimension,1)
% p1 is a column vector: (dimension, 1)
  
% t1Norm = norm (t1);
% expT1_P1 = cos (t1Norm) * p1 + sin (t1Norm) * t1 / max (1e-10, t1Norm);


% t1_all is a matrix: (number of histogram bins X number of mesh nodes)
% p1_all is a matrix: (number of histogram bins X number of mesh nodes)

t1Norm_all = sqrt (sum (t1_all .^ 2));
% t1Norm_all is a matrix: (1 X number of mesh nodes)

term1 = bsxfun (@times, p1_all, cos (t1Norm_all));

term2 = bsxfun (@times, t1_all, sin (t1Norm_all) ./ max (1e-10, t1Norm_all));

expT1_P1_all = term1 + term2;

return