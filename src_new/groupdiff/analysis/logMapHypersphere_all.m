function logP2_P1_all = logMapHypersphere_all (p1_all, p2_all)

% p1 is a column vector: (dimension,1)
% p2 is a matrix: (dimension, numOfPoints)

% innerProduct = max (-1, min (1, p2' * p1));
% 
% logP2_P1 ...
%     =  ( p2 - p1 * innerProduct' ) ...
%     .* ( ones (size (p2,1), 1, 'single') * ...
%     ( acos (innerProduct) ./ ...
%     max (1e-10, sqrt (1 - innerProduct.^2)) ...
%     )' ...
%     );

% p1_all is a matrix: (number of histogram bins, number of mesh nodes)
% p2_all is a matrix: (number of histogram bins, num of subjects, number of mesh nodes)

numOfBins      = size (p2_all, 1);
numOfSubjects  = size (p2_all, 2);
numOfMeshNodes = size (p2_all, 3);

p2_all_trans = permute (p2_all, [2,1,3]);
% p2_all_trans  is a matrix: (num of subjects, number of histogram bins, number of mesh nodes)
p1_all_reshape = reshape (p1_all, numOfBins, 1, numOfMeshNodes);
% p1_all_reshape is a matrix: (number of histogram bins, 1, number of mesh nodes)
innerProduct_all = max (-1, min (1, (multiprod_sparse (p2_all_trans, p1_all_reshape))));
% innerProduct_all is a matrix: (num of subjects, 1, number of mesh nodes)

innerProduct_all_permute = permute (innerProduct_all,[2,1,3]);
factor1 = p2_all - multiprod_sparse (p1_all_reshape, innerProduct_all_permute);
% numer1 is a matrix: (number of histogram bins, num of subjects, number of mesh nodes)

factor2 = acos (innerProduct_all_permute) ./ max (1e-10, sqrt (1 - innerProduct_all_permute.^2));

logP2_P1_all = factor1 .* multiprod_sparse (ones (numOfBins, 1, numOfMeshNodes), factor2);
% logP2_P1_all is a matrix: (number of histogram bins, num of subjects, number of mesh nodes)

return