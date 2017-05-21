function [ cov_all t_all ] = covHypersphere_all (mu_all, p_all)

% p is a matrix: (dimension, numOfPoints)
% mu is the frechet mean: (dimension, 1)

% t = logMapHypersphere (mu, p);
% cov = (t * t') / size (t, 2);

% p_all is a matrix: (number of histogram bins X number of subjects X number of mesh nodes)
% mu_all is the matrix: (number of histogram bins X number of mesh nodes)

t_all = logMapHypersphere_all (mu_all, p_all);
% t_all is a matrix: (number of histogram bins X number of subjects X number of mesh nodes)

% cov_all = multiprod (t_all, permute (t_all, [2 1 3])) / size (t_all, 2);
cov_all = multiprod_sparse (t_all, permute (t_all, [2 1 3])) / size (t_all, 2);
% cov_all is a matrix: (number of histogram bins X number of histogram bins X number of mesh nodes)

% regularize
cov_all = bsxfun (@plus, cov_all, 1e-3 * eye (size (cov_all, 1)));

return