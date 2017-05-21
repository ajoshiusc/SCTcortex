function [ cov t] = covHypersphere (mu, p)

% p is a matrix: (dimension, numOfPoints)
% mu is the frechet mean: (dimension, 1)

t = logMapHypersphere (mu, p);
cov = (t * t') / size (t, 2);

return