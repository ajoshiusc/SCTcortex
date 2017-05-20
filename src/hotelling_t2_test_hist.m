function Tsquared = hotelling_t2_test_hist (groupP, groupQ)

% groupP is a matrix: (dimension, numOfPoints)
% groupQ is a matrix: (dimension, numOfPoints)

muP = meanHypersphere (groupP);
muQ = meanHypersphere (groupQ);

covP = covHypersphere (muP, groupP) + 1e-15 * eye (size (groupP, 1));
covQ = covHypersphere (muQ, groupQ) + 1e-15 * eye (size (groupQ, 1));

logQ_muP = logMapHypersphere (muP, groupQ);
logP_muQ = logMapHypersphere (muQ, groupP);

mahalaQ_P = mean (sum (logQ_muP .* (covP \ logQ_muP)));
mahalaP_Q = mean (sum (logP_muQ .* (covQ \ logP_muQ)));

Tsquared = (mahalaQ_P + mahalaP_Q) / 2;

return