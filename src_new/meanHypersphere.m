function mu = meanHypersphere (p)

% p is a matrix: (dimension, numOfPoints)
% mu is the frechet mean: (dimension, 1)

% intiialize mean
mu = mean (p, 2);
mu = mu / max (1e-5, norm (mu));

maxIter = 20;
for iterNum = 1 : maxIter
   % negative gradient
   t = logMapHypersphere (mu, p);
   avgT = mean (t, 2);
   % gradient descent update
   muNew = expMapHypersphere (mu, 0.5 * avgT);
   % stopping criterion
   if norm (muNew - mu) < 5e-3
       mu = muNew;
       break
   end
   mu = muNew;
end
if iterNum >= maxIter
    '!!! maximum iterations reached in computing frechet mean !!!'
end

return