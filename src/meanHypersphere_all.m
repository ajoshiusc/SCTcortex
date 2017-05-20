function mu_all = meanHypersphere_all (p_all)

% p_all is a matrix : (number of histogram bins, num of subjects in group 1, number of mesh nodes)
% mu_all is a matrix: (number of histogram bins, number of mesh nodes)

% intiialize means
mu_all = squeeze (mean (p_all, 2));
mu_norm_all = sqrt (sum (mu_all.^2, 1));
mu_all = bsxfun (@times, mu_all, 1 ./ max (1e-5, mu_norm_all));

maxIter = 20;
for iterNum = 1 : maxIter
%     '-----------------------'
%     iterNum
   % negative gradient
   t_all = logMapHypersphere_all (mu_all, p_all);
   
%    t1 = logMapHypersphere (mu_all (:,1), p_all (:,:,1));
%    t2 = logMapHypersphere (mu_all (:,2), p_all (:,:,2));
%    t1 (:,1:5)
%    t_all (:,1:5,1)
%    t2 (:,1:5)
%    t_all (:,1:5,2)
%    pause
   
   avgT_all = squeeze (mean (t_all, 2));
   % gradient descent update
   muNew_all = expMapHypersphere_all (mu_all, 0.5 * avgT_all);
   % stopping criterion
   if norm (muNew_all - mu_all, 'fro') / size (mu_all, 2) < 1e-3
       mu_all = muNew_all;
       break
   end
   mu_all = muNew_all;
%    pause
end
if iterNum >= maxIter
    '!!! maximum iterations reached in computing frechet mean !!!'
end

return