function Tsquared_all = hotelling_t2_test_hist_kernel_all (groupP_all, groupQ_all)
  
% groupP_all is a matrix: (number of histogram bins, num of subjects in group 1, number of mesh nodes)
% groupQ_all is a matrix: (number of histogram bins, num of subjects in group 2, number of mesh nodes)
  
  flagDebug = 0;
  
  % put histograms on hypersphere
  
  groupP_all_sphere = sqrt (groupP_all); clear groupP_all
  groupQ_all_sphere = sqrt (groupQ_all); clear groupQ_all
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  % compute Gram matrices at all nodes
  
  numOfMeshNodes = size (groupP_all_sphere, 3);
  numOfSubjectsP = size (groupP_all_sphere, 2);
  numOfSubjectsQ = size (groupQ_all_sphere, 2);
  
  % rng (0); groupP_all = mvnrnd ([0 0]', [1 0; 0 2], 500)'; groupQ_all = mvnrnd ([2 2]', [2 0; 0 1], 1000)';
  
  % Gram matrix, using vMF kernel
  conc = 1.5; % 3; % parameter
  subsetSize = 256;
  numOfSubsets = ceil (numOfMeshNodes / subsetSize);
  %
  for jj = 1 : numOfSubsets
      indicesVertices = [ 1+(jj-1)*subsetSize : min(jj*subsetSize,numOfMeshNodes) ];
      groupP_all_Gram (:,:,indicesVertices) = exp (conc * multiprod (permute (groupP_all_sphere (:,:,indicesVertices), [2 1 3]), groupP_all_sphere (:,:,indicesVertices)));
  end
  %
  for jj = 1 : numOfSubsets
      indicesVertices = [ 1+(jj-1)*subsetSize : min(jj*subsetSize,numOfMeshNodes) ];
      groupQ_all_Gram (:,:,indicesVertices) = exp (conc * multiprod (permute (groupQ_all_sphere (:,:,indicesVertices), [2 1 3]), groupQ_all_sphere (:,:,indicesVertices)));
  end
  %
  for jj = 1 : numOfSubsets
      indicesVertices = [ 1+(jj-1)*subsetSize : min(jj*subsetSize,numOfMeshNodes) ];
      groupPQ_all_Gram (:,:,indicesVertices) = exp (conc * multiprod (permute (groupP_all_sphere (:,:,indicesVertices), [2 1 3]), groupQ_all_sphere (:,:,indicesVertices)));
  end
  % takes a lot of memory
  %  groupP_all_Gram  = exp (conc * multiprod (permute (groupP_all_sphere, [2 1 3]), groupP_all_sphere));
  %  groupQ_all_Gram  = exp (conc * multiprod (permute (groupQ_all_sphere, [2 1 3]), groupQ_all_sphere));
  %  groupPQ_all_Gram = exp (conc * multiprod (permute (groupP_all_sphere, [2 1 3]), groupQ_all_sphere));
  
  %
  %groupP_all_Gram  = multiprod (permute (groupP_all_sphere, [2 1 3]), groupP_all_sphere);
  %groupQ_all_Gram  = multiprod (permute (groupQ_all_sphere, [2 1 3]), groupQ_all_sphere);
  %groupPQ_all_Gram = multiprod (permute (groupP_all_sphere, [2 1 3]), groupQ_all_sphere);
  %
  clear groupP_all_sphere
  clear groupQ_all_sphere
  
  % centering
  for i = 1 : numOfMeshNodes
    meanRowsP = mean (groupP_all_Gram (:,:,i), 2);
    meanRowsQ = mean (groupQ_all_Gram (:,:,i), 2);
    meanColsP = mean (groupP_all_Gram (:,:,i), 1);
    meanColsQ = mean (groupQ_all_Gram (:,:,i), 1);
    meanAllP  = mean (meanColsP);
    meanAllQ  = mean (meanColsQ);
    % subtract mean of rows
    groupP_all_Gram_centered (:,:,i) = bsxfun (@minus, groupP_all_Gram          (:,:,i) , meanRowsP);
    groupQ_all_Gram_centered (:,:,i) = bsxfun (@minus, groupQ_all_Gram          (:,:,i) , meanRowsQ);
    % subtract mean of columns
    groupP_all_Gram_centered (:,:,i) = bsxfun (@minus, groupP_all_Gram_centered (:,:,i), meanColsP);
    groupQ_all_Gram_centered (:,:,i) = bsxfun (@minus, groupQ_all_Gram_centered (:,:,i), meanColsQ);
    % add mean of Gram matrix
    groupP_all_Gram_centered (:,:,i) = groupP_all_Gram_centered (:,:,i) + meanAllP;
    groupQ_all_Gram_centered (:,:,i) = groupQ_all_Gram_centered (:,:,i) + meanAllQ;
  end
  
  % compute eigen analysis for all Gram matrices
  topKP = min (20, numOfSubjectsP-1); % 2, 6, 20
  topKQ = min (20, numOfSubjectsQ-1); % 2, 6, 20
  opts.issym  = true;
  opts.isreal = true;
  
  for i = 1 : numOfMeshNodes
%       if mod (i,100) == 0
%           i
%       end
    %%% group P
    [ alpha lambda ] = eigs (double (groupP_all_Gram_centered (:,:,i)) / numOfSubjectsP, topKP, 'LM', opts);
    alpha  = single (real (alpha));
    lambda = single (real (lambda));
    % normalize eigenvectors
    squaredNorms = diag (alpha' * groupP_all_Gram_centered (:,:,i) * alpha);
    alpha = bsxfun (@times, alpha, 1 ./ max (1e-5, (sqrt (squaredNorms'))));
    % save
    groupP_all_eigVecAlpha (:,:,i) = alpha;
    groupP_all_eigValues     (:,i) = diag (lambda);
    
    %%% group Q
    [ alpha lambda ] = eigs (double (groupQ_all_Gram_centered (:,:,i)) / numOfSubjectsQ, topKQ, 'LM', opts);
    alpha  = single (real (alpha));
    lambda = single (real (lambda));
    % normalize eigenvectors
    squaredNorms = diag (alpha' * groupQ_all_Gram_centered (:,:,i) * alpha);
    alpha = bsxfun (@times, alpha, 1 ./ max (1e-5, (sqrt (squaredNorms'))));
    % save
    groupQ_all_eigVecAlpha (:,:,i) = alpha;
    groupQ_all_eigValues     (:,i) = diag (lambda);
  end
  clear groupP_all_Gram_centered;
  clear groupQ_all_Gram_centered;
  
  % compute test statistic
  betaP = ones (numOfSubjectsP, 1, 'single') / numOfSubjectsP;
  betaQ = ones (numOfSubjectsQ, 1, 'single') / numOfSubjectsQ;
  for i = 1 : numOfMeshNodes
    G_P  = groupP_all_Gram  (:,:,i);
    G_PQ = groupPQ_all_Gram (:,:,i);
    %
    term_1 = (betaP' * G_P   - mean (G_P (:))) * groupP_all_eigVecAlpha (:,:,i);
    term_2 = (betaQ' * G_PQ' - mean (G_PQ(:))) * groupP_all_eigVecAlpha (:,:,i);
    term = term_1 - term_2;
    %
    LambdaP = diag (1 ./ max (1e-3, groupP_all_eigValues (:,i)));
    %
    Tsquared_all (i,1) = term * LambdaP * term';
    
    G_Q  = groupQ_all_Gram  (:,:,i);
    %
    term_1 = (betaQ' * G_Q  - mean (G_Q (:))) * groupQ_all_eigVecAlpha (:,:,i);
    term_2 = (betaP' * G_PQ - mean (G_PQ(:))) * groupQ_all_eigVecAlpha (:,:,i);
    term = term_1 - term_2;
    %
    LambdaQ = diag (1 ./ max (1e-3, groupQ_all_eigValues (:,i)));
    %
    Tsquared_all (i,1) = Tsquared_all (i,1) + term * LambdaQ * term';
  end
  
  return
  