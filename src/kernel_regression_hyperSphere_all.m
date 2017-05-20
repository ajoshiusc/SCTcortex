function [ valuesPredicted_all, predictionErrors_all, predictionRMSE_all, conc_all ] ...
    = kernel_regression_hyperSphere_all (histsSqrtCrossExp_all, indicesVertex, values)

% independent variable = hists (hist-dimension X number of subjects X number of vertices)
% dependent variable = values (1 X number of subjects)

% for testing
% rng (0); dim = 125; mu = rand (dim,1); mu = mu / norm (mu); 
% N = 8e2; a = bsxfun (@plus, 0.3 * randn (dim,N), mu); a = bsxfun (@times, a', 1./sqrt(sum(a.^2)'))';
% values = 1 + 10 * exp (0.3 * mu' * a); values = values - min(values); values = values / max (values); values = values + 0.2 * randn (size (values));
% values = sin (12 * a (1,:));
% [ valuesPredicted predictionErrors predictionRMSE conc ] = kernel_regression_hyperSphere (a, values); conc
% we should expect an optimal value of conc to be around 0.5
% we should expect optimal valuesPredicted to be close to 'values'
% norm (values - valuesPredicted') / norm (values)
% plot (values, valuesPredicted, 'ro'); grid on; axis equal tight; pause, close

% map hists to the unit hypersphere
% hists = max (0, hists);
% hists = bsxfun (@multiply, hists', 1./sum(hists)')';
% hists = hists.^(0.5);

% kernel = von Mises Fisher; exp (conc * x' * y)

numOfSubjects = size (histsSqrtCrossExp_all,1);
numOfVertices = length (indicesVertex);

% store the indices of the diagonal elements in histsCross_all
% will be used later in LOO cross validation for tuning 
indicesDiag = zeros (numOfSubjects * numOfVertices, 1, 'uint32');
for jj = 1 : numOfVertices
    indicesDiag (1+(jj-1)*numOfSubjects : jj*numOfSubjects) = [ 1 : numOfSubjects+1 : numOfSubjects^2 ] + (jj-1)*numOfSubjects^2;
end

flagFindOptBandwidth = 1;
if flagFindOptBandwidth == 1
    %
    weightsFactor = histsSqrtCrossExp_all (:,:,indicesVertex);
    
    % tune concentration parameter using leave-one-out cross validation
    concValues = [ [0:1:20] [22:2:40] [44:4:60] ]; % 0 : 0.01 : 1; % use a better heuristic to choose the search domain
    valuesPredictedTemp_all = nan (numOfSubjects, numOfVertices, length (concValues));
    for i = 1 : length (concValues)
        %
        concTemp = concValues (i);
        % perform regression using the parameter value
        if concTemp == 0
            weights_all = ones (size (weightsFactor), 'single');
            % leave one out
            weights_all (indicesDiag) = 0;
        else
            if concTemp == 22 | concTemp == 44
                weightsFactor = weightsFactor .^ 2;
            end
            weights_all = weights_all .* weightsFactor;
        end
        % normalize weights (columns)
        weights_all = bsxfun (@rdivide, weights_all, sum (weights_all,1));
        % regress
        valuesPredictedTemp_all (:,:,i) = squeeze (multiprod (values, weights_all));
    end
    
    % find the optmial concentration parameter
    predictionErrorsRMSDtemp_all = sqrt (sum ((bsxfun (@minus, valuesPredictedTemp_all, values')).^2) / numOfSubjects);
    [ minError_all minErrorIndex_all ] = min (predictionErrorsRMSDtemp_all, [], 3);
    clear valuesPredictedTemp_all predictionErrorsRMSDtemp_all
    % optimal conc parameter
    conc_all = concValues (minErrorIndex_all);
else
    conc = 50
end

% regularize a bit more
% conc_all = conc_all / 4;

% predict again (without the leave-one-out style)
weights_all = bsxfun (@power, histsSqrtCrossExp_all (:,:,indicesVertex), reshape (conc_all, [1 size(conc_all)] ));
weights_all = bsxfun (@rdivide, weights_all, sum (weights_all,1));

% optimal predictions for the optimal conc
valuesPredicted_all = squeeze (multiprod (values, weights_all));

% optimal prediction errors for the optimal conc
predictionErrors_all = bsxfun (@minus, valuesPredicted_all, values');
predictionRMSE_all = sqrt (sum (predictionErrors_all.^2) / numOfSubjects);

return