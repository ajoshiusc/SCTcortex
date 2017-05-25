function [ valuesPredicted predictionErrors predictionRMSE conc ] = kernel_regression_hyperSphere (hists, values)

% independent variable = hists (each column is a hist)
% dependent variable = values (row vector)

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

histsCross = hists' * hists;
numOfHists = size (hists,2);

flagFindOptBandwidth = 1;
if flagFindOptBandwidth == 1
    % tune concentration parameter using leave-one-out cross validation
    concValues = [ [0:1:29] [30:2:60] ]; % 0 : 0.01 : 1; % use a better heuristic to choose the search domain
    valuesPredictedTemp = nan (numOfHists, length (concValues));
    for i = 1 : length (concValues)
        concTemp = concValues (i);
        % perform regression using the parameter value
        weights = exp (concTemp * histsCross);
        % leave one out
        % weights = weights - diag (diag (weights));
        weights (1 : numOfHists+1 : numOfHists*numOfHists) = 0;
        % normalize weights
        weights = bsxfun (@times, weights, 1./ max (1e-10, sum (weights,2)));
        % regress
        valuesPredictedTemp (:,i) = weights * values';
    end
    % find the optmial concentration parameter
    predictionErrorsRMSDtemp = sqrt(sum ((bsxfun (@minus, valuesPredictedTemp, values')).^2) / numOfHists);
    [ minError minErrorIndex ] = min (predictionErrorsRMSDtemp);
    % optimal conc parameter
    conc = concValues (minErrorIndex);
else
    conc = 50
end

% predict again (without the leave-on-out style)
weights = exp (conc * histsCross);
weights = bsxfun (@times, weights, 1./ max (1e-10, sum (weights,2)));
% optimal predictions for the optimal conc
valuesPredicted = weights * values';

% optimal prediction errors for the optimal conc
predictionErrors = valuesPredicted - values';
predictionRMSE = sqrt (sum (predictionErrors.^2) / numOfHists);

return