clear all
close all
rng (0);

flagShow = 0

numOfHistBins = 125
numOfSubjects = 80
numOfMeshNodes = 5e3

'create groupA ...'
groupA_all= zeros (numOfHistBins, numOfSubjects, numOfMeshNodes);
for n = 1 : numOfMeshNodes
    muRand = ones (numOfHistBins, 1) / sqrt (numOfHistBins); muRand (1) = muRand (1) + 0.5 * rand; muRand = muRand / norm (muRand);
    a = max (0, bsxfun (@plus, 0.1 * randn (numOfHistBins, numOfSubjects), muRand)); a = bsxfun (@times, a, 1 ./ sqrt (sum (a.^2)));
    groupA_all (:,:,n) = a;
end

'create groupB ...'
groupB_all= zeros (numOfHistBins, numOfSubjects, numOfMeshNodes);
for n = 1 : numOfMeshNodes
    muRand = ones (numOfHistBins, 1) / sqrt (numOfHistBins); muRand (2) = muRand (2) - 0.5 * rand; muRand = muRand / norm (muRand);
    b = max (0, bsxfun (@plus, 0.1 * randn (numOfHistBins, numOfSubjects), muRand)); b = bsxfun (@times, b, 1 ./ sqrt (sum (b.^2)));
    groupB_all (:,:,n) = b;
end

if flagShow == 1
    figure; hold on;
    n1 = 1
    n2 = 2
    plot3 (groupA_all(1,:,n1), groupA_all(2,:,n1), groupA_all(3,:,n1), 'ro', groupA_all(1,:,n2), groupA_all(2,:,n2), groupA_all(3,:,n2), 'mo');
    plot3 (groupB_all(1,:,n1), groupB_all(2,:,n1), groupB_all(3,:,n1), 'bo', groupB_all(1,:,n2), groupB_all(2,:,n2), groupB_all(3,:,n2), 'co');
    axis equal tight; grid on; pause, close
end

'Hotelling test ...'
profile off; profile on;
%tic
Tsquared_all = hotelling_t2_test_hist_all (groupA_all, groupB_all);
%toc
profile viewer;
