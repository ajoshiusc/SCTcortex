function Tsquared_all = hotelling_t2_test_hist_all (groupP_all, groupQ_all)

% groupP_all is a matrix: (number of histogram bins, num of subjects in group 1, number of mesh nodes)
% groupQ_all is a matrix: (number of histogram bins, num of subjects in group 2, number of mesh nodes)

flagDebug = 0;

groupP_all_sphere = sqrt (groupP_all); clear groupP_all
groupQ_all_sphere = sqrt (groupQ_all); clear groupQ_all

% compute means for both groups, at each mesh node
muP_all = meanHypersphere_all (groupP_all_sphere);
muQ_all = meanHypersphere_all (groupQ_all_sphere);

% debug
%Tsquared_all = squeeze (sqrt (sum ( (muP_all - muQ_all).^2 )))';
%return

% compute test statistics
logQ_muP_all = logMapHypersphere_all (muP_all, groupQ_all_sphere);
logP_muQ_all = logMapHypersphere_all (muQ_all, groupP_all_sphere);
% logQ_muP_all is a matrix: (number of histogram bins, num of subjects in group 2, number of mesh nodes)
% logP_muQ_all is a matrix: (number of histogram bins, num of subjects in group 1, number of mesh nodes)
%
%covP_inv_all = multinv (covP_all);
%covQ_inv_all = multinv (covQ_all);
%
%mahalaQ_P_all = squeeze (mean (sum (logQ_muP_all .* multiprod (covP_inv_all, logQ_muP_all), 1), 2));
%mahalaP_Q_all = squeeze (mean (sum (logP_muQ_all .* multiprod (covQ_inv_all, logP_muQ_all), 1), 2));

if flagDebug == 1
    % for debugging
    logP_muP_all = logMapHypersphere_all (muP_all, groupP_all_sphere);
    logQ_muQ_all = logMapHypersphere_all (muQ_all, groupQ_all_sphere);
    tar_orig=readdfs('C:\Users\ajoshi\Documents\coding_ground\svreg-matlab/BrainSuiteAtlas1\mri.right.mid.cortex.dfs');
    tar=reducepatch(tar_orig,10000);
    [~,ia,indx_reduced]=intersect(tar.vertices,tar_orig.vertices,'rows','stable');
    tarsm=smooth_cortex_fast(tar,.5,200);
    %
    mean_hist_diff = squeeze (sqrt (sum ( (mean(groupP_all_sphere,2) - mean(groupQ_all_sphere,2)).^2 )));
    patch('vertices',tarsm.vertices,'faces',tar.faces,'facevertexcdata',mean_hist_diff,'edgecolor','none','facecolor','interp');axis equal;axis off;camlight;material dull;lighting phong;colormap jet; colorbar;
    close
    %
    diff = squeeze (sqrt (sum ( (muP_all - muQ_all).^2 )));
    patch('vertices',tarsm.vertices,'faces',tar.faces,'facevertexcdata',diff','edgecolor','none','facecolor','interp');axis equal;axis off;camlight;material dull;lighting phong;colormap jet; colorbar;
    close
    %
    diff = squeeze (sqrt (sum ( (mean(logP_muQ_all,2) - mean(logQ_muQ_all,2)).^2 )));
    patch('vertices',tarsm.vertices,'faces',tar.faces,'facevertexcdata',diff,'edgecolor','none','facecolor','interp');axis equal;axis off;camlight;material dull;lighting phong;colormap jet; colorbar;
    close
    %
    diff = squeeze (sqrt (sum ( (mean(logP_muP_all,2) - mean(logQ_muP_all,2)).^2 )));
    patch('vertices',tarsm.vertices,'faces',tar.faces,'facevertexcdata',diff,'edgecolor','none','facecolor','interp');axis equal;axis off;camlight;material dull;lighting phong;colormap jet; colorbar;
    close
end

% compute covariances for both groups, at each mesh node
covP_all = covHypersphere_all (muP_all, groupP_all_sphere);
covQ_all = covHypersphere_all (muQ_all, groupQ_all_sphere);
% covP_all is a matrix: (number of histogram bins X number of histogram bins X number of mesh nodes)
% covQ_all is a matrix: (number of histogram bins X number of histogram bins X number of mesh nodes)

% for debugging
if flagDebug == 1
    traceP_all = squeeze (sum (sum (covP_all.^2, 1), 2));
    traceQ_all = squeeze (sum (sum (covQ_all.^2, 1), 2));
    patch('vertices',tarsm.vertices,'faces',tar.faces,'facevertexcdata',traceP_all,'edgecolor','none','facecolor','interp');axis equal;axis off;camlight;material dull;lighting phong;colormap jet; colorbar;
    close
    patch('vertices',tarsm.vertices,'faces',tar.faces,'facevertexcdata',traceQ_all,'edgecolor','none','facecolor','interp');axis equal;axis off;camlight;material dull;lighting phong;colormap jet; colorbar;
    close
end

if flagDebug == 1
    reducedDimension = 3;
    covP_all_reduced = zeros (size (covP_all), 'single');
    covQ_all_reduced = zeros (size (covQ_all), 'single');
    for i=1:size(covP_all,3)
        [V D] = eigs (double (covP_all (:,:,i)), reducedDimension);
        covP_all_reduced (:,:,i) = single (V * D * V');
        [V D] = eigs (double (covQ_all (:,:,i)), reducedDimension);
        covQ_all_reduced (:,:,i) = single (V * D * V');
    end
    covP_all_reduced = bsxfun (@plus, covP_all_reduced, 1e-3 * eye (size (covP_all_reduced, 1)));
    covQ_all_reduced = bsxfun (@plus, covQ_all_reduced, 1e-3 * eye (size (covQ_all_reduced, 1)));
end

mahalaQ_P_all = squeeze (mean (sum (logQ_muP_all .* multi_backslash (covP_all, logQ_muP_all), 1), 2));
mahalaP_Q_all = squeeze (mean (sum (logP_muQ_all .* multi_backslash (covQ_all, logP_muQ_all), 1), 2));
if flagDebug == 1
    % for debugging
    patch('vertices',tarsm.vertices,'faces',tar.faces,'facevertexcdata',mahalaQ_P_all,'edgecolor','none','facecolor','interp');axis equal;axis off;camlight;material dull;lighting phong;colormap jet; colorbar;
    close
    patch('vertices',tarsm.vertices,'faces',tar.faces,'facevertexcdata',mahalaP_Q_all,'edgecolor','none','facecolor','interp');axis equal;axis off;camlight;material dull;lighting phong;colormap jet; colorbar;
    close
end

% mahalaQ_P_all is a vector of length: number of mesh nodes
% mahalaP_Q_all is a vector of length: number of mesh nodes
%
Tsquared_all = (mahalaQ_P_all + mahalaP_Q_all) / 2;

return