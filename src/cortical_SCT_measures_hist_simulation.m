function [SCT_hist,SCT_hist_atlas]= cortical_SCT_measures_hist_simulation (sub_surf,atlas_surf,NRingsNbr,attrib)

flagShow = 0;

s=readdfs(sub_surf);
s=smooth_cortex_fast(s,.1,1000);
'compute principal curvatures, mean curvature, gaussian curvature ...'
[Cmean,Cgaussian,Dir1,Dir2,Lambda1,Lambda2]=patchcurvature(s,1);
%
tmp1 = max (Lambda1, Lambda2);
tmp2 = min (Lambda1, Lambda2);
Lambda1 = tmp1;
Lambda2 = tmp2;
p=prctile(Lambda1,[1,99]);   Lambda1=max(p(1),min(Lambda1,p(2)));
p=prctile(Lambda2,[1,99]);   Lambda2=max(p(1),min(Lambda2,p(2)));
p=prctile(Cmean,[1,99]);     Cmean=max(p(1),min(Cmean,p(2)));
p=prctile(Cgaussian,[1,99]); Cgaussian=max(p(1),min(Cgaussian,p(2)));

if flagShow == 1
    a=bipolar;
    figure;
    patch('faces',s.faces,'vertices',s.vertices,'facevertexcdata',-Cmean,'edgecolor','none','facecolor','interp');
    caxis([-.25,.25]); axis equal;camlight;material dull;view(-90,0);camlight;colormap(a);colorbar;axis off;
    figure;
    patch('faces',s.faces,'vertices',s.vertices,'facevertexcdata',Cgaussian,'edgecolor','none','facecolor','interp');
    caxis([-.0015,.0015]); axis equal;camlight;material dull;view(-90,0);camlight;colormap(a);colorbar;axis off;
end

% compute shape index and curvedness
S=(2/pi)*atan((Lambda2+Lambda1)./(Lambda2-Lambda1));
if sum (isnan (S) > 0)
    error ('S has a Nan');
end
C=0.5*(Lambda2.^2+Lambda1.^2).^0.5 ;
if sum (isnan (C) > 0)
    error ('C has a Nan');
end
% get thickness values
T=s.attributes;
if sum (isnan (T) > 0)
    error ('T has a Nan');
end

if flagShow == 1
    a=bipolar;
    figure;
    patch('faces',s.faces,'vertices',s.vertices,'facevertexcdata',S,'edgecolor','none','facecolor','interp');
    caxis([-1,1]); axis equal;camlight;material dull;view(-90,0);camlight;colormap(a);colorbar;axis off;
    figure;
    patch('faces',s.faces,'vertices',s.vertices,'facevertexcdata',C,'edgecolor','none','facecolor','interp');
    caxis([0,.5]); axis equal;camlight;material dull;view(-90,0);camlight;colorbar;axis off;
end

'find neighborhoods ...'
tri=triangulation(s.faces,s.vertices);
ed=tri.edges;
Adj=sparse(ed(:,1),ed(:,2),1,length(s.vertices),length(s.vertices));
Adj = Adj+Adj';
switch NRingsNbr
    case 4
        AdjN = Adj^2 > 0;
        AdjN = AdjN^2 > 0;
    case 8
        AdjN = Adj^2 > 0;
        AdjN = AdjN^2 > 0;
        AdjN = AdjN^2 > 0;
    case 16
        AdjN = Adj^2 > 0;
        AdjN = AdjN^2 > 0;
        AdjN = AdjN^2 > 0;
        AdjN = AdjN^2 > 0;
    case 32
        AdjN = Adj^2 > 0;
        AdjN = AdjN^2 > 0;
        AdjN = AdjN^2 > 0;
        AdjN = AdjN^2 > 0;
        AdjN = AdjN^2 > 0;
    otherwise
        error ('');
end
% figure;
% patch('faces',s.faces,'vertices',s.vertices,'facevertexcdata',double(full(AdjN(:,1))),'facecolor','interp','edgecolor','none')
% axis off;axis equal;camlight;material dull;

% S=attrib.*S;
SmeanPositive = mean (S (S > 0));
SmeanNegative = mean (S (S < 0));
SmedPositive = median (S (S > 0));
SmedNegative = median (S (S < 0));

S (S > 0) = (S (S > 0) - SmedPositive) .* attrib (S > 0) + SmedPositive;
S (S < 0) = (S (S < 0) - SmedNegative) .* attrib (S < 0) + SmedNegative;

C=attrib.*C;
T=attrib.*T;
SCT=[S,C,T]';

% joint histogram of S,C,T
numOfBins = 5;
% shape indices in neighborhood
edgesS = linspace(-1,1,numOfBins+1);
% shape indices in neighborhood
minCurvedness = 0;
maxCurvedness = 0.6;
edgesC = linspace(minCurvedness,maxCurvedness,numOfBins+1);
% shape indices in neighborhood
minThickness = 1;
maxThickness = 8;
edgesT = linspace(minThickness,maxThickness,numOfBins+1);
% clamping
SCT (1,:) = max (edgesS(1)+1e-5, SCT (1,:)); SCT (1,:) = min (edgesS(end)-1e-5, SCT (1,:));
SCT (2,:) = max (edgesC(1)+1e-5, SCT (2,:)); SCT (2,:) = min (edgesC(end)-1e-5, SCT (2,:));
SCT (3,:) = max (edgesT(1)+1e-5, SCT (3,:)); SCT (3,:) = min (edgesT(end)-1e-5, SCT (3,:));
%
SCT_hist=zeros(numOfBins, numOfBins, numOfBins, length(s.vertices));
for jj=1:length(s.vertices)
    nbr = AdjN(:,jj);
    h = histcn(SCT(:,nbr)',edgesS,edgesC,edgesT);
    SCT_hist (:,:,:,jj) = h;
end
% vectorize the histograms
SCT_hist = reshape (SCT_hist, [ numOfBins^3 length(s.vertices) ]);

% read atlas
at=readdfs(atlas_surf);

% warp SCT to atlas coordinate space
SCT_hist = SCT_hist'; % dimensions: #vertices x #bins
%
SCT_hist_atlas = map_data_flatmap(s,SCT_hist,at);
% reshape: % dimensions: #bins x #vertices
SCT_hist = SCT_hist';
SCT_hist_atlas = SCT_hist_atlas';
% normalize
SCT_hist       = min (1, max (0, bsxfun (@times, SCT_hist',       1./sum(SCT_hist)')'));
SCT_hist_atlas = min (1, max (0, bsxfun (@times, SCT_hist_atlas', 1./sum(SCT_hist_atlas)')'));

return
