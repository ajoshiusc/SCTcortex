function show_surf_val(so,val,clim,cmap)
if ~exist('clim','var')
    clim=[min(val),max(val)];
end
if ~exist('cmap','var')
    cmap=jet;
end
so.attributes=val;
s=close_surf(so);
val=s.attributes;
figure;
patch('faces',s.faces,'vertices',s.vertices,'facevertexcdata',val,'edgecolor','none','facecolor','interp');
axis equal;axis off tight;material dull;view(90,0);camlight;lighting phong;colormap(cmap); caxis(clim);
colorbar
