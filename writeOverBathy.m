function ZZ2 = writeOverBathy(p1x, p1y, L, b, ht, ZZ)

% % note: x is x axis of the image, not in xbeach (which is cross-shore) 
p2x = p1x+L;
p2y = p1y; 
p3x = p1x; 
p3y = p1y+b; 
p4x = p2x; 
p4y = p3y;
px = [p1x p2x p3x p4x];
py = [p1y p2y p3y p4y];

% figure 
% imagesc(ZZ)
% hold on 
% colorbar
% clim([-15, 0])
% impixelinfo
% scatter(px,py,10,'r+')

% depth at p1 :
p1_dep = ZZ(p1y,p1x);
c_dep = p1_dep + ht; % crest depth 

% set depths inside rectangle (if else statements)
ZZ2 = ZZ; 
ZZ2(p1y:p3y, p1x:p2x) = c_dep; 

% %% plots to check new bathymetry 
% figure 
% imagesc(ZZ2-ZZ)
% colorbar
% impixelinfo
% 
% figure 
% imagesc(ZZ2)
% colorbar
% impixelinfo
% clim([-15, 0])
