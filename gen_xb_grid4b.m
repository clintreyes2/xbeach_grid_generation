% generate xbeach grid from nc file 
% and generate xbeach offshore boundary file 
% note: that varying grid is based on one transect of choice, usually at
% project site 

% ===================== note on previous data =====================
% previous grid was overwritten, 
% to recover old data - rerun with 3c version ! 
% ===============================================================

% added from 3c, we add breawater bathymetry using a subroutine - 
% first test is no subroutine (gen_xb_grid4)

% % include openearthtools
pwd1=pwd;
cd '/home/clint/xbeach/OET/matlab'; % openearthtools directory 
addpathfast(pwd);
cd(pwd1)

%% Load data 
nc = '../Job.nc'; % from QGIS
gridname='xb_3m_E';

Z=ncread(nc,'Band1');
Zmin=-30; 

x=ncread(nc,'x');
y=ncread(nc,'y');
X=meshgrid(x,y)';
Y=meshgrid(y,x);

% Create xb grid 
% xori= 637152;
% yori= 2358338;
% xori= 639473;
% yori= 2360800;
xori= 632000; %637101, 2364610
yori= 2369000;
res_x=10; % offshore resolution 
res_y=10;
dist_x=10000;
dist_y=9000;

rot= 0;

ix2=631;

xx=xori:res_x:xori+dist_x; 
yy=yori:res_y:yori+dist_y;

DEP=rotateGrid(xx,yy,xori,yori,rot,X,Y,Z,Zmin);

%
% get the transect depths 
dep_c=DEP(:,ix2);

% set depth limits for grid determination 

% xx=xori:res_x:xori+dist_x; 
dep_lims=[20 10 5 3];
dxs=     [10 5 3 1];

% find indices
i1=zeros(length(dep_lims),1); 

for i=1:length(dep_lims)
    [~,i1(i)]=min(abs(-dep_lims(i)-dep_c));
end

%% interpolate at each grid division - make into a subroutine 
% i0 = 1; 
for i=1:length(dep_lims)
    if i==1
        xx1=xori:dxs(i):xori+dxs(i)*(i1(i)-2);
%         dep1=dep_c(1:i1(i));

    elseif i==length(dep_lims)
        x0=xx1(end)+dxs(i-1);
        % xx1=[ xx1 x0:dxs(i):x0+dxs(i)*((i1(i)-i1(1))*dxs(i-1)/dxs(i)-2) ];
        xx1=[ xx1 x0:dxs(i):xori+dist_x ];
    else
        x0=xx1(end)+dxs(i-1);
        xx1=[ xx1 x0:dxs(i):x0+dxs(i)*((i1(i)-i1(1))*dxs(i-1)/dxs(i)-2) ];
        % no need to interp every step, just build xfinal 
%         dep2=interp1() 
%         depf=[ dep1 dep1(end)]
    end
end

% xx1(end)
dep_c(end)
dep1=interp1(xx, dep_c, xx1);

%% plot transect - with new varied spacing 
figure 
plot(xx,dep_c)
ylim([-30,5])
hold on 
scatter(xx(i1),dep_c(i1))
scatter(xx(i1),-dep_lims,'r+')
grid on 

plot(xx1,dep1,'ko','markersize',2)
saveas(gcf, 'transect1.png')

%% obtain new rotated grid 
[~,Xqrs,Yqrs,Xq,Yq,ZZ,R]=rotateGrid(xx1,yy,xori,yori,rot,X,Y,Z,Zmin);
% %% check difference of smoothed and non-smoothed 
% figure
% imagesc(DEP-ZZ)
% impixelinfo
% colorbar

% !!!!!! updates with gen_xb_grid4b from v4 !!!!!!
% ==== remove/ignore DEP(smoothed ZZ) from rotateGrid ====

% ==== take ZZ, convert to ZZ with writeOverBathy ====

% plot the (unrotated) depth 
figure 
imagesc(ZZ)
hold on 
scatter(ix2, ix1, 'r', 'filled')
axis on
colorbar
clim([-15, 0])
impixelinfo
title('ZZ')

% ============ take user input for Breakwater origin ============

% [p1x,p1y]=ginput(1);
% %scatter(p1x,p1y,'r')
% p1x = round(p1x); 
% p1y = round(p1y);

%override for now
p1x = 600; % indices
p1y = 840; 
% scatter(p1x,p1y,10,'r+')

% other user input 
% currently in grid spacing, convert to meters later 

L = 20; % breakwater length (extend to the right)  
b = 7; % breakwater width (extend shoreward) 
ht = 1; % height of breakwater from toe (match elev with initial point - for a constant crest structure) 

%  ============ end user input ====================================

ZZ2 = writeOverBathy(p1x, p1y, L, b, ht, ZZ); 

scatter(px,py,10,'r+')

% plots to check new bathymetry 
figure 
imagesc(ZZ2-ZZ)
colorbar
impixelinfo
title('ZZ2-ZZ')

figure 
imagesc(ZZ2)
colorbar
impixelinfo
clim([-15, 0])
title('ZZ2')

% ==== apply smooth2a to obtain DEP ====
DEP = smooth2a(ZZ2, 1, 1);
% check difference of smoothed and non-smoothed 
figure
imagesc(DEP-ZZ2)
impixelinfo
colorbar
title('DEP-ZZ2')

% plot DEP
figure 
imagesc(DEP)
colorbar
impixelinfo
clim([-15, 0])
title('DEP')


%%  take a transect and plot 1d 
figure 
plot(830:870, DEP(830:870,610))
hold on 
plot(830:870, ZZ2(830:870,610))
plot(830:870, ZZ(830:870,610))
legend('smoothed', 'wBW','orig')

% DEP is then output to a file ...

%%  ----- new in version 4 ---- 
% add breakwater - rectangular section 
% input taken: rectangular crest coordinates 
%              height of structure from the seabed 
% smooth with bed using a smoothing function - smooth2a 

% crest coordinates: 
% pof1 = ; % offshore point 1
% pof2 = ; % offshore point 2  
% pon1 = ; % onshore point 1 
% pon2 = ; % onshore point 2

% try to do by selecting, imagesc? 

%% Get offshore locations for SWAN
% plot nesting points 

n=50; %spacing, no of points

XXX=zeros(1,size(Xqrs,2));
YYY=zeros(1,size(Yqrs,2));

for cc=1:size(Yqrs,2)
    for rr=1:size(Xqrs,1)
        if ZZ(rr,cc)>Zmin
            XXX(cc)=Xqrs(rr,cc); 
            YYY(cc)=Yqrs(rr,cc);
            break
%         else
        end
    end
end

% make to SWAN resolution
XXX=XXX(1:n:end);
YYY=YYY(1:n:end);

% cleanup unpopulated areas
XXX(XXX==0)=[];
YYY(YYY==0)=[];

length(XXX)

% [] other points, too shallow/land regions

% write to location file
dlmwrite([gridname,'_c.loc'],[XXX;YYY]','delimiter','\t','precision','%.3f');

% !! location to plot !! not related to transect ! be careful 
% this is just a reference for the rough location of the site 
xloc1=636910.97;
yloc1=2358765.61;
% [a,ix1]=min(abs(xloc1-Xqrs(:,:)));
% [b,iy1]=min(abs(yloc1-Yqrs(:,:)));

XYq2 = [xloc1-xori yloc1-yori]; 
XY2 = XYq2 * inv(R'); 
[a,ix1] = min(abs(XY2(1) - Xq(1,:))); 
[b,ix2] = min(abs(XY2(2) - Yq(:,1))); 

% check
figure
surf(Xqrs,Yqrs,ZZ);shading interp; view(2); clim([Zmin 0])
hold on
%xbdomain
plot3(XXX,YYY,500.*ones(size(XXX)),'*r')
axis equal
colorbar
scatter(xori,yori,'r','filled')

scatter(xloc1,yloc1,'r','filled')
% scatter(Xqrs(1,ix1),Yqrs(ix2,1),'b','filled')
% xloc1-Xqrs(1,ix1)
% yloc1-Yqrs(iy1,1)

% scatter(Xqrs(1,iy1),Yqrs(ix1,1),'b','filled')

saveas(gcf,'nesting_points2b.jpg')
% 

%%
% % figure 
% % surf(ZZ,'FaceAlpha',0.85);shading interp; view(2); clim([Zmin 0])
% % % hold on
% % set(gca, 'ydir', 'normal')
% % set(gca, 'xdir', 'reverse')
% % xlim([0,length(Xqrs)])
% % ylim([0,length(Xqrs)])
% % grid on
% 
%% write depth file 

DEP1=DEP;
DEP1(end+1,:)=999.*ones(1,length(DEP1(1,:)));
DEP1(:,end+1)=999.*ones(length(DEP1(:,1)),1); %999 delft3d file reqt

%%
% figure 
% imagesc(flipud(DEP1))
% impixelinfo
% colorbar
% clim([-30 0])
%%
% DEP1=flipud(DEP1);
% 
% %% make erodible layer 
% ne_DEP1=zeros(size(DEP1));
% ne_DEP1(:,:)=1; 
% %%
% 
outdir='../XB2';
if ~exist(outdir)
    mkdir(outdir);
end
% cd(outdir);
wldep('write',[outdir '/' gridname ,'.dep'],'', -DEP1); % positive downward depth 

% write non erodible depth file 
% ne_name='reef';
% % write ne_layer
% wldep('write',[ne_name,'.dep'],'', ne_DEP1); % positive downward depth 

% write grid file (.grd file - openearth tools)
wlgrid('write','FileName',[outdir '/' gridname],'X',Xqrs,'Y',Yqrs,'AutoEnclosure');

%
fprintf(['grid created in ' outdir '!\n'] )
%% 
function [DEP,Xqrs,Yqrs,Xq,Yq,ZZ,R] = rotateGrid(xx,yy,xori,yori,rot,X,Y,Z,Zmin)

[XX,YY]=meshgrid(xx,yy); % use this for plotting ? 

Xq=XX-xori;
Yq=YY-yori;
XY=[Xq(:) Yq(:)]; %create matrix of vectors
theta=-(rot-270); %to rotate clockwise by x degrees
R=[cosd(theta) -sind(theta); sind(theta) cosd(theta)]; %create the matrix
rotXY=XY*R'; %multiply vectors by the rot matrix
Xqr=reshape(rotXY(:,1),size(Xq,1),[]);
Yqr=reshape(rotXY(:,2),size(Yq,1),[]);
%shifting
Xqrs=Xqr+xori;
Yqrs=Yqr+yori;
Xqrs=Xqrs';
Yqrs=Yqrs';

%get boundiung box
boun_x=[Xqrs(1,1); Xqrs(1,end); Xqrs(end,end); Xqrs(end,1); Xqrs(1,1)];
boun_y=[Yqrs(1,1); Yqrs(1,end); Yqrs(end,end); Yqrs(end,1); Yqrs(1,1)];

%check for shadow zone
%plot xb grid over Hs

% interpolate
F=griddedInterpolant(X,Y,Z);
ZZ=F(Xqrs,Yqrs);
ZZ(ZZ<=Zmin)=Zmin;
%replot

% !!!!!!! might be better to write over before rotateGrid function !!!!!!!!!!!!!!!!!

% smooth bathy
DEP = smooth2a(ZZ,1,1);

end