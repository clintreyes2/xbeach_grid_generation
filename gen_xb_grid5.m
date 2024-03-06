% generate xbeach grid from nc file 
% and generate xbeach offshore boundary file 

% take swan output and gen xb grid with some user input - no iterations
% with swan runs
% 
% shows bathymetry 
%

% version 5:
% includes function getgridinput - user input to determine rotation angle
% and x, y distances of xbeach grid without rerunning SWAN 
% extract offshore boundary points, (less complicated than previous approach)
% ===============================================================

% % include openearthtools
% pwd1=pwd;
% cd '/home/clint/xbeach/Libraries/openearthtools/matlab'; % openearthtools directory 
% addpathfast(pwd);
% cd(pwd1)

%% Load data 
nc = '../Job.nc'; % from QGIS
outdir='../XB3'; % xb grid file output directory 
gridname='xb_3m_E'; % xb grid output filename

Z=ncread(nc,'Band1');
Zmin=-30; 

x=ncread(nc,'x');
y=ncread(nc,'y');
X=meshgrid(x,y)';
Y=meshgrid(y,x);

% Create xb grid 
res_x=10;
res_y=10;

fswan2d='swan_out.mat'; % 2d swan output filename - need H, MWD (Dir), and bathymetry information
[xori, yori, dist_x, dist_y, rot]=getgridinput(fswan2d);

% write grid information to file 
fid = fopen('xb_grid_params.dat','w');
fprintf(fid,'x0\ty0\tdist_x\tdist_y\trot\n');
fprintf(fid,'%f\t%f\t%f\t%f\t%f\n', xori, yori, dist_x, dist_y, rot);
fclose(fid); 

% rotation angle must be determined by user based on swan wave output 
% xori= 636500; 
% yori= 2367000;
% dist_x=6000;
% dist_y=9000;
% rot= 47;

ix2=631; % x index location of interest (site loc)

xx=xori:res_x:xori+dist_x; 
yy=yori:res_y:yori+dist_y;

% rotate the grid 
DEP=rotateGrid(xx,yy,xori,yori,rot,X,Y,Z,Zmin);

% get the transect depths 
dep_c=DEP(:,ix2);

% set depth limits for grid determination 

% xx=xori:res_x:xori+dist_x; 
dep_lims=[10 5 3 1 -3 -5 -10];
dxs=     [10 5 3 1 3  5  10];

% find indices
i1=zeros(length(dep_lims),1); 

for i=1:length(dep_lims)
    [~,i1(i)]=min(abs(-dep_lims(i)-dep_c));
end

%% interpolate at each grid division 
% based on main transect ( transect at region of interest - our site) 
% i0 = 1; 

% setup x coordinates
for i=1:length(dep_lims)
    if i==1
        xx1=xori:dxs(i):xori+dxs(i)*(i1(i)-2);
%         dep1=dep_c(1:i1(i));

%     elseif i==length(dep_lims)
%         x0=xx1(end)+dxs(i-1);
%         xx1=[ xx1 x0:dxs(i):x0+dxs(i)*((i1(i)-i1(1))*dxs(i-1)/dxs(i)-2) ];
    else
        x0=xx1(end)+dxs(i-1);
        xx1=[ xx1 x0:dxs(i):x0+dxs(i)*((i1(i)-i1(1))*dxs(i-1)/dxs(i)-2) ];
        % no need to interp every step, just build xfinal 
%         dep2=interp1() 
%         depf=[ dep1 dep1(end)]
    end
end

dep1=interp1(xx, dep_c, xx1);

%% plot - interpolated values at transect 
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
[DEP,Xqrs,Yqrs,Xq,Yq,ZZ,R]=rotateGrid(xx1,yy,xori,yori,rot,X,Y,Z,Zmin);


%% Get offshore locations for SWAN
% v5 - changed to just the offshore boundary 

% n=50; % no of points
% % dy = dist_y/n; % spacing of boundary points
% dy = floor(length(Xqrs(1,:))/n);
dy = 50;
ny = length(Xqrs(1,:)); 

%store x, y coordinates of boundary points 
XXX = Xqrs(1,1:dy:ny);
YYY = Yqrs(1,1:dy:ny);

% write to location file
dlmwrite([gridname,'_c.loc'],[XXX;YYY]','delimiter','\t','precision','%.3f');

%% plot location of interest (site) 
% get indices 

% override location! 
xloc1=636910.97;
yloc1=2358765.61;
% [a,ix1]=min(abs(xloc1-Xqrs(:,:)))
% [b,iy1]=min(abs(yloc1-Yqrs(:,:)))

XYq2 = [xloc1-xori yloc1-yori]; 
XY2 = XYq2 * inv(R'); 
[a,ix1] = min(abs(XY2(1) - Xq(1,:))); 
[b,ix2] = min(abs(XY2(2) - Yq(:,1))); 

%% check
figure
surf(Xqrs,Yqrs,ZZ);shading interp; view(2); clim([Zmin 0])
hold on
%xbdomain
plot3(XXX,YYY,500.*ones(size(XXX)),'*r')
axis equal
colorbar
scatter(xori,yori,'r','filled')

scatter(XXX,YYY,'r')

scatter(xloc1,yloc1,'r','filled')
% scatter(Xqrs(1,ix1),Yqrs(ix2,1),'b','filled')
% xloc1-Xqrs(1,ix1)
% yloc1-Yqrs(iy1,1)

% scatter(Xqrs(1,iy1),Yqrs(ix1,1),'b','filled')

saveas(gcf,'nesting_points2b.jpg')

%% write depth file 

DEP1=DEP;
DEP1(end+1,:)=999.*ones(1,length(DEP1(1,:)));
DEP1(:,end+1)=999.*ones(length(DEP1(:,1)),1); %999 delft3d file reqt

%% make erodible layer 
ne_DEP1=zeros(size(DEP1));
ne_DEP1(:,:)=1; 
%%

if ~exist(outdir)
    mkdir(outdir);
end
% cd(outdir);
wldep('write',[outdir '/' gridname ,'.dep'],'', -DEP1); % positive downward depth 

ne_name='reef';
% write ne_layer
wldep('write',[outdir '/' ne_name,'.dep'],'', ne_DEP1); % positive downward depth 

% write grid file (.grd file - openearth tools)
wlgrid('write','FileName',[outdir '/' gridname],'X',Xqrs,'Y',Yqrs,'AutoEnclosure');


function [DEP,Xqrs,Yqrs,Xq,Yq,ZZ,R] = rotateGrid(xx,yy,xori,yori,rot,X,Y,Z,Zmin)

% rotate at the origin xori yori
% and make grid with 1d vectors xx and yy


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

% smooth bathy
DEP = smooth2a(ZZ,1,1);

end