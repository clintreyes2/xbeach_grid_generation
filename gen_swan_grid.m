close all; clear all; clc

%% Load data
nc='Job.nc';

% Z=ncread('job.nc','Band1');
% x=ncread('job.nc','x');
% y=ncread('job.nc','y');
% [Y,X]=meshgrid(y,x);

x_temp=ncread(nc,'x');
y_temp=ncread(nc,'y');
X=meshgrid(x_temp,y_temp)';
Y=meshgrid(y_temp,x_temp);
Z=ncread(nc,'Band1');

%--
figure
surf(X,Y,Z)
view(2)
shading interp
colorbar
clim([-30 0])
% title('SWAN Domain');
% xlabel('[UTM]');
% ylabel('[UTM]');

%% Build SWAN grid
MaxY= 2355664.5+13116;
MinY= 2355664.5;
MaxX= 632329.5+14565;
MinX= 632329.5;

% x array
xx=MinX:50:MaxX;

% y array
yy=MinY:50:MaxY;

% mesh grid
[XX,YY]=meshgrid(xx,yy);

F=griddedInterpolant(X,Y,Z);
ZZ=F(XX,YY);
ZZ_min=-30; % buoy depth
ZZ(ZZ<ZZ_min)=ZZ_min;

figure()
surf(XX,YY,ZZ)
view(2);shading interp;clim([-30,0]);
% title('SWAN Domain');
% xlabel('[UTM]');
% ylabel('[UTM]');
saveas(gcf,'SWAN_Domain.jpg')

%% Convert into SWAN
Xlen=MaxX-MinX;
Ylen=MaxY-MinY;

%export
% cd('/home/charlotte/waimanalonewcode1/SWAN')
swan_io_bot('write','swan_bot.txt',ZZ')


