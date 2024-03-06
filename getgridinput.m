function [xori, yori, dist_x, dist_y, rot] = getgridinput(fswan2d)

% this function reads swan 2d output 
% and give input params to generate xb grid 
% xori, yori, dist_x, dist_y, rot

% uncomment some lines for user input 

load(fswan2d);

%%
theta_N = Dir_00000100_000000; % MWD in degrees - nautical  
theta_cart = 90-theta_N; % degrees - cartesian 
theta_cart(theta_cart<0) = 360-theta_cart(theta_cart<0);
theta_cart(theta_cart>360) = theta_cart(theta_cart>360) - 360; 

figure('position',[500 500 500 1000]);
subplot(2,1,1)
imagesc(theta_N)
impixelinfo
colorbar
clim([0 360])
subplot(2,1,2)
imagesc(theta_cart)
impixelinfo
colorbar
clim([0 360])

%% get x and y 
qx = -cos(theta_cart*pi/180);
qy = -sin(theta_cart*pi/180);

quivspcg=20;

fsize=1000;
figure('position',[500 500 fsize fsize])
h=pcolor(Xp,Yp,Hsig_00000100_000000);
set(h, 'edgecolor', 'none')
colorbar
clim([0 4])
hold on
quiver(Xp(1:quivspcg:end,1:quivspcg:end),Yp(1:quivspcg:end,1:quivspcg:end),...
    qx(1:quivspcg:end,1:quivspcg:end),qy(1:quivspcg:end,1:quivspcg:end),'k')
contour(Xp,Yp,Botlev,[-10:10:30],'k-','ShowText','on'); 
title('Hs')

% user input: xb-grid origin (uncomment for user input)
% ups=zeros(3,2); % x, y 
% [ups(1,1),ups(1,2)]=ginput(1); 
% scatter(ups(1,1),ups(1,2),'red')
% % user input: y point 2 - to det rotation angle and y dist (xb domain)
% [ups(2,1),ups(2,2)]=ginput(1); 
% scatter(ups(2,1),ups(2,2),'red')
% % user input: x point 2 - to det x dist (xb domain)
% [ups(3,1),ups(3,2)]=ginput(1); 
% scatter(ups(3,1),ups(3,2),'red')

% override - with known points 
ups = [636010.654109589	2367857.89163498;
643672.243150685	2362305.61787072;
631321.921232877	2363103.54942966];

% write ups (user points) to file 
fileID = fopen('user_points.txt','w');
fprintf(fileID, '#x\ty\n');
for i =1:3
    fprintf(fileID, '%f\t%f\t%f\n', ups(i,1), ups(i,2)); 
end
fclose(fileID);

scatter(ups(1,1),ups(1,2),'red')
scatter(ups(2,1),ups(2,2),'red')
scatter(ups(3,1),ups(3,2),'red')
plot(ups(1:2,1),ups(1:2,2),'red')
plot(ups(1:2:3,1),ups(1:2:3,2),'red')
saveas(gcf,'xbgrid_show.png')

xori=ups(1,1); 
yori=ups(1,2); 
dist_y = dist(ups(1,1),ups(1,2),ups(2,1),ups(2,2)); %input: x1,y1,x2,y2
dist_x = dist(ups(1,1),ups(1,2),ups(3,1),ups(3,2));
% obtain rotation angle, rot
x1=ups(1,1);
x2=ups(2,1);
y1=ups(1,2);
y2=ups(2,2);
dx21 = x2-x1;
dy21 = y2-y1; 
ang = atan(abs(dy21)/abs(dx21)); 
if x2>x1 && y2<y1
    rot=ang; % rotation angle 
elseif x2<x1 && y2<y1
    rot=180-ang; 
elseif x2<x1 && y2>y1
    rot=180+ang;
else
    rot=360-ang; 
end
rot = rot/pi*180;

end

function d=dist(x1,y1,x2,y2)
    d=sqrt((x1-x2)^2 + (y1-y2)^2);
end
