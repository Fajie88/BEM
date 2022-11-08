
%  ========================================================================
%      BEM for the 2D Laplace equation
%                  ---- Constant element (geometry; physics)
%                  ---- Gauss integration 
% 
%      Author:  Fajie Wang          All Rights Reserved. 
%      Date: 2019-12-22
%  --------------------------------------------------------------------------
%  Qingdao University 
%  National Engineering Research Center for Intelligent Electrical Vehicle Power System  
%  E-mail: wfj1218@126.com
%  ========================================================================
clear; clc; format short e;
t0 = clock;   % Time Start  

iboundary=2;   % boundary geometry: 1-circle; 2-irregular circle; 2-rectangle
GN=8;          % number of Gaussian integrating points
[GP,GW] = GaussPoints(GN);
GP=GP'; GW=GW';

% ====== Boundary knots and boundary nodes
if iboundary==1   %---circle
nb=1000;
r=2.0;thetab=linspace(0,2*pi,nb+1);
x=r.*cos(thetab);y=r.*sin(thetab);      % Boundary knots 
xp=(x(1:nb)+x(2:nb+1))/2; yp=(y(1:nb)+y(2:nb+1))/2;        % Boundary nodes    
elseif iboundary==2   %---irregular circle
nb=1000;
thetab=linspace(0,2*pi,nb+1);
r=1.0*(cos(4*thetab)+sqrt(18/5-sin(4*thetab).^2)).^(1/3);
x=r.*cos(thetab);y=r.*sin(thetab);      % Boundary knots  
xp=(x(1:nb)+x(2:nb+1))/2; yp=(y(1:nb)+y(2:nb+1))/2;        % Boundary nodes  
elseif iboundary==3   %---rectangle
x_scale=[-1 1]; y_scale=[-1 1];
nx=401; ny=nx; nb=2*(nx+ny)-4; 
xx=linspace(x_scale(1),x_scale(2),nx); yy=linspace(y_scale(1),y_scale(2),ny);
x=[xx x_scale(2)*ones(1,ny-1) fliplr(xx(1:nx-1)) x_scale(1)*ones(1,ny-2)];
y=[y_scale(1)*ones(1,nx) yy(2:ny) y_scale(2)*ones(1,nx-1) fliplr(yy(2:ny-1))];
x(nb+1)=x(1); y(nb+1)=y(1);       % bottom,right,top,left
xp=(x(1:nb)+x(2:nb+1))/2;yp=(y(1:nb)+y(2:nb+1))/2;   
end
x=x'; y=y'; xp=xp';yp=yp';
plot(x,y,'r*',xp,yp,'b.');

% ======= Interior points (for testing)
ni=50;
theta=linspace(0,2*pi,ni);
xi=0.5*cos(theta);yi=0.5*sin(theta); xi=xi'; yi=yi';
% figure(1);h=fill(xp,yp,'k'); hold on; plot(xp,yp,'r.',xi,yi,'b.');legend('Domain', 'Bounddary nodes','Test points');
% set(h,'edgealpha',0,'facealpha',0.3)

% ======== Exact solution & Boundary condition
%u_exact = @(x,y) cos(x).*cosh(y)+sin(x).*sinh(y);
u_exact = @(x,y) cos(x).*cosh(y)+sin(x).*sinh(y)+exp(x).*cos(y)+exp(y).*sin(x)+x.^2-y.^2+2*x+3*y+1;
b = u_exact(xp, yp);

% ======== Matrix computation of BEM
a1=(x(2:nb+1)-x(1:nb))/2;a2=(x(2:nb+1)+x(1:nb))/2;
b1=(y(2:nb+1)-y(1:nb))/2;b2=(y(2:nb+1)+y(1:nb))/2;
yacb=sqrt((x(2)-x(1))^2+(y(2)-y(1))^2)/2;   % Jacobi
xk=a1*GP+a2; yk=b1*GP+b2;
xk=permute(xk, [3, 1, 2]); yk=permute(yk, [3, 1, 2]);
GW=permute(GW, [1, 3, 2]);
BEM=log(sqrt((xp - xk).^2+(yp - yk).^2));
BEM=GW.*BEM;
A=-1.0/(2*pi)*yacb.*sum(BEM,3);
A(logical(eye(size(A))))=yacb/pi*(1-log(yacb));  %奇异积分解析计算(矩阵对角线)
coeff=A\b;

% ===== Calculate interior points
BEMi=log(sqrt((xi - xk).^2+(yi - yk).^2));
BEMi=GW.*BEMi;
BEMi=-1.0/(2*pi)*yacb.*sum(BEMi,3);
u_BEM=BEMi*coeff; 

% ===== Calculate Errors 
err=abs((u_exact(xi, yi)-u_BEM)./u_exact(xi, yi));
error_max = max(abs(u_exact(xi, yi)-u_BEM)); % Maximum absolue errors
error_global = sqrt(sum((u_exact(xi, yi)-u_BEM).^2)/sum(u_exact(xi, yi).^2)); % Global error

% ===== Figures 
figure(2)
subplot(1,2,1)
plot(theta,u_exact(xi, yi),theta,u_BEM,'o');title('BEM solution');legend('Exact','BEM');xlabel('\theta');ylabel('u')
subplot(1,2,2)
plot(theta,err,'-.');title('Absolute error');xlabel('\theta');ylabel('Absolute error')

% ===== Output
disp('---- BEM -  Laplace equation ----')
disp(['Boundary element No.: ',num2str(nb)])
disp(['Tested nodes No.: ',num2str(ni)])
disp(['Elapsed total_time is: ',num2str(etime(clock,t0)),' seconds.'])
disp(['Max_Error: ',num2str(error_max)]);
disp(['Global_Error: ',num2str(error_global)]);
 