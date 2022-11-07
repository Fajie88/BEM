%  ==========================================
%      Artificial Neural Networks for Solving ODE and PDE
%      <<<<<<<<<<<<<<   coupled ODEs   >>>>>>>>>>>>>>>>
%           Problem 4
%           
%-------------------------------------------------------------------------
%         Author:  Fajie Wang          All Rights Reserved. 
%         Qingdao University             
%         wfj1218@126.com 
%         2021-07-12
%===========================================
clc; clear; format long
t0 = clock; % Time Start
global x m A1 A2 f1_right f2_right f1_phi1 f2_phi1 f1_phi2 f2_phi2
n=41;
m=10;
iexample=1;
[x_left, x_right, A1, A2, u1_exact, f1_right, u2_exact, f2_right, f1_phi1, f2_phi1, f1_phi2, f2_phi2] = example(iexample);
x=linspace(x_left,x_right,n); x=x';
a0=rand(6*m,1)+0.1; 
[a,val,k]=bfgs('Loss','Grad',a0);  %bfgs
v1=a(1:m); 
u1=a(m+1:2*m); u1=u1';
w1=a(2*m+1:3*m); w1=w1';
v2=a(3*m+1:4*m); 
u2=a(4*m+1:5*m); u2=u2';
w2=a(5*m+1:6*m); w2=w2';

sig1=logsig(x*w1+u1);
N1=sig1*v1;
N1_x=((sig1.*(1-sig1)).*w1)*v1;
sig2=logsig(x*w2+u2);
N2=sig2*v2;
N2_x=((sig2.*(1-sig2)).*w2)*v2;
u1n=A1+(x-x(1)).*N1;
u2n=A2+(x-x(1)).*N2;
u1e=u1_exact(x);
u2e=u2_exact(x);
%=================================
error1_max=max(abs(u1e-u1n));
error1_ave=sum(abs(u1e-u1n))/n;
error2_max=max(abs(u2e-u2n));
error2_ave=sum(abs(u2e-u2n))/n;
figure(1)
subplot(1,2,1), plot(x,u1e,'b-', x,u1n,'ro'); xlabel('x'); ylabel('solution'); legend('Exact', 'ANN', 'Location', 'best')
subplot(1,2,2), plot(x,u2e,'b-', x,u2n,'ro'); xlabel('x'); ylabel('solution'); legend('Exact', 'ANN', 'Location', 'best')
figure(2)
subplot(1,2,1), plot(x,abs(u1e-u1n)); xlabel('x'); ylabel('absolute error')
subplot(1,2,2), plot(x,abs(u2e-u2n)); xlabel('x'); ylabel('absolute error')
%======================Output========
disp(['Total number of points: ',num2str(n)])
disp(['Hiden units: ',num2str(m)])
disp(['Elapsed total_time is: ',num2str(etime(clock,t0)),' seconds.'])
disp(['Max_Error1: ',num2str(error1_max)]);
disp(['Ave_Error1: ',num2str(error1_ave)]);
disp(['Max_Error2: ',num2str(error2_max)]);
disp(['Ave_Error2: ',num2str(error2_ave)]);
