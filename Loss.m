function F=Loss(a)
global x m A1 A2 f1_right f2_right
v1=a(1:m); 
u1=a(m+1:2*m); u1=u1';
w1=a(2*m+1:3*m); w1=w1';
v2=a(3*m+1:4*m); 
u2=a(4*m+1:5*m); u2=u2';
w2=a(5*m+1:6*m); w2=w2';
% w=rand(1,m);
% v=rand(m,1);
% u=rand(1,m);

sig1=logsig(x*w1+u1);
N1=sig1*v1;
N1_x=((sig1.*(1-sig1)).*w1)*v1;
phi1=A1+(x-x(1)).*N1;
phi1_x=N1+(x-x(1)).*N1_x;

sig2=logsig(x*w2+u2);
N2=sig2*v2;
N2_x=((sig2.*(1-sig2)).*w2)*v2;
phi2=A2+(x-x(1)).*N2;
phi2_x=N2+(x-x(1)).*N2_x;

F=sum((phi1_x-f1_right(x,phi1,phi2)).^2)+sum((phi2_x-f2_right(x,phi1,phi2)).^2);
end
