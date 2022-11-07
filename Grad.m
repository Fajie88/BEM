function G = Grad(a)
global x m A1 A2 f1_right f2_right f1_phi1 f2_phi1 f1_phi2 f2_phi2
v1=a(1:m); 
u1=a(m+1:2*m); u1=u1';
w1=a(2*m+1:3*m); w1=w1';
v2=a(3*m+1:4*m); 
u2=a(4*m+1:5*m); u2=u2';
w2=a(5*m+1:6*m); w2=w2';

sig1=logsig(x*w1+u1);
sig1z=sig1.*(1-sig1);
sig1zz=sig1.*(1-sig1).*(1-sig1)-sig1.*sig1.*(1-sig1);
sig2=logsig(x*w2+u2);
sig2z=sig2.*(1-sig2);
sig2zz=sig2.*(1-sig2).*(1-sig2)-sig2.*sig2.*(1-sig2);

N1=sig1*v1;
N1_x=((sig1.*(1-sig1)).*w1)*v1;
phi1=A1+(x-x(1)).*N1;
phi1_x=N1+(x-x(1)).*N1_x;

sig2=logsig(x*w2+u2);
N2=sig2*v2;
N2_x=((sig2.*(1-sig2)).*w2)*v2;
phi2=A2+(x-x(1)).*N2;
phi2_x=N2+(x-x(1)).*N2_x;

N1_v = sig1;
N1_w = x.*sig1z.*v1';
N1_u =sig1z.*v1';

N2_v = sig2;
N2_w = x.*sig2z.*v2';
N2_u =sig2z.*v2';

N1_x_v =sig1z.*w1;
N1_x_w = x.*sig1zz.*w1.*v1'+sig1z.*v1';
N1_x_u =sig1zz.*w1.*v1';

N2_x_v =sig2z.*w2;
N2_x_w = x.*sig2zz.*w2.*v2'+sig2z.*v2';
N2_x_u =sig2zz.*w2.*v2';

phi1_v=(x-x(1)).*N1_v;
phi1_u=(x-x(1)).*N1_u;
phi1_w=(x-x(1)).*N1_w;
phi2_v=(x-x(1)).*N2_v;
phi2_u=(x-x(1)).*N2_u;
phi2_w=(x-x(1)).*N2_w;

phi1_x_v = N1_v + (x-x(1)).*N1_x_v;
phi1_x_w = N1_w + (x-x(1)).*N1_x_w;
phi1_x_u = N1_u + (x-x(1)).*N1_x_u;

phi2_x_v = N2_v + (x-x(1)).*N2_x_v;
phi2_x_w = N2_w + (x-x(1)).*N2_x_w;
phi2_x_u = N2_u + (x-x(1)).*N2_x_u;

F_v1=2*sum((phi1_x-f1_right(x,phi1,phi2)).*(phi1_x_v-f1_phi1(phi1,phi2).*phi1_v))+2*sum((phi2_x-f2_right(x,phi1,phi2)).*(-f2_phi1(phi1,phi2).*phi1_v));
F_w1=2*sum((phi1_x-f1_right(x,phi1,phi2)).*(phi1_x_w-f1_phi1(phi1,phi2).*phi1_w))+2*sum((phi2_x-f2_right(x,phi1,phi2)).*(-f2_phi1(phi1,phi2).*phi1_w));
F_u1=2*sum((phi1_x-f1_right(x,phi1,phi2)).*(phi1_x_u-f1_phi1(phi1,phi2).*phi1_u))+2*sum((phi2_x-f2_right(x,phi1,phi2)).*(-f2_phi1(phi1,phi2).*phi1_u));

F_v2=2*sum((phi2_x-f2_right(x,phi1,phi2)).*(phi2_x_v-f2_phi2(phi1,phi2).*phi2_v))+2*sum((phi1_x-f1_right(x,phi1,phi2)).*(-f1_phi2(phi1,phi2).*phi2_v));
F_w2=2*sum((phi2_x-f2_right(x,phi1,phi2)).*(phi2_x_w-f2_phi2(phi1,phi2).*phi2_w))+2*sum((phi1_x-f1_right(x,phi1,phi2)).*(-f1_phi2(phi1,phi2).*phi2_w));
F_u2=2*sum((phi2_x-f2_right(x,phi1,phi2)).*(phi2_x_u-f2_phi2(phi1,phi2).*phi2_u))+2*sum((phi1_x-f1_right(x,phi1,phi2)).*(-f1_phi2(phi1,phi2).*phi2_u));
%F=sum((phi1_x-f1_right(x,phi1)).^2)+sum((phi2_x-f2_right(x,phi2)).^2);

G= [F_v1 F_u1 F_w1 F_v2 F_u2 F_w2];
G=G';
end

