function [x_left, x_right, u1_x0, u2_x0, u1_exact, f1_right, u2_exact, f2_right, f1_u1, f2_u1, f1_u2, f2_u2] = example(iexample)
if iexample==1
    x_left=0;  x_right=3;
    u1_x0=0;  u2_x0=1;
    u1_exact=@(x) sin(x);
    u2_exact=@(x) 1+x.^2;
    f1_right=@(x ,u1, u2) cos(x)+u1.^2+u2-(1+x.^2+sin(x).*sin(x));
    f2_right=@(x, u1, u2) 2*x-(1+x.^2).*sin(x)+u1.*u2;
    f1_u1 = @ (u1, u2) 2*u1;
    f2_u1 = @ (u1, u2) u2;
    f1_u2 = @ (u1, u2) 1;
    f2_u2 = @ (u1, u2) u1;
elseif iexample==2
    x_left=0;  x_right=3;
    u1_x0=0;  u2_x0=1;
    u1_exact=@(x) sin(x);
    u2_exact=@(x) 1+x.^2;
    f1_right=@(x ,u1, u2) cos(x)+u1.^2+u2-(1+x.^2+sin(x).*sin(x));
    f2_right=@(x, u1, u2) 2*x-(1+x.^2).*sin(x)+u1.*u2;
elseif iexample==3
    x_left=0;  x_right=3;
    u1_x0=0;  u2_x0=1;
    u1_exact=@(x) sin(x);
    u2_exact=@(x) 1+x.^2;
    f1_right=@(x ,u1, u2) cos(x)+u1.^2+u2-(1+x.^2+sin(x).*sin(x));
    f2_right=@(x, u1, u2) 2*x-(1+x.^2).*sin(x)+u1.*u2;
elseif iexample==4
    x_left=0;  x_right=3;
    u1_x0=0;  u2_x0=1;
    u1_exact=@(x) sin(x);
    u2_exact=@(x) 1+x.^2;
    f1_right=@(x ,u1, u2) cos(x)+u1.^2+u2-(1+x.^2+sin(x).*sin(x));
    f2_right=@(x, u1, u2) 2*x-(1+x.^2).*sin(x)+u1.*u2;
end
end

