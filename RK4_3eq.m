function [xs,ys,zs] = RK4_3eq(a,b,h,x0,y0,z0,funcF,funcG,funcH)
% % RUNGE-KUTTA ORDE 4 DENGAN 3 PERSAMAAN 

% % membangun fungsi F,G,H--------------------------------
syms F(t,x,y,z) G(t,x,y,z) H(t,x,y,z)
F(t,x,y,z) = funcF;
G(t,x,y,z) = funcG;
H(t,x,y,z) = funcH;

% % membangun array ts,xs,ys,zs--------------------------------
n = (b-a)/h;
ts = zeros(1,n+1) + a;
xs = zeros(1,n+1); xs(1) = x0;
ys = zeros(1,n+1); ys(1) = y0;
zs = zeros(1,n+1); zs(1) = z0;

% loop
for i = 2:n+1
    ts(i) = ts(i) + (i-1)*h;                      % t baru
    
    % konstanta
    K1x = F( ts(i-1) , xs(i-1) , ys(i-1) , zs(i-1) );
    K1y = G( ts(i-1) , xs(i-1) , ys(i-1) , zs(i-1) );
    K1z = H( ts(i-1) , xs(i-1) , ys(i-1) , zs(i-1) );

    K2x = F( ts(i-1)+h/2 , xs(i-1)+K1x/2*h , ys(i-1)+K1y/2*h , zs(i-1)+K1z/2*h );
    K2y = G( ts(i-1)+h/2 , xs(i-1)+K1x/2*h , ys(i-1)+K1y/2*h , zs(i-1)+K1z/2*h );
    K2z = H( ts(i-1)+h/2 , xs(i-1)+K1x/2*h , ys(i-1)+K1y/2*h , zs(i-1)+K1z/2*h );

    K3x = F( ts(i-1)+h/2 , xs(i-1)+K2x/2*h , ys(i-1)+K2y/2*h , zs(i-1)+K2z/2*h );
    K3y = G( ts(i-1)+h/2 , xs(i-1)+K2x/2*h , ys(i-1)+K2y/2*h , zs(i-1)+K2z/2*h );
    K3z = H( ts(i-1)+h/2 , xs(i-1)+K2x/2*h , ys(i-1)+K2y/2*h , zs(i-1)+K2z/2*h );

    K4x = F( ts(i-1)+h , xs(i-1)+K3x*h , ys(i-1)+K3y*h , zs(i-1)+K3z*h );
    K4y = G( ts(i-1)+h , xs(i-1)+K3x*h , ys(i-1)+K3y*h , zs(i-1)+K3z*h );
    K4z = H( ts(i-1)+h , xs(i-1)+K3x*h , ys(i-1)+K3y*h , zs(i-1)+K3z*h );


    % x,y,z baru
    xs(i) = xs(i-1) + h/6 * ( K1x + 2*K2x + 2*K3x + K4x );
    ys(i) = ys(i-1) + h/6 * ( K1y + 2*K2y + 2*K3y + K4y );
    zs(i) = zs(i-1) + h/6 * ( K1z + 2*K2z + 2*K3z + K4z );

end


end