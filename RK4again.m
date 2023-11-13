function [xs,ys,zs,ws] = RK4again(a,b,h,x0,y0,z0,w0,funcF1,funcF2,funcF3,funcF4)

% % RUNGE-KUTTA ORDE 4 DENGAN 4 PERSAMAAN 

% % membangun fungsi F,G,H--------------------------------
syms F1(x,y) F2(x,y) F3(x,y) F4(x,y)
F1(x,y) = funcF1;
F2(x,y) = funcF2;
F3(x,y) = funcF3;
F4(x,y) = funcF4;

% % membangun array xs,ys,zs,ws------------------------------
n = (b-a)/h;
xs = zeros(1,n+1); xs(1) = x0;
ys = zeros(1,n+1); ys(1) = y0;
zs = zeros(1,n+1); zs(1) = z0;
ws = zeros(1,n+1); ws(1) = w0;

% loop
for i = 2:n+1
    
    % konstanta
    K1x = F1( xs(i-1) , ys(i-1) );
    K1y = F2( xs(i-1) , ys(i-1) );
    K1z = F3( xs(i-1) , ys(i-1) );
    K1w = F4( xs(i-1) , ys(i-1) );

    K2x = F1( xs(i-1)+K1x/2*h , ys(i-1)+K1y/2*h );
    K2y = F2( xs(i-1)+K1x/2*h , ys(i-1)+K1y/2*h );
    K2z = F3( xs(i-1)+K1x/2*h , ys(i-1)+K1y/2*h );
    K2w = F4( xs(i-1)+K1x/2*h , ys(i-1)+K1y/2*h );

    K3x = F1( xs(i-1)+K2x/2*h , ys(i-1)+K2y/2*h );
    K3y = F2( xs(i-1)+K2x/2*h , ys(i-1)+K2y/2*h );
    K3z = F3( xs(i-1)+K2x/2*h , ys(i-1)+K2y/2*h );
    K3w = F4( xs(i-1)+K2x/2*h , ys(i-1)+K2y/2*h );

    K4x = F1( xs(i-1)+K3x*h , ys(i-1)+K3y*h );
    K4y = F2( xs(i-1)+K3x*h , ys(i-1)+K3y*h );
    K4z = F3( xs(i-1)+K3x*h , ys(i-1)+K3y*h );
    K4w = F4( xs(i-1)+K3x*h , ys(i-1)+K3y*h );


    % x,y,z,w baru
    xs(i) = xs(i-1) + h/6 * ( K1x + 2*K2x + 2*K3x + K4x );
    ys(i) = ys(i-1) + h/6 * ( K1y + 2*K2y + 2*K3y + K4y );
    zs(i) = zs(i-1) + h/6 * ( K1z + 2*K2z + 2*K3z + K4z );
    ws(i) = ws(i-1) + h/6 * ( K1w + 2*K2w + 2*K3w + K4w );

end



end