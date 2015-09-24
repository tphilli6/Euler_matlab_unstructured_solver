function func = setup_mms_crossterm(mms)

% Supersonic manufactured solution

l = 1;

switch mms
    case(0) 
        rho = [1,zeros(1,9)];
        u = [1,zeros(1,9)];
        v = [1,zeros(1,9)];
        p = [1,zeros(1,9)];
        
        func.rho = @(x) base_function(x, rho, l);
        func.u = @(x) base_function(x, u, l);
        func.v = @(x) base_function(x, v, l);
        func.p = @(x) base_function(x, p, l);

    case (1)
    rho = [1, ...
         -0.15, 2, 1/3, ...
         0.1, 1, -1/5, ...
         -0.05, 1, 1/4];

    u = [800, ...
         50, 2/3, 1/2, ...
         30, 1.5, -1/8, ...
         -60, 1.5, 1/8];

    v = [800, ...
         -75, 5/3, -1/2, ...
         40, 3, 0, ...
         60, 1.5, 1/8];

    p = [100000, ...
         20000, 1/2, 1/2, ...
         50000, 1, -1/2, ...
         -10000, 3/4, -1/4];
     
        func.rho = @(x) base_function(x, rho, l);
        func.u = @(x) base_function(x, u, l);
        func.v = @(x) base_function(x, v, l);
        func.p = @(x) base_function(x, p, l);

    case (2)
      rho = [1, ...
           0.15, 2, 1/3, ...
           -0.1, 1, -1/5, ...
           -0.05, 1, 1/4];

      u = [80, ...
           5, 2/3, 1/2, ...
           -3, 1.5, -1/8, ...
           -6, 1.5, 1/8];

      v = [80, ...
           -7.5, 5/3, -1/2, ...
           4.0, 3, 0, ...
           6.0, 1.5, 1/8];

      p = [100000, ...
           -20000, 1/2, 1/2, ...
           50000, 1, -1/2, ...
           -10000, 3/4, -1/4];

    func.rho = @(x) base_function(x, rho, l);
    func.u = @(x) base_function(x, u, l);
    func.v = @(x) base_function(x, v, l);
    func.p = @(x) base_function(x, p, l);
    
    case(3)
        rhoinf = 1;
        pinf = 100000;
        minf = 0.5;
        vinv = sqrt(1.4*pinf/rhoinf)*minf;
        alpha = 0;
        mu = complex(-0.0475,0);
        te_angle = -2*(0.594689181*(1/2*0.298222773/sqrt(1) - 0.127125232 - 2*0.357907906*1 + 3*0.291984971*1*1 - 4*0.105174606*1*1*1));
        
      func.rho = @(x) karman_trefftz_airfoil(x, 1, rhoinf, pinf, vinv, alpha, mu, te_angle);
      func.u   = @(x) karman_trefftz_airfoil(x, 2, rhoinf, pinf, vinv, alpha, mu, te_angle);
      func.v   = @(x) karman_trefftz_airfoil(x, 3, rhoinf, pinf, vinv, alpha, mu, te_angle);
      func.p   = @(x) karman_trefftz_airfoil(x, 4, rhoinf, pinf, vinv, alpha, mu, te_angle);
end

    func.intx = @(x, c) intx_base_function(x, c, l);
    func.cex = [rho; u; v; p];

    func.drdx = @(x) bf_x(x, rho, l);
    func.drdy = @(x) bf_y(x, rho, l);
    
    func.drdxx = @(x) bf_xx(x, rho, l);
    func.drdxy = @(x) bf_xy(x, rho, l);
    func.drdyy = @(x) bf_yy(x, rho, l);
    
    func.drdxxx = @(x) bf_xxx(x, rho, l);
    func.drdxxy = @(x) bf_xxy(x, rho, l);
    func.drdxyy = @(x) bf_xyy(x, rho, l);
    func.drdyyy = @(x) bf_yyy(x, rho, l);
    
    func.drdxxxx = @(x) bf_xxxx(x, rho, l);
    func.drdxxxy = @(x) bf_xxxy(x, rho, l);
    func.drdxxyy = @(x) bf_xxyy(x, rho, l);
    func.drdxyyy = @(x) bf_xyyy(x, rho, l);
    func.drdyyyy = @(x) bf_yyyy(x, rho, l);
    
    
    func.dudx = @(x) bf_x(x, u, l);
    func.dudy = @(x) bf_y(x, u, l);
    
    func.dudxx = @(x) bf_xx(x, u, l);
    func.dudxy = @(x) bf_xy(x, u, l);
    func.dudyy = @(x) bf_yy(x, u, l);
    
    func.dudxxx = @(x) bf_xxx(x, u, l);
    func.dudxxy = @(x) bf_xxy(x, u, l);
    func.dudxyy = @(x) bf_xyy(x, u, l);
    func.dudyyy = @(x) bf_yyy(x, u, l);
    
    func.dudxxxx = @(x) bf_xxxx(x, u, l);
    func.dudxxxy = @(x) bf_xxxy(x, u, l);
    func.dudxxyy = @(x) bf_xxyy(x, u, l);
    func.dudxyyy = @(x) bf_xyyy(x, u, l);
    func.dudyyyy = @(x) bf_yyyy(x, u, l);
    
    
    func.dvdx = @(x) bf_x(x, v, l);
    func.dvdy = @(x) bf_y(x, v, l);
    
    func.dvdxx = @(x) bf_xx(x, v, l);
    func.dvdxy = @(x) bf_xy(x, v, l);
    func.dvdyy = @(x) bf_yy(x, v, l);
    
    func.dvdxxx = @(x) bf_xxx(x, v, l);
    func.dvdxxy = @(x) bf_xxy(x, v, l);
    func.dvdxyy = @(x) bf_xyy(x, v, l);
    func.dvdyyy = @(x) bf_yyy(x, v, l);
    
    func.dvdxxxx = @(x) bf_xxxx(x, v, l);
    func.dvdxxxy = @(x) bf_xxxy(x, v, l);
    func.dvdxxyy = @(x) bf_xxyy(x, v, l);
    func.dvdxyyy = @(x) bf_xyyy(x, v, l);
    func.dvdyyyy = @(x) bf_yyyy(x, v, l);
    
    
    func.dpdx = @(x) bf_x(x, p, l);
    func.dpdy = @(x) bf_y(x, p, l);
    
    func.dpdxx = @(x) bf_xx(x, p, l);
    func.dpdxy = @(x) bf_xy(x, p, l);
    func.dpdyy = @(x) bf_yy(x, p, l);
    
    func.dpdxxx = @(x) bf_xxx(x, p, l);
    func.dpdxxy = @(x) bf_xxy(x, p, l);
    func.dpdxyy = @(x) bf_xyy(x, p, l);
    func.dpdyyy = @(x) bf_yyy(x, p, l);
    
    func.dpdxxxx = @(x) bf_xxxx(x, p, l);
    func.dpdxxxy = @(x) bf_xxxy(x, p, l);
    func.dpdxxyy = @(x) bf_xxyy(x, p, l);
    func.dpdxyyy = @(x) bf_xyyy(x, p, l);
    func.dpdyyyy = @(x) bf_yyyy(x, p, l);
    
end

function f = base_function(x,c,l)

% p=4;
% x0 = [0.5,0.5];
% f = c(1) + c(2)*(x(1)-x0(1)).^p + c(5)*(x(2)-x0(2)).^p;% + c(8)*x(1)*x(2);
% f = c(1) + c(2)*x(1)*x(2) + c(5)*x(1).^2*x(2).^2 + c(8)*x(1).^3*x(2);% + c(8)*x(1)*x(2);

f = c(1) + c(2)*sin(c(3)*pi*x(1)/l + c(4)*pi) ...
        + c(5)*sin(c(6)*pi*x(2)/l + c(7)*pi);% ...
%         + c(8)*sin(c(9)*pi*x(1)*x(2)/l^2 + c(10)*pi);
end

% For exact solution calculation
function f = intx_base_function(x,c,l)

% p=4;
% x0 = [0.5,0.5];
% f = c(1).*x(1) + c(2)*(x(1)-x0(1)).^(p+1)/(p+1) + c(5)*(x(2)-x0(2)).^p.*x(1);% + c(8)*x(1)*x(2);
% f = c(1) + c(2)*x(1)*x(2) + c(5)*x(1).^2*x(2).^2 + c(8)*x(1).^3*x(2);% + c(8)*x(1)*x(2);

f = c(1).*x(1) - c(2)*l/(pi*c(3))*cos(c(3)*pi*x(1)/l + c(4)*pi) ...
        + c(5)*sin(c(6)*pi*x(2)/l + c(7)*pi).*x(1);% ...
%         - c(8)*l^2./(pi*c(9)*x(2)).*cos(c(9)*pi*x(1)*x(2)/l^2 + c(10)*pi);
    
    
end

%df/dx
function f = bf_x(x,c,l)
    f = (c(3)*pi/l).*c(2)*cos(c(3)*pi*x(1)/l + c(4)*pi) ...
    + (c(9)*pi*x(2)/l^2).*c(8)*cos(c(9)*pi*x(1)*x(2)/l^2 + c(10)*pi);
end

%df/dy
function f = bf_y(x,c,l)

f = (c(6)*pi/l).*c(5)*cos(c(6)*pi*x(2)/l + c(7)*pi) ...
        + (c(9)*pi*x(1)/l^2).*c(8)*cos(c(9)*pi*x(1)*x(2)/l^2 + c(10)*pi);
end





%d2f/dx2
function f = bf_xx(x,c,l)
    f = -(c(3)*pi/l).^2.*c(2)*sin(c(3)*pi*x(1)/l + c(4)*pi) ...
    - (c(9)*pi*x(2)/l^2).^2.*c(8)*sin(c(9)*pi*x(1)*x(2)/l^2 + c(10)*pi);
end

%d2f/dxdy
function f = bf_xy(x,c,l)
    f = (c(9)*pi/l^2).*c(8)*cos(c(9)*pi*x(1)*x(2)/l^2 + c(10)*pi)...
       - x(1).*x(2).*(c(9)*pi/l^2).^2.*c(8)*sin(c(9)*pi*x(1)*x(2)/l^2 + c(10)*pi);
end

%df/dy
function f = bf_yy(x,c,l)

f = -(c(6)*pi/l).^2.*c(5)*sin(c(6)*pi*x(2)/l + c(7)*pi) ...
        - (c(9)*pi*x(1)/l^2).^2.*c(8)*sin(c(9)*pi*x(1)*x(2)/l^2 + c(10)*pi);
end





%d3f/dx3
function f = bf_xxx(x,c,l)
    f = -(c(3)*pi/l).^3.*c(2)*cos(c(3)*pi*x(1)/l + c(4)*pi) ...
    - x(2).^3*(c(9)*pi/l^2).^3.*c(8)*cos(c(9)*pi*x(1)*x(2)/l^2 + c(10)*pi);
end

%d3f/dx2dy
function f = bf_xxy(x,c,l)
    a = c(8);
    b = c(9)*pi/l^2;
    f = -2*x(2).*b.^2.*a*sin(b*x(1)*x(2) + c(10)*pi)...
        -x(1).*x(2).^2.*b.^3.*a*cos(b*x(1)*x(2) + c(10)*pi);
end

%d3f/dxdy2
function f = bf_xyy(x,c,l)
    a = c(8);
    b = c(9)*pi/l^2;
    f = -2.*x(1)*b.^2.*a*sin(b*x(1)*x(2) + c(10)*pi)...
        -x(1).^2.*x(2).*b.^3.*a*cos(b*x(1)*x(2) + c(10)*pi);
end

%d3f/dy3
function f = bf_yyy(x,c,l)

f = -(c(6)*pi/l).^3.*c(5)*cos(c(6)*pi*x(2)/l + c(7)*pi) ...
        - (c(9)*pi*x(1)/l^2).^3.*c(8)*cos(c(9)*pi*x(1)*x(2)/l^2 + c(10)*pi);
end





%d4f/dx4
function f = bf_xxxx(x,c,l)
    f = (c(3)*pi/l).^4.*c(2)*sin(c(3)*pi*x(1)/l + c(4)*pi) ...
      + (c(9)*pi*x(2)/l^2).^4.*c(8)*sin(c(9)*pi*x(1)*x(2)/l^2 + c(10)*pi);
end

%d4f/dx3dy
function f = bf_xxxy(x,c,l)
  a = c(8);
  b = c(9)*pi/l^2;
%         f = - x(2).^3*b.^3.*a*cos( b*x(1)*x(2) + c(10)*pi);
%         f = -3*x(2).^2*b.^3.*a*cos(b*x(1)*x(2) + c(10)*pi)...
%             +x(1)*x(2).^3*b.^4*a*sin(b*x(1)*x(2) + c(10)*pi);
    f = a*b^4*x(2)^3*x(1)*sin(b*x(1)*x(2)+c(10)*pi)...
  - 3*a*b^3*x(2)^2*cos(b*x(1)*x(2)+c(10)*pi);
end

%d4f/dx2dy2
function f = bf_xxyy(x,c,l)
    a = c(8);
    b = c(9)*pi/l^2;
%     f = (c(9)*pi*x(2)/l^2).^2.*(c(9)*pi*x(1)/l^2).^2.*c(8)*sin(c(9)*pi*x(1)*x(2)/l^2 + c(10)*pi);
    f = a*b^4*x(1)^2*x(2)^2*sin(b*x(1)*x(2)+c(10)*pi)...
      - 4*a*b^3*x(1)*x(2)*cos(b*x(1)*x(2)+c(10)*pi)...
      - 2*a*b^2*sin(b*x(1)*x(2)+c(10)*pi);
end

%d4f/dxdy3
function f = bf_xyyy(x,c,l)
    a = c(8);
    b = c(9)*pi/l^2;
%     f = (c(9)*pi*x(2)/l^2).*(c(9)*pi*x(1)/l^2).^3.*c(8)*sin(c(9)*pi*x(1)*x(2)/l^2 + c(10)*pi);
    f = a*b^4*x(1)^3*x(2)*sin(b*x(1)*x(2)+c(10)*pi)...
      - 3*a*b^3*x(1)^2*cos(b*x(1)*x(2)+c(10)*pi);
end

%d4f/dy4
function f = bf_yyyy(x,c,l)
  f = (c(6)*pi/l).^4.*c(5)*sin(c(6)*pi*x(2)/l + c(7)*pi) ...
    + (c(9)*pi*x(1)/l^2).^4.*c(8)*sin(c(9)*pi*x(1)*x(2)/l^2 + c(10)*pi);
end


function derivative_tests
    l = 1;
        
    coef = [1, ...
         -0.15, 2, 1/3, ...
         0.1, 1, -1/5, ...
         -0.05, 1, 1/4];
     
    f = @(x) base_function(x, coef, l);
    
    x =[0.1,-0.1];
    dx0 = 1;
    
    dfdx(1,:) = bf_x(x,coef,l);
    dfdx(2,:) = bf_y(x,coef,l);
    
    dfdx(3,:) = bf_xx(x,coef,l);
    dfdx(4,:) = bf_xy(x,coef,l);
    dfdx(5,:) = bf_yy(x,coef,l);
    
    dfdx(6,:) = bf_xxx(x,coef,l);
    dfdx(7,:) = bf_xxy(x,coef,l);
    dfdx(8,:) = bf_xyy(x,coef,l);
    dfdx(9,:) = bf_yyy(x,coef,l);
    
    dfdx(10,:) = bf_xxxx(x,coef,l);
    dfdx(11,:) = bf_xxxy(x,coef,l);
    dfdx(12,:) = bf_xxyy(x,coef,l);
    dfdx(13,:) = bf_xyyy(x,coef,l);
    dfdx(14,:) = bf_yyyy(x,coef,l);
    
    for i = 1:10;
      dx = [dx0,0]; dfdx_test(1,i) = (f(x+dx)-f(x-dx))./(2*dx0);
      dx = [0,dx0]; dfdx_test(2,i) = (f(x+dx)-f(x-dx))./(2*dx0);
        
      dx = [dx0,0]; dfdx_test(3,i) = (f(x+dx)-2*f(x)+f(x-dx))./(dx0.^2);
        dfdx_test(4,i) = (f(x+[dx0,dx0])-f(x+[dx0,-dx0]) - f(x+[-dx0,dx0]) + f(x+[-dx0,-dx0]))./(4*dx0^2);
      dx = [0,dx0]; dfdx_test(5,i) = (f(x+dx)-2*f(x)+f(x-dx))./(dx0.^2);
      
      dx = [dx0,0]; dfdx_test(6,i) = (-1/2*f(x-2*dx) + f(x-dx) - f(x+dx) +1/2*f(x+2*dx))./(dx0^3);
        d2fdx2_p1 = bf_xx(x+[0,dx0],coef,l);
        d2fdx2_m1 = bf_xx(x+[0,-dx0],coef,l);
      dfdx_test(7,i) = (d2fdx2_p1-d2fdx2_m1)/(2*dx0);
        d2fdx2_p1 = bf_yy(x+[dx0,0],coef,l);
        d2fdx2_m1 = bf_yy(x+[-dx0,0],coef,l);
      dfdx_test(8,i) = (d2fdx2_p1-d2fdx2_m1)/(2*dx0);
      dx = [0,dx0]; dfdx_test(9,i) = (-1/2*f(x-2*dx) + f(x-dx) - f(x+dx) +1/2*f(x+2*dx))./(dx0^3);

        d2fdx2_p1 = bf_xxx(x+[dx0,0],coef,l);
        d2fdx2_m1 = bf_xxx(x+[-dx0,0],coef,l);
      dfdx_test(10,i) = (d2fdx2_p1-d2fdx2_m1)/(2*dx0);
        d2fdx2_p1 = bf_xxx(x+[0,dx0],coef,l);
        d2fdx2_m1 = bf_xxx(x+[0,-dx0],coef,l);
      dfdx_test(11,i) = (d2fdx2_p1-d2fdx2_m1)/(2*dx0);
        d2fdx2_p1 = bf_xx(x+[0,dx0],coef,l);
        d2fdx2_n  = bf_xx(x,coef,l);
        d2fdx2_m1 = bf_xx(x+[0,-dx0],coef,l);
      dfdx_test(12,i) = (d2fdx2_p1 - 2*d2fdx2_n + d2fdx2_m1)/(dx0^2);
        d2fdx2_p1 = bf_yyy(x+[dx0,0],coef,l);
        d2fdx2_m1 = bf_yyy(x+[-dx0,0],coef,l);
      dfdx_test(13,i) = (d2fdx2_p1-d2fdx2_m1)/(2*dx0);
        d2fdx2_p1 = bf_yyy(x+[0,dx0],coef,l);
        d2fdx2_m1 = bf_yyy(x+[0,-dx0],coef,l);
      dfdx_test(14,i) = (d2fdx2_p1-d2fdx2_m1)/(2*dx0);
        
      dx0=dx0/2;
    end
        
    err = dfdx_test - repmat(dfdx,[1,size(dfdx_test,2)]);
    
    p = log( err(:,1:end-1)./err(:,2:end) )/log(2);
    
    r=fliplr(0:size(p,2)-1);
    h = 2.^r;
    
    semilogx(h,p,'-o');
    legend('dfdx','df/dy','d2f/dx2','d2f/dxdy','d2f/dy2',...
           'd3f/dx3','d3f/dx2dy','d3f/dxdy2','d3f/dy3',...
           'd4f/dx4','d4f/dx3dy','d4f/dx2dy2','d4f/dxdy3','d4f/dy4');
    
end
