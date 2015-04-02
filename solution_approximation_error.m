function [soln_aprx_err] = solution_approximation_error(cell, face, vertex, ...
                        kexact_order, ...
                        kexact_type, ...
                        fit_type,...
                        flux_integral_order,...
                        analytic_soln,...
                        flux)
                    poly_test
% Solution approximation error
% Estimate the solution approximation error and return the adaption metric
% Original test is setup to estimate the solution approximation error for a
% third-order accurate solution
% 
% Following the work of Pagnutti
% The aniso metric is based off of the Hessian 
%        | d2f/dx2  d2f/dxdy |
%    H = | d2fdxdy  d2f/dy2  |
%
% This is based on a Taylor series approximation for the error in a
%  third-order solution
% error(dx, dy) = | 1/6 d3f/dx3 (dx)^3 + 1/2 d3f/dx2y (dx^2 dy) + 1/2 d3f/dxdy2 (dx dy^2) + 1/6 d3f/dy3 (dy)^3
%
% The node spacing is given by
% d(dx,dy) = sqrt( [dx,dy] M [dx; dy] )
%
% The Taylor series error is approximated by a polynomial of degree p
% i.e. P = ax^3 + bx^2y + cxy^2 + dy^3
% error(x,y) = 1/2 sum_{i=1}^n |P_i|
% i is the number of equations (i.e. L1 norm of reconstruction error)
% Normalization due to magnitude differences?
%
% P is integrated around a circle to result in a Fourier series function
% (f(theta))^(2/3) = a0 + a2 cos(2theta) + b2 sin(2theta)
%
% Where the adaption metric is 
%           | a0/2 + a2  b2        |
%      M  = | b2         a0/2 - a2 |
%
%  And a0 is limited to a0 = max( a0, 2sqrt( a2^2 + b2^2) ) so that M is
%  positive definate

  
  % Normalize solution
  rho_inf=1;
  u_inf = 800;
  v_inf = 800;
  p_inf = 100000;
  cell.soln(:,1) = cell.soln(:,1)/rho_inf;
  cell.soln(:,2) = cell.soln(:,2)/u_inf;
  cell.soln(:,3) = cell.soln(:,3)/v_inf;
  cell.soln(:,4) = cell.soln(:,4)/p_inf;


[regular_reconstruction, cell.lhs_set] = reconstruct_solution(cell, fit_type);
regular_reconstruction_param = cell.reconstruction_param;


eq = {'mass','xmtm', 'ymtm', 'nrgy'};
base_order = kexact_order;
base_type  = kexact_type;



  % Step 1: compute a lsq reconstruction
  % - setup_reconstruction sets up the polynomial for the reconstruction
  kexact_order = base_order + 1;
  fit_type = 'lsq';%base_fit_type;
 
  setup_reconstruction;
  higher_order_recon_param = cell.reconstruction_param;
  build_kexact_stencil;

  % ******* Important output ********* higher_order_recon
  % reconstruct the higher order solution
  % Sets up higher_order_recon(n).coef  % coefficients of higher order reconstruction

  
  
  cell.lhs_set = 0; % Reset stencil
  [higher_order_recon, ~] = reconstruct_solution(cell, fit_type);
%    cell.reconstruction = higher_order_recon;
   
   % higher order terms for k=3
   % p_high = [x^3 + x^2 y + xy^2 + y^3];
   p_high = [4,7,9,10];
   Tcoef = [1/6,1/2,1/2,1/6];
   for n = 1:cell.ncells
      coefs_high = higher_order_recon(n).coef;
      px_high = cell.reconstruction_param.px;
      py_high = cell.reconstruction_param.py;
      
      coefs = regular_reconstruction(n).coef;
      px = regular_reconstruction_param.px;
      py = regular_reconstruction_param.py;
      
      xc = cell.xc(n,1:2);

      poly = @(dx,dy,coef,px,py) ((dx+xc(1)).^px.*(dy+xc(2)).^py)*coef;

      p1 = @(dx,dy) poly(dx,dy,coefs(:,1),px,py) - poly(dx,dy,coefs_high(:,1),px_high,py_high);
      p2 = @(dx,dy) poly(dx,dy,coefs(:,2),px,py) - poly(dx,dy,coefs_high(:,2),px_high,py_high);
      p3 = @(dx,dy) poly(dx,dy,coefs(:,3),px,py) - poly(dx,dy,coefs_high(:,3),px_high,py_high);
      p4 = @(dx,dy) poly(dx,dy,coefs(:,4),px,py) - poly(dx,dy,coefs_high(:,4),px_high,py_high);
%       coef_high = coefs(p_high,:);

%       xc = [0,0];

%       p1 = @(dx,dy) ((dx+xc(1)).^px.*(dy+xc(2)).^py).*Tcoef*coef_high(:,1);
%       p2 = @(dx,dy) ((dx+xc(1)).^px.*(dy+xc(2)).^py).*Tcoef*coef_high(:,2);
%       p3 = @(dx,dy) ((dx+xc(1)).^px.*(dy+xc(2)).^py).*Tcoef*coef_high(:,3);
%       p4 = @(dx,dy) ((dx+xc(1)).^px.*(dy+xc(2)).^py).*Tcoef*coef_high(:,4);
%       p1_coef = p1(0,0);
%       p2_coef = p2(0,0);
%       p3_coef = p3(0,0);
%       p4_coef = p4(0,0);
%       p1 = @(dx,dy) ((dx+xc(1)).^px.*(dy+xc(2)).^py).*Tcoef*coef_high(:,1)-p1_coef;
%       p2 = @(dx,dy) ((dx+xc(1)).^px.*(dy+xc(2)).^py).*Tcoef*coef_high(:,2)-p2_coef;
%       p3 = @(dx,dy) ((dx+xc(1)).^px.*(dy+xc(2)).^py).*Tcoef*coef_high(:,3)-p3_coef;
%       p4 = @(dx,dy) ((dx+xc(1)).^px.*(dy+xc(2)).^py).*Tcoef*coef_high(:,4)-p4_coef;

      te      = @(dx,dy) te_fun(dx,dy,coefs,px,py,poly);
      te_high = @(dx,dy) te_fun(dx,dy,coefs_high,px_high,py_high,poly);
      
%       Perr_cart1 = @(dx,dy) (abs(p1(dx,dy)) + abs(p2(dx,dy)) + abs(p3(dx,dy)) + abs(p4(dx,dy)) )/4;
      Perr_cart1 = @(dx,dy) sum( abs(te(dx,dy) - te_high(dx,dy) ) )/4; 
%       Perr_cart1 = @(dx,dy) sum( abs(te(dx,dy) - te_high(dx,dy) ) )/4; 
      Perr_cyln  = @(dr,th) Perr_cart1( dr*cos(th), dr*sin(th) );

      Perr_cart1_pagnutti = @(x,y) ( abs(-2/3*x^3 - 10*x^2*y - 8*x*y^2 + y^3/3) + abs(x^3 + 5*x^2*y + x*y^2 - 5/3*y^3) )/2;
      Perr_cyln_pagnutti  = @(dr,th) Perr_cart1_pagnutti( dr*cos(th), dr*sin(th) );

      th = linspace(0,2*pi,100);
      dth = th(2)-th(1);
      r=1;
      for ii = 1:length(th)
          err(ii) = Perr_cyln(1,th(ii));
          err_pagnutti(ii) = Perr_cyln_pagnutti(1,th(ii));
          err1(ii) = p1(r*cos(th(ii)), r*sin(th(ii)));
          err2(ii) = p2(r*cos(th(ii)), r*sin(th(ii)));
          err3(ii) = p3(r*cos(th(ii)), r*sin(th(ii)));
          err4(ii) = p4(r*cos(th(ii)), r*sin(th(ii)));
          errte(ii,:) = te( r*cos(th(ii)), r*sin(th(ii)) );
          errte_high(ii,:) = te_high( r*cos(th(ii)), r*sin(th(ii)) );
      end
%       err = err.^(2/3);
      % numerical integration using trapezoidal rule
      a0 = 1/pi*sum( dth/2*( err(1:end-1) + err(2:end) ) );
      a2 = 1/pi*sum( dth/2*( err(1:end-1).*cos(2*th(1:end-1)) + err(2:end).*cos(2*th(2:end)) ) );
      b2 = 1/pi*sum( dth/2*( err(1:end-1).*sin(2*th(1:end-1)) + err(2:end).*sin(2*th(2:end)) ) );
      
      fth_23 = @(th) (a0/2 + a2*cos(2*th) + b2*sin(2*th));
      figure(3)
      subplot(3,1,1)
      plot(th,err,'r-o',th,fth_23(th),'k-*',th,err1,'r-*',th,err2,'g-*',th,err3,'k-*',th,err4,'b-*')
      subplot(3,1,2)
      plot(th,errte(:,1),'r-v',th,errte(:,2),'g-v',th,errte(:,3),'k-v',th,errte(:,4),'b-v')
      subplot(3,1,3)
      plot(th,errte_high(:,1),'r-v',th,errte_high(:,2),'g-v',th,errte_high(:,3),'k-v',th,errte_high(:,4),'b-v')
%       axis([0,2*pi,0,5.5])
      
      M = [a0/2 + a2  b2
           b2  a0/2 - a2];
       
      dfun = @(x,y) M(1,1).*x.^2 + (M(1,2)+M(2,1)).*x.*y + M(2,2).*y.^2;
      fprintf('d(x,y) = ( %4.3f x^2 + %5.4f xy + %4.3f y^2)^(1/2)\n', M(1,1), (M(1,2)+M(2,1)), M(2,2))
      xtest = linspace(-1,1,25);
      [yt, xt] = meshgrid(xtest,xtest);
      figure(1)
      surf(xt,yt,dfun(xt,yt));
      
      figure(2)
      for ii = 1:numel(xt)
          err_plot(ii) = Perr_cart1(xt(ii),yt(ii));
      end
      err_plot = reshape(err_plot,size(xt));
      surf(xt,yt,err_plot);
%       Perr_cart1 = @(dx,dy) 

   end
   
function poly_test

  xc = [0,0];
  poly = @(dx,dy,coef,px,py) ((dx+xc(1)).^px.*(dy+xc(2)).^py)*coef;

  dx = 0.2;
  dy = 0.3;
  
  px = [0,4,0,4];
  py = [0,0,4,4];
  coefs = ones(size(px))';

  dpx = max(px-1,0);
  dpy = max(py-1,0);

  cpx = px;
  cpy = py;

  r   =  poly(dx,dy,coefs,px,py);
  drx =  poly(dx,dy,coefs.*cpx',dpx,py);
  dry =  poly(dx,dy,coefs.*cpy',px,dpy);
  
  dpx_t = [3,3];
  dpy_t = [0,4];
  cpx_t = [4,4];
  drx_t = sum(cpx_t.*dx.^dpx_t.*dy.^dpy_t);
  
  dpx_t = [0,4];
  dpy_t = [3,3];
  cpy_t = [4,4];
  dry_t =  sum(cpy_t.*dx.^dpx_t.*dy.^dpy_t);

  

   
   
function te = te_fun(x,y,coefs,px,py,poly)
  gamma = 1.4;
  
  dpx = max(px-1,0);
  dpy = max(py-1,0);

  cpx = px;
  cpy = py;

  r = @(dx,dy) poly(dx,dy,coefs(:,1),px,py);
  u = @(dx,dy) poly(dx,dy,coefs(:,2),px,py);
  v = @(dx,dy) poly(dx,dy,coefs(:,3),px,py);
  p = @(dx,dy) poly(dx,dy,coefs(:,4),px,py);
  e = @(dx,dy) p(dx,dy)./( r(dx,dy).*(gamma-1 ) ) + (u(dx,dy).^2 + v(dx,dy).^2)/2;
  
  drx = @(dx,dy) poly(dx,dy,coefs(:,1).*cpx',dpx,py);
  dux = @(dx,dy) poly(dx,dy,coefs(:,2).*cpx',dpx,py);
  dvx = @(dx,dy) poly(dx,dy,coefs(:,3).*cpx',dpx,py);
  dpx = @(dx,dy) poly(dx,dy,coefs(:,4).*cpx',dpx,py);
  dex = @(x,y) 1./(r(x,y).*(gamma-1)).*dpx(x,y) - p(x,y)./(r(x,y).^2*(gamma-1)).*drx(x,y) ...
       + u(x,y).*dux(x,y) + v(x,y).*dvx(x,y);
  
  dry = @(dx,dy) poly(dx,dy,coefs(:,1).*cpy',px,dpy);
  duy = @(dx,dy) poly(dx,dy,coefs(:,2).*cpy',px,dpy);
  dvy = @(dx,dy) poly(dx,dy,coefs(:,3).*cpy',px,dpy);
  dpy = @(dx,dy) poly(dx,dy,coefs(:,4).*cpy',px,dpy);
  dey = @(x,y) 1./(r(x,y).*(gamma-1)).*dpy(x,y) - p(x,y)./(r(x,y).^2*(gamma-1)).*dry(x,y) ...
       + u(x,y).*duy(x,y) + v(x,y).*dvy(x,y);
    
  te1_f = @(x,y) r(x,y).*dux(x,y) + drx(x,y).*u(x,y) ...
              + r(x,y).*dvx(x,y) + drx(x,y).*v(x,y);
          
  te2_f = @(x,y) dpx(x,y) ...
                + r(x,y).*u(x,y).*( dux(x,y) + dvy(x,y) )...
                + u(x,y).*( u(x,y).*drx(x,y) + v(x,y).*dry(x,y) )...
                + r(x,y).*( u(x,y).*dux(x,y) + v(x,y).*duy(x,y) );
            
  te3_f = @(x,y) dpy(x,y) ...
                + r(x,y).*v(x,y).*( dux(x,y) + dvy(x,y) )...
                + v(x,y).*( u(x,y).*drx(x,y) + v(x,y).*dry(x,y) )...
                + r(x,y).*( u(x,y).*dvx(x,y) + v(x,y).*dvy(x,y) );
 
  te4_f = @(x,y) r(x,y)*e(x,y)*( dux(x,y) + dvy(x,y) ) ...
             + r(x,y)*( u(x,y)*dex(x,y) + v(x,y)*dey(x,y) ) ...
             + e(x,y)*( u(x,y)*drx(x,y) + v(x,y)*dry(x,y) ) ...
             + p(x,y)*( dux(x,y) + dvy(x,y) ) ...
             + ( u(x,y)*dpx(x,y) + v(x,y)*dpy(x,y) );
         
  te(1) = te1_f(x,y);
  te(2) = te2_f(x,y);
  te(3) = te3_f(x,y);
  te(4) = te4_f(x,y);
  
