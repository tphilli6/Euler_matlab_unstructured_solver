function [ti, vertex] = generate_kt_airfoil(mu, te_angle)
%% Generate KT airfoil

% mu = complex(-0.0475,0);
% xc = real(mu);
% te_angle = -(0.594689181*(1/2*0.298222773/sqrt(1) - 0.127125232 - 2*0.357907906*1 + 3*0.291984971*1*1 - 4*0.105174606*1*1*1));
n = 2 - te_angle/pi;
r = sqrt( (1-real(mu)).^2 + imag(mu).^2 );

imax = 80;
ds = 2*pi*r/imax;
dth = 2*pi/imax;
gf = 1.3;
nmax = 1200;
rmax = floor(nmax/imax);

R(1) = r;
R(2) = R(1)+ds;
for i = 3:rmax
   R(i) = (R(i-1)-R(i-2))*gf + R(i-1);
end

Theta = 0:dth:2*pi-dth;


cnt = 1;
for rr = 1:length(R)
    for tt = 1:length(Theta)
      x(cnt) = R(rr)*cos(Theta(tt)) + real(mu);
      y(cnt) = R(rr)*sin(Theta(tt)) + imag(mu);
      cnt = cnt + 1;
    end
    Theta = Theta+dth/2;
end

xi_new = [x',y']; 
[ti, xi] = delaunay_triangulation( xi_new(:,1),xi_new(:,2) );

zeta = complex(xi(:,1),xi(:,2));
% n=2;
z = n*(  (1+1./zeta).^n + (1-1./zeta).^n )./(  (1+1./zeta).^n - (1-1./zeta).^n );
zscale = z(1:imax);
% z = (z-min(real(zscale)))./(max(real(zscale))-min(real(zscale)));

vertex = [real(z),imag(z)];
ti = ti(imax-1:end,:);
% plot_cells(ti,vertex)
% axis equal

% rho = karman_trefftz_airfoil(vertex, 1, 1, 100000, 400, 0, mu, te_angle);
% 
% u = karman_trefftz_airfoil(vertex, 2, 1, 100000, 400, 0, mu, te_angle);
% 
% v = karman_trefftz_airfoil(vertex, 3, 1, 100000, 400, 0, mu, te_angle);
% 
% p = karman_trefftz_airfoil(vertex, 4, 1, 100000, 400, 0, mu, te_angle);
% 
% plot3( vertex(:,1), vertex(:,2), u' )