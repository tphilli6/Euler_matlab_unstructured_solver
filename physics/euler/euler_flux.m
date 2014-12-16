function flux = euler_flux(V, normal)
% Computes the Euler analytic flux for a given function
% Inputs:
%       rho : Density
%       u     X-velocity
%       v     Y-Velocity
%       p     Pressure
%    normal : face normal [2x1]
%
% Outputs:
%      flux : analytic flux
% HARDWIRE: gamma = 1.4
rho = V(1);
u = V(2);
v = V(3);
p = V(4);

gamma = 1.4;


v_norm = u*normal(1) + v*normal(2);
ht = gamma/(gamma-1)*p/rho + 0.5*(u^2 + v^2);

flux(1) = rho*v_norm;
flux(2) = rho*u*v_norm + normal(1)*p;
flux(3) = rho*v*v_norm + normal(2)*p;
flux(4) = rho*ht*v_norm;
