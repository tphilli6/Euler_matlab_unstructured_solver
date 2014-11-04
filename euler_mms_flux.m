function flux = euler_mms_flux(x, fun, normal)
% Computes the Euler analytic flux for a given function
% Inputs:
%         x : array with x and y location to evaluate functions
%       fun : Density,    fun.rho(x,y)  
%             X-velocity, fun.u(x,y)
%             Y-Velocity, fun.v(x,y)
%             Pressure,   fun.p(x,y)
%    normal : face normal [2x1]
%
% Outputs:
%      flux : analytic flux
% HARDWIRE: gamma = 1.4
gamma = 1.4;

rho = fun.rho(x);
u = fun.u(x);
v = fun.v(x);
p = fun.p(x);

v_norm = u*normal(1) + v*normal(2);
ht = gamma/(gamma-1)*p/rho + 0.5*(u^2 + v^2);

flux(1) = rho*v_norm;
flux(2) = rho*u*v_norm + normal(1)*p;
flux(3) = rho*v*v_norm + normal(2)*p;
flux(4) = rho*ht*v_norm;

