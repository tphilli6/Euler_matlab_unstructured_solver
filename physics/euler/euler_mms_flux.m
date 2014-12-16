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

rho = fun.rho(x);
u = fun.u(x);
v = fun.v(x);
p = fun.p(x);

flux = euler_flux([rho,u,v,p],normal);
