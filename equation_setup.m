%% Setup the inputs required for solving the Euler equations

addpath('./multiquad');

% Set the equations to solve and parameters
equation = 'euler';
imax = 9;
jmax = 9;
iterations = 1000;
CFL = 0.15;
toler = 1e-10;
glb_dt = 0; %global time step, 1 = true, 0=false
mms_number=1; %1=supersonic, 2=subsonic
grid_type = 0; %quad = 0, triangles = 1, mixed = 2 (predecided mix)
restart = 0;

source_term_order=6;
exact_order=6;

kexact_order = 2;
flux_integral_order = 10;
kexact_type = 'kexact'; %'kexact', 'kexact_extended'
%kexact_type = 'kexact_extended';
% kexact : the sum of the exponents <= k
%          i.e. p = a + bx + cy
% kexact_extended : the maximum exponent for each dimension is k
%                 : i.e. p = a + bx + cy + dxy
fit_type = 'lsq'; %'lsq', 'kexact'
% lsq : nstencil = nunknowns + 1, conservation of the mean for only cell i
% lsq-true : nstencil = nunknowns + 1, conservation of the mean for only cell i
% kexact : nstencil = nunknowns, conservation of the mean for all cells

if (equation == 'euler')
  % Add to the path the subfolder with functions
  addpath('./physics/euler');
  gamma = 1.4; %HARDWIRED elsewhere
  neq  = 4;
  
  % Initial conditions, reference values, and used for catastrophic limiter
  rho_inf=1;
  u_inf = 800;
  v_inf = 800;
  p_inf = 100000;
  var_inf = [rho_inf, u_inf, v_inf, p_inf];
  
  var_lim = [1, 0, 0, 1];
  
  % Labels for output variables
  var = {'rho','u','v','p'};
  eq = {'mass','xmtm', 'ymtm', 'nrgy'};
  
  % Sets the discrete flux function
  flux = @(ul, ur, normal) flux_vanleer(ul, ur, normal);
  
  % Sets up the exact function
  exact_fun = setup_mms_crossterm(mms_number);
  
  % Sets up the exact solution
  analytic_soln = @(x) euler_analytic_solution(x,exact_fun); 

  % Analytic euler flux
  eq_flux = @(V, normal) euler_flux(V,normal);

  % Sets up the exact flux
  exact_flux = @(x, normal) euler_mms_flux(x, exact_fun, normal);
 
  % Sets up the time step
  local_time_step = @(vertex, cell, face, CFL, glb) euler_time_step(vertex, cell, face, CFL, glb);

end
