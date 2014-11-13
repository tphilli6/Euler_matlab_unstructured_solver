%% Script to run an unstructured solver
%clc
clear all

imax = 17;
jmax = 17;
iterations = 100;
CFL = 0.2;
toler = 1e-10;
glb_dt = 0; %global time step, 1 = true, 0=false
mms_number=1; %1=supersonic, 2=subsonic
grid_type = 0; %quad = 0, triangles = 1, mixed = 2 (predecided mix)

gamma = 1.4; %HARDWIRED elsewhere
neq  = 4;

% Initial conditions, reference values, and used for catastrophic limiter
rho_inf=1;
u_inf = 800;
v_inf = 800;
p_inf = 100000;

var = {'rho','u','v','p'};
eq = {'mass','xmtm', 'ymtm', 'nrgy'};

flux = @(ul, ur, normal) flux_vanleer(ul, ur, normal);

% Generates the simple mesh for the solver, also saves it as a vtk grid file. 
% While not currently operational, can add a vtk input for a more general solution
[vertex, face, cell] = generate_mesh(imax, jmax, grid_type, neq);

% Setup exact solution --------------------------
exact_fun = setup_mms_crossterm(mms_number);

% Compute exact solution
exact = analytic_solution(vertex, cell, exact_fun);
write_vtk_solution( vertex, cell, exact, 'exact', 'grid.vtk','a',var )

% Initialize the solution variables -------------
cell.soln = exact;

cell.soln(:,1) = rho_inf * ones(size(exact(:,1)));
cell.soln(:,2) = u_inf * ones(size(exact(:,1)));
cell.soln(:,3) = v_inf * ones(size(exact(:,1)));
cell.soln(:,4) = p_inf * ones(size(exact(:,1)));

% Write the MMS source term to file -------------
cell.mms_source = analytic_flux(vertex, cell, face, exact_fun);
write_vtk_solution( vertex, cell, cell.mms_source, 'source', 'grid.vtk','a',eq )

% Setup dirichlet bc
% assigns the exact solution to apply a dirichlet bc computed later
for n = 1 : length(face)
  face(n).func = exact_fun;
end


write_vtk_solution( vertex, cell, cell.soln, 'soln', ['soln-',num2str(0),'.vtk'],'w',var )

% TIME LOOP -----------------------------------------------------------

l2_to_normalize = ones(1,neq);
for iter = 1:iterations

  %Compute time step
  dt = local_time_step(vertex, cell, face, CFL, glb_dt);
 
  %Compute left and right face states
  % Not effecient but can only do this in matlab
 
  face_out = compute_left_and_right_state(vertex, cell, face);
  face = face_out;
 
  %Compute residual
  resid = compute_residual( cell, face, flux );

  % Check residual convergence
  [l2norms, converged, l2_to_normalize] = check_convergence( resid, iter, l2_to_normalize, toler );
  fprintf('%6.0f %12.6e %12.6e %12.6e %12.6e \n', iter, l2norms);
  if (converged)
    break
  end

%Compute updated solution
U = V_to_U(cell.soln, gamma);
constant = dt'./cell.volume';

for i = 1:neq
  update(:,i) = constant.*resid(:,i);
end

Unew = U - update;
soln = U_to_V(Unew, gamma);

cell.soln = soln;
cell.soln(:,1) = max(soln(:,1), 0.01*rho_inf);
cell.soln(:,4) = max(soln(:,4), 0.01*p_inf);

write_vtk_solution( vertex, cell, cell.soln, 'soln', ['soln-',num2str(iter),'.vtk'],'w',var )
write_vtk_solution( vertex, cell, abs(resid), 'resid', ['soln-',num2str(iter),'.vtk'],'a',eq )
 

end

write_vtk_solution( vertex, cell, cell.soln, 'soln', ['soln-',num2str(iter),'.vtk'],'w',var )
write_vtk_solution( vertex, cell, resid, 'resid', ['soln-',num2str(iter),'.vtk'],'a',eq )











% Info on cell data type
% Cell data type is struct of array so indexing is cell.variable(i)
% cell.
%      ncells = number of cells
%      volume = cell volume [ncells]
%      cell_type = 0 for quad, 1, for triangle
%      nodes = node indices [ nvertex, verte1, vertex2, ...]
%      soln = cell solution [ncells x soln vars]
%      nfaces = number of faces per cell [ncells]
%      faces = list of face indices for each cell [ ncells, nfaces]

% Info on face data type
% Face data type is array of struct so indexing is face(i).variable
% face.
%      cell_neg  = neighboring cell index
%      cell_plus = cell index
%      normal    = face normal
%      area      = face area
%      nodes     = node list making up face [2]
