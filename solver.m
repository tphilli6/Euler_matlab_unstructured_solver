%% Script to run an unstructured solver
clc
clear all

imax = 17;
jmax = 17;
neq  = 4;
iterations = 100;
grid_type = 0; %quad = 0, triangles = 1, mixed = 2 (predecided mix)
CFL = 0.6;
toler = 1e-10;
gamma = 1.4; %HARDCODED elsewhere
glb_dt = 0; %global time step, 1 = true, 0=false

flux = @(ul, ur, normal) flux_vanleer(ul, ur, normal);

% Generates the simple mesh for the solver, also saves it as a vtk grid file. 
% While not currently operational, can add a vtk input for a more general solution

[vertex, face, cell] = generate_mesh(imax, jmax, grid_type, neq);

% Setup exact solution --------------------------
exact_fun = setup_mms_crossterm;

% Compute exact solution
exact = analytic_solution(vertex, cell, exact_fun);
write_vtk_data(exact(:,1), 'rho')
write_vtk_data(exact(:,2), 'u')
write_vtk_data(exact(:,3), 'v')
write_vtk_data(exact(:,4), 'p')

% Initialize the solution variables -------------
cell.soln = exact;

% Write the MMS source term to file -------------
cell.mms_source = analytic_flux(vertex, cell, face, exact_fun);
write_vtk_data(cell.mms_source(:,1), 'mass')
write_vtk_data(cell.mms_source(:,2), 'xmtm')
write_vtk_data(cell.mms_source(:,3), 'ymtm')
write_vtk_data(cell.mms_source(:,4), 'nrgy')

% Setup dirichlet bc
% assigns the exact solution to apply a dirichlet bc computed later
for n = 1 : length(face)
  face(n).func = exact_fun;
end


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
cell.soln = U_to_V(Unew, gamma);


end


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
