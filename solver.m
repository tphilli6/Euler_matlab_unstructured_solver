%% Script to run an unstructured solver
%clc
clear all

equation_setup

% Generates the simple mesh for the solver, also saves it as a vtk grid file. 
% While not currently operational, can add a vtk input for a more general solution
[vertex, face, cell] = generate_mesh(imax, jmax, grid_type, neq);
% Computes the cell mapping required for quadratures
setup_mapping;
setup_reconstruction;
build_kexact_stencil;

% Setup dirichlet bc -------------------------------------------------------------------------------
% assigns the exact solution to apply a dirichlet bc computed later
% Sets all faces to have an exact function, only the ones at boundaries are used
for n = 1 : length(face)
  face(n).func = exact_fun;
end

% Write the MMS source term to file ----------------------------------------------------------------
cell.mms_source = analytic_flux(vertex, cell, face, exact_flux, neq,...
                                source_term_order);
write_vtk_solution( vertex, cell, cell.mms_source, 'source', 'grid.vtk','a',eq )

% Setup and Compute the Exact Solution --------------------------------------------------------------
exact = analytic_solution(vertex, cell, exact_order, analytic_soln);
write_vtk_solution( vertex, cell, exact, 'exact', 'grid.vtk','a',var )

%Compute Exact TE ----------------------------------------------------------------------------------
% Initializes the exact solution and computes the discrete residual to get the exact te
cell.soln = exact;
kexact_order
face_out = compute_left_and_right_state(vertex, cell, face,...
                                        kexact_order, analytic_soln);
face = face_out;

te = compute_residual( cell, face, flux );
write_vtk_solution( vertex, cell, te, 'exact-te', 'grid.vtk','a',eq )

% Initialize the solution variables and writes the initial conditions-------------------------------
for n = 1:neq
  cell.soln(:,n) = var_inf(n) * ones([size(exact,1),1]);
end
write_vtk_solution( vertex, cell, cell.soln, 'soln', ['soln-',num2str(0),'.vtk'],'w',var )





% TIME LOOP ----------------------------------------------------------------------------------------

l2_to_normalize = ones(1,neq);
for iter = 1:iterations

  %Compute time step
  dt = local_time_step(vertex, cell, face, CFL, glb_dt);
 
  %Compute left and right face states
  % Not effecient but can only do this in matlab
  face_out = compute_left_and_right_state(vertex, cell, face,...
                                          kexact_order, analytic_soln);
  face = face_out;
 
  %Compute residual
  resid = compute_residual( cell, face, flux );

  % Check residual convergence
  [l2norms, converged, l2_to_normalize] = check_convergence( resid, iter, l2_to_normalize, toler );
  fprintf('%6.0f %12.6e %12.6e %12.6e %12.6e \n', iter, l2norms);
  if (converged)
    fprintf('Solution has converged!!\n');
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
  for n = 1:neq
    if (var_lim(n) == 1)
      cell.soln(:,n) = max( soln(:,n), 0.01*var_inf(n) );
    end
  end
 
  % Write out new solution
  write_vtk_solution( vertex, cell, cell.soln, 'soln', ['soln-',num2str(iter),'.vtk'],'w',var )
  write_vtk_solution( vertex, cell, abs(resid), 'resid', ['soln-',num2str(iter),'.vtk'],'a',eq )
 

end

% Write out solution
write_vtk_solution( vertex, cell, cell.soln, 'soln', ['soln-',num2str(iter),'.vtk'],'w',var )
write_vtk_solution( vertex, cell, resid, 'resid', ['soln-',num2str(iter),'.vtk'],'a',eq )

% Compute discretization error and write to file
de = cell.soln - exact;
write_vtk_solution( vertex, cell, de, 'exact-de', 'grid.vtk','a',eq )

write_norms(de, 'exact-de')



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
%      nnbr = number of neighboring cells
%      nbrs  = list of cell indicies of neighboring cells
% Info on face data type
% Face data type is array of struct so indexing is face(i).variable
% face.
%      cell_neg  = neighboring cell index
%      cell_plus = cell index
%      normal    = face normal
%      area      = face area
%      nodes     = node list making up face [2]
