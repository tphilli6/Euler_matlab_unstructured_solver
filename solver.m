%% Script to run an unstructured solver
clc
% clear all
clear resid update

equation_setup

% Generates the simple mesh for the solver, also saves it as a vtk grid file. 
% While not currently operational, can add a vtk input for a more general solution
[vertex_grid, face_grid, cell_grid, vertex, face, cell, cell_tri, cell_tri_to_cv] = generate_mesh(imax, jmax, grid_type, neq, vertex_centered,refinement, grid_in);
% load('vc_mesh3.mat')

% save('vc_mesh3.mat','vertex_grid','face_grid','cell_grid','vertex','face','cell','cell_tri','cell_tri_to_cv')
% load('vc_mesh4.mat')
% Computes the cell mapping required for quadratures
setup_reconstruction;
write_sten=1;
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
if vertex_centered;   data_write = vertex_to_cell_average(cell.mms_source,cell_grid); else data_write=cell.mms_source; end;                            
write_vtk_solution( vertex_grid, cell_grid, data_write, 'source', 'grid.vtk', 'a', eq )

% Setup and Compute the Exact Solution --------------------------------------------------------------
cell.exact = compute_analytic_exact(cell_tri, vertex, exact_order, analytic_soln, cell_tri_to_cv);

% HARDWIRE ----

% soln_file_in='naca0012-soln.vtk';
% [soln, soln_points] = read_vtk(soln_file_in);
% soln(:,1) = soln(:,1)*r_inf;
% soln(:,2) = soln(:,2)*c_inf;
% soln(:,3) = soln(:,3)*c_inf;
% soln(:,4) = soln(:,4)*(r_inf*c_inf^2);
% 
% % extract solution corresponding to vertex centered discretization
% for i = 1:cell.ncells
%     dsv = soln_points(:,1:2) - repmat(cell.xc(i,:), [size(soln_points,1),1]);
%     ds = sqrt( dsv(:,1).^2 + dsv(:,2).^2 );
%     [~,imin] = min(ds);
%     cell.soln(i,:) = soln(imin,:);
% end
% cell.exact = cell.soln;
% END HARDWIRE

if vertex_centered;   data_write = vertex_to_cell_average(cell.exact,cell_grid, cell); else data_write=cell.exact; end;
write_vtk_solution( vertex_grid, cell_grid, data_write, 'exact', 'grid.vtk','a',var )

%Compute Exact TE ----------------------------------------------------------------------------------
% Initializes the exact solution and computes the discrete residual to get the exact te
cell.soln = cell.exact;


reconstruction = reconstruct_solution(cell, fit_type, 1);
cell.reconstruction = reconstruction;
face_out = compute_left_and_right_state(vertex, cell, face,...
                                        analytic_soln);
face = face_out;

cell.te_exact = compute_residual( cell, face, flux );
% for j = 1:length(cell.te_exact(:,1))
%    f = cell.nbrs(j,1:cell.nnbr(j));
%    t = cell.te_exact(j,:)/2;
%    for i = 1
%        t = t + cell.te_exact(f(i),:)/(2*cell.nnbr(j));
%    end
%    te_new(j,:) = t;
% end
% cell.te_exact=te_new;
if vertex_centered;   data_write = vertex_to_cell_average(cell.te_exact,cell_grid, cell); else data_write=cell.te_exact; end;
write_vtk_solution( vertex_grid, cell_grid, data_write, 'exact-te', 'grid.vtk','a',eq )

% Initialize the solution variables and writes the initial conditions-------------------------------
for n = 1:neq
  cell.soln(:,n) = var_inf(n) * ones([size(cell.exact,1),1]);
end




% Setup initialization parameters
% clears parameters
if restart == 0
  cell.iteration = 1;
  cell.l2_to_normalize = ones(1,neq);

  cell.soln = cell.exact;

  if vertex_centered;   data_write = vertex_to_cell_average(cell.soln,cell_grid, cell); else data_write=cell.soln; end;
  write_vtk_solution( vertex_grid, cell_grid, data_write, 'soln', ['soln-',num2str(0),'.vtk'],'w',var )

elseif restart == 1
  load('output_solution.mat');

end



if dc_estimate == 1
    load('te_estimate.mat')
     if vertex_centered;   data_write = vertex_to_cell_average(te_smooth,cell_grid, cell); else data_write=te_smooth; end;
     write_vtk_solution( vertex_grid, cell_grid, data_write, 'te-for-dc', ['grid.vtk'],'w',eq )
     
    for ij = 1:size(te_smooth,2)
       te_smooth(:,ij) = te_smooth(:,ij)./cell.volume';
    end
    cell.mms_source = cell.mms_source + te_smooth;
    cell.primal = primal_soln;
    


end





% TIME LOOP ----------------------------------------------------------------------------------------

l2_to_normalize = ones(1,neq);
for iter = cell.iteration:cell.iteration + iterations - 1

  %Compute time step
  dt = local_time_step(vertex, cell, face, CFL, glb_dt);

  % Compute the reconstruction
%  cell.soln = exact;
  [reconstruction, cell.lhs_set] = reconstruct_solution(cell, fit_type, 1);
  cell.reconstruction = reconstruction;

  %Compute left and right face states
  % Not effecient but can only do this in matlab
  face_out = compute_left_and_right_state(vertex, cell, face,...
                                          analytic_soln);
  face = face_out;
 
  %Compute residual
  resid = compute_residual( cell, face, flux );
  
%   disp(vertex([1;2;11;10],1:2))
%   disp([face(1).ul(:,1), face(2).ul(:,1), face(4).ul(:,1),face(8).ul(:,1)])
%   disp([face(1).ur(:,1), face(2).ur(:,1), face(4).ur(:,1),face(8).ur(:,1)])
  
  % Check residual convergence
  [l2norms, converged, cell.l2_to_normalize] = check_convergence( resid, iter, cell.l2_to_normalize, toler );
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
 
%  % Write out new solution


    if vertex_centered;   data_write = vertex_to_cell_average(cell.soln,cell_grid, cell); else data_write=cell.soln; end;
    write_vtk_solution( vertex_grid, cell_grid, data_write, 'soln', ['soln-',num2str(iter),'.vtk'],'w',var )
    
    if vertex_centered;   data_write = vertex_to_cell_average(resid,cell_grid, cell); else data_write=resid; end;
    write_vtk_solution( vertex_grid, cell_grid, data_write, 'resid', ['soln-',num2str(iter),'.vtk'],'a',eq )
 

end

% Write out solution
if vertex_centered;   data_write = vertex_to_cell_average(cell.soln,cell_grid, cell); else data_write=cell.soln; end;
write_vtk_solution( vertex_grid, cell_grid, data_write, 'soln', ['soln-',num2str(iter),'.vtk'],'w',var )

if vertex_centered;   data_write = vertex_to_cell_average(resid,cell_grid, cell); else data_write=resid; end;
write_vtk_solution( vertex_grid, cell_grid, data_write, 'resid', ['soln-',num2str(iter),'.vtk'],'a',eq )

% Compute discretization error and write to file
if vertex_centered;   data_write = vertex_to_cell_average(cell.soln,cell_grid, cell); else data_write=cell.soln; end;
write_vtk_solution( vertex_grid, cell_grid, data_write, 'soln', 'grid.vtk','a',var )

de = cell.soln - cell.exact;

if vertex_centered;   data_write = vertex_to_cell_average(de,cell_grid, cell); else data_write=de; end;
write_vtk_solution( vertex_grid, cell_grid, data_write, 'exact-de', 'grid.vtk','a',var )

if dc_estimate == 1
    de_dc = cell.primal - cell.soln;
    
    if vertex_centered;   data_write = vertex_to_cell_average(de_dc,cell_grid, cell); else data_write=de_dc; end;
    write_vtk_solution( vertex_grid, cell_grid, data_write, 'dc-de', 'grid.vtk','a',var )
end


write_norms(de, 'exact-de')
cell.iteration = iter+1;
save('output_solution.mat','cell','vertex','face')

% soln_aprx_err = solution_approximation_error(cell, face, vertex, ...
%                         kexact_order, ...
%                         kexact_type, ...
%                         fit_type,...
%                         flux_integral_order,...
%                         analytic_soln,...
%                         flux);


if (dc_estimate==0)
%     compare_derivatives(cell, face, vertex, ...
%                         kexact_order, ...
%                         kexact_type, ...
%                         fit_type, ...
%                         kexact_order+1,...
%                         flux_integral_order,...
%                         analytic_soln,...
%                         flux,...
%                         vertex_centered,...
%                         cell_grid,...
%                         face_grid,...
%                         vertex_grid,...
%                         cell_tri_to_cv,...
%                         exact_flux,...
%                         exact_fun);
    
% te_smooth_grid(cell, face, vertex, ...
%                     kexact_order, ...
%                     kexact_type, ...
%                     fit_type, ...
%                     kexact_order+2,...
%                     flux_integral_order,...
%                     analytic_soln,...
%                     flux,...
%                     vertex_centered,...
%                     cell_grid,...
%                     face_grid,...
%                     vertex_grid,...
%                     cell_tri_to_cv,...
%                     exact_flux,...
%                     exact_fun);

% te = te_kexact_test(cell, face, vertex, ..., ...
%             kexact_type, ...
%             fit_type, ...
%             kexact_order+1,...
%             flux_integral_order,...
%             analytic_soln,...
%             flux,...
%             eq_flux,...
%             vertex_centered,...
%             cell_grid);
            
% te = te_kexact_estimate(cell, face, vertex, ..., ...
%                                 kexact_type, ...
%                                 fit_type, ...
%                                 kexact_order+1,...
%                                 flux_integral_order,...
%                                 analytic_soln,...
%                                 flux,...
%                                 eq_flux,...
%                                 vertex_centered,...
%                                 cell_grid);
                            
% te = te_kexact_estimate_lsq(cell, face, vertex, ..., ...
%                                 kexact_type, ...
%                                 fit_type, ...
%                                 kexact_order+1,...
%                                 flux_integral_order,...
%                                 analytic_soln,...
%                                 flux,...
%                                 eq_flux,...
%                                 vertex_centered,...
%                                 cell_grid);
%                             
% te = te_kexact_estimate_lsq_point(cell, face, vertex, ..., ...
%                                 kexact_type, ...
%                                 fit_type, ...
%                                 kexact_order+1,...
%                                 flux_integral_order,...
%                                 analytic_soln,...
%                                 flux,...
%                                 eq_flux,...
%                                 vertex_centered,...
%                                 cell_grid); 
end





% Info on cell data type
% Cell data type is struct of array so indexing is cell.variable(i)
% cell.
%      ncells = number of cells
%      volume = cell volume [ncells]
%      cell_type = 0 for quad, 1, for triangle
%      nodes = node indices [ nvertex, verte1, vertex2, ...]
%      soln =  cell solution [ncells x soln vars]
%      exact = exact cell solution [ncells x soln vars]
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
