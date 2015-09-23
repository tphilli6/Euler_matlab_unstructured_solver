%% Compute the smooth grid TE from a given solution

clc
clear all




equation_setup
mms_number=1;
grid_type = -1; %quad = 0, triangles = 1, mixed = 2, from file = -1 (predecided mix)
grid_in = '0012-1397.mesh';
soln_file_in='naca0012-soln.vtk';

% Generates the simple mesh for the solver, also saves it as a vtk grid file. 
[vertex_grid, face_grid, cell_grid, vertex, face, cell, cell_tri, cell_tri_to_cv] = generate_mesh(imax, jmax, grid_type, neq, vertex_centered,refinement, grid_in);

setup_reconstruction;
write_sten=1;
build_kexact_stencil;


% Remove boundary te effects. Sets outside face equal to inside face;
for i = 1:length(face)
    if (face(i).cell_neg<0)
        face(i).cell_neg=-2; %Turns off mms boundary stuff and just sets ul=ur
    end
end

% vtk file for triangles 
% [soln, soln_points] = read_vtk(soln_file_in);
% 
% soln(:,1) = soln(:,1)*r_inf;
% soln(:,2) = soln(:,2)*c_inf;
% soln(:,3) = soln(:,3)*c_inf;
% soln(:,4) = soln(:,4)*(r_inf*c_inf^2);

% extract solution corresponding to vertex centered discretization
% for i = 1:cell.ncells
%     dsv = soln_points(:,1:2) - repmat(cell.xc(i,:), [size(soln_points,1),1]);
%     ds = sqrt( dsv(:,1).^2 + dsv(:,2).^2 );
%     [~,imin] = min(ds);
%     cell.soln(i,:) = soln(imin,:);
% end

cell.soln = compute_analytic_exact(cell_tri, vertex, exact_order, analytic_soln, cell_tri_to_cv);

% Write the MMS source term to file ----------------------------------------------------------------
cell.mms_source = analytic_flux(vertex, cell, face, exact_flux, neq,...
                                source_term_order);

reconstruction = reconstruct_solution(cell, fit_type, 1);
cell.reconstruction = reconstruction;
face_out = compute_left_and_right_state(vertex, cell, face,...
                                        analytic_soln);
face = face_out;

cell.te_exact = compute_residual( cell, face, flux );

te = cell.te_exact;
 for i = 1:4
     te(:,i) = te(:,i)./cell.volume';
 end

if vertex_centered;   data_write = vertex_to_cell_average(te,cell_grid); else data_write=te; end;                            
write_vtk_solution( vertex_grid, cell_grid, data_write, 'te', 'grid.vtk', 'a', eq ) 




if vertex_centered;   data_write = vertex_to_cell_average(cell.soln,cell_grid, cell); else data_write=cell.soln; end;
write_vtk_solution( vertex_grid, cell_grid, data_write, 'soln', ['grid.vtk'],'a',var )
    
% cell.mms_source = zeros(size(cell.soln));



te = te_kexact_estimate(cell, face, vertex, ..., ...
                                kexact_type, ...
                                fit_type, ...
                                kexact_order+1,...
                                flux_integral_order,...
                                analytic_soln,...
                                flux,...
                                eq_flux,...
                                vertex_centered,...
                                cell_grid);
                            

 for i = 1:4
     te(:,i) = te(:,i)./cell.volume';
 end
if vertex_centered;   data_write = vertex_to_cell_average(te,cell_grid); else data_write=te; end;                            
write_vtk_solution( vertex_grid, cell_grid, data_write, 'te-kexact', 'grid.vtk', 'a', eq )                            




te = te_smooth_grid(cell, face, vertex, ...
                        kexact_order, ...
                        kexact_type, ...
                        fit_type, ...
                        kexact_order+2,...
                        flux_integral_order,...
                        analytic_soln,...
                        flux,...
                        vertex_centered,...
                        cell_grid,...
                        face_grid,...
                        vertex_grid,...
                        cell_tri_to_cv,...
                        exact_flux);

for i = 1:4
    te(:,i) = te(:,i)./cell.volume';
end
if vertex_centered;   data_write = vertex_to_cell_average(te,cell_grid); else data_write=te; end;                            
write_vtk_solution( vertex_grid, cell_grid, data_write, 'te-smooth', 'grid.vtk', 'a', eq ) 