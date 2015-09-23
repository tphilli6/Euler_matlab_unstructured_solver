%% Derivative Tests

clc
clear all


equation_setup
mms_number=1;
grid_type = -1; %quad = 0, triangles = 1, mixed = 2, from file = -1 (predecided mix)
grid_in = 'equilateral.mesh';


% Generates the simple mesh for the solver, also saves it as a vtk grid file. 
[vertex_grid, face_grid, cell_grid, vertex, face, cell, cell_tri, cell_tri_to_cv] = generate_mesh(imax, jmax, grid_type, neq, vertex_centered,refinement, grid_in);

setup_reconstruction;
write_sten=1;
build_kexact_stencil;