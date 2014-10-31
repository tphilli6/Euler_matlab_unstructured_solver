%% Script to run an unstructured solver


imax = 9;
jmax = 9;
grid_type = 0; %quad = 0, triangles = 1, mixed = 2 (predecided mix)

% Generates the simple mesh for the solver, also saves it as a vtk grid file. 
% While not currently operational, can add a vtk input for a more general solution
[vertex, face, cell] = generate_mesh(imax, jmax, grid_type);


% Setup exact solution
exact_fun = setup_mms_crossterm;

% Compute exact solution
exact = analytic_solution(vertex, cell, exact_fun);
write_vtk_data(exact(:,1), 'rho')
write_vtk_data(exact(:,2), 'u')
write_vtk_data(exact(:,3), 'v')
write_vtk_data(exact(:,4), 'p')

mms_source = analytic_flux(vertex, face, exact_fun);
write_vtk_data(mms_source(:,1), 'rho')
write_vtk_data(mms_source(:,2), 'u')
write_vtk_data(mms_source(:,3), 'v')
write_vtk_data(mms_source(:,4), 'p')
