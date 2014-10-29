%% Script to run an unstructured solver


imax = 17;
jmax = 17;
grid_type = 0; %quad = 0, triangles = 1, mixed = 2 (predecided mix)

% Generates the simple mesh for the solver, also saves it as a vtk grid file. 
% While not currently operational, can add a vtk input for a more general solution
[vertex, face, cell] = generate_mesh(imax, jmax, grid_type);
