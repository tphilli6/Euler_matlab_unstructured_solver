function [vertex, face, cell] = generate_mesh(imax, jmax, grid_type, neq)

perturb = 0.0;
seed = 0;

x0 = 0;
xL = 1;
y0 = 0;
yL = 1;

% Generate grid points  -------------------------------------------------------
x = zeros(imax, jmax);
y = zeros(imax, jmax);
% Cartesian
xi = linspace(x0,xL,imax);
eta = linspace(y0,yL,jmax);

for n = 1:imax
  x(n,:) = xi(n);
  y(n,:) = eta;
end

rng(seed);
dx = 1/(imax-1)*perturb;
dy = 1/(jmax-1)*perturb;
xdiff = 2*(rand( imax-2, jmax-2)-0.5)*dx;
ydiff = 2*(rand( imax-2, jmax-2)-0.5)*dy;
x(2:end-1,2:end-1) = x(2:end-1,2:end-1) + xdiff;
y(2:end-1,2:end-1) = y(2:end-1,2:end-1) + ydiff;

% Nozzle
%eta = linspace(0,2*pi,jmax);
%xi = linspace(x0,xL,imax);
%
%for n = 1:imax
%  x(n,:) = xi(n);
%  ymax = 1 + 0.5*cos(xi(n));
%  y(n,:) = linspace(-ymax,ymax,jmax);
%end



%% Build unstructured grid ----------------------------------------------------
[vertex, cell, face] = compute_grid_derived_data(x,y, grid_type);

%% Write grid to file %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fid = fopen('grid.vtk','w');
write_vtk(vertex, cell, fid);
fclose(fid);

fid = fopen('grid.vtk','a');
fprintf(fid,'CELL_DATA %8.0f\n',length(cell.volume));
write_vtk_data(cell.volume, 'volume', fid);
fclose(fid);


for n = 1:length(face)
  face(n).ul = zeros(1,neq);
  face(n).ur = zeros(1,neq);
end

