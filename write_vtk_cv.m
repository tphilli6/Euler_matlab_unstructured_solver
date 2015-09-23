function write_vtk_cv(vertex, cell, fid, stencil)
%% Writes vtk unstructured grid to file for use with paraview
% Inputs: 
%          vertex : 3 x n list of [x,y,z] node locations
%          cell.ncells     : number of cells (scalar)
%          cell.vtk_size   : cummulative size of cell.nodes (scalar)
%          cell.cell_type  : type of cell (0=quad, 1=triangle) (vector [ncells x 1])
%          cell.nodes      : node connectivity
%                            array [ncells x I] where I = [nnodes, n1, n2, n3,...]
%% Write unstructured grid to file %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin==3
    stencil = 1:cell.ncells;
end


% Write header
fprintf(fid,'# vtk DataFile Version 2.0\n');
fprintf(fid,'Unstructured Grid\n');
fprintf(fid,'ASCII\n\n');
% fprintf(fid,'DATASET UNSTRUCTURED_GRID\n');

fprintf(fid,'DATASET POLYDATA\n');

% Write out data points, doesn't matter what cell type
fprintf(fid,'POINTS %8.0f double\n',size(vertex,1));
for i = 1:size(vertex,1)
   fprintf(fid,'%23.15e %23.15e %23.15e\n',vertex(i,:));
end


% Write out connectivity

% fprintf(fid,'LINES%12.0f %12.0f\n',cell.ncells, cell.vtk_size+cell.ncells);
vtk_size = 2*numel(stencil) + sum(cell.nodes(stencil,1));
fprintf(fid,'LINES%12.0f %12.0f\n',numel(stencil), vtk_size);
for i = stencil
   fprintf(fid,'%8.0f', cell.nodes(i,1)+1 );
   for j = 1:cell.nodes(i,1)
     fprintf(fid, '%8.0f', cell.nodes(i,1+j)-1 );
   end
   fprintf(fid, '%8.0f', cell.nodes(i,2)-1 );
   fprintf(fid,'\n'); 
end


