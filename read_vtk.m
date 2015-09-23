function [soln, vertex] = read_vtk(file)
%% Writes vtk unstructured grid to file for use with paraview
% Inputs: 
%          vertex : 3 x n list of [x,y,z] node locations
%          cell.ncells     : number of cells (scalar)
%          cell.vtk_size   : cummulative size of cell.nodes (scalar)
%          cell.cell_type  : type of cell (0=quad, 1=triangle) (vector [ncells x 1])
%          cell.nodes      : node connectivity
%                            array [ncells x I] where I = [nnodes, n1, n2, n3,...]
%% Write unstructured grid to file %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% file='naca0012-soln.vtk';

fid = fopen(file,'r');

% READ header
% # vtk DataFile Version 2.0
fgetl(fid);
% Unstructured Grid
fgetl(fid);
% ASCII
fgetl(fid);
% DATASET UNSTRUCTURED_GRID
fgetl(fid);


% Write out data points, doesn't matter what cell type
% fprintf(fid,'POINTS %8.0f double\n',size(vertex,1));
points_data = fgetl(fid);
points = sscanf(points_data(8:end-6), '%i');
point_cnt = 0;
vertex = [];
while point_cnt < points
    node_str = fgetl(fid);
    v_temp = sscanf(node_str, '%f')';
    nrow = length(v_temp)/3;
    vertex = [vertex; reshape(v_temp,[3,nrow])' ];
    point_cnt = point_cnt + nrow;
%    fprintf(fid,'%23.15e %23.15e %23.15e\n',vertex(i,:));
end

% Read in connectivity
% fprintf(fid,'CELLS %12.0f %12.0f\n',ncells, vtk_size);
cell_data_str = fgetl(fid);
cell_data = sscanf(cell_data_str(6:end), '%i');
ncells = cell_data(1);
for i = 1:ncells
    nodes = sscanf( fgetl(fid), '%i')';
    nodes(2:end) = nodes(2:end)+1;
    ti(i,:) = nodes(2:end);
    cell.nodes(i,:) = nodes;
end

fgetl(fid); % read space
fgetl(fid); % 'CELL_TYPES #'
for i = 1:ncells
    fgetl(fid); % read cell type (not needed)
end


fgetl(fid); % read space
fgetl(fid); % read 'POINT_DATA #' #=points from earlier
field_str = fgetl(fid); % read 'FIELD ' #=points from earlier
nvar = sscanf( field_str(16:end),'%i');
soln = [];
for i = 1:nvar
    var = [];
    fgetl(fid); % read '$\var$ 1 # double'

    cell_cnt = 0;
    while cell_cnt < points
        data = sscanf( fgetl(fid), '%f' );
        var = [var; data];
        cell_cnt = cell_cnt + length(data);
    end

    soln = [soln, var];

end
% 
% 
% % Write out cell type
% fprintf(fid,'CELL_TYPES %12.0f\n',ncells);
% for i = 1:ncells
%    if (cell.cell_type(stencil(i)) == 0)
%      fprintf(fid,'9\n'); % For quad vtk cell type id
%    elseif (cell.cell_type(stencil(i)) == 1)
%      fprintf(fid,'5\n'); % For triangles vtk cell type id
%    end
% end
% 
