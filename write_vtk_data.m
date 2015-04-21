function write_vtk_data(data, name, fid)
%% Writes vtk unstructured grid to file for use with paraview
% Inputs: 
%          vertex : 3 x n list of [x,y,z] node locations
%          cell.ncells     : number of cells (scalar)
%          cell.vtk_size   : cummulative size of cell.nodes (scalar)
%          cell.cell_type  : type of cell (0=quad, 1=triangle) (vector [ncells x 1])
%          cell.nodes      : node connectivity
%                            array [ncells x I] where I = [nnodes, n1, n2, n3,...]
%% Write unstructured grid to file %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% Write out data points, doesn't matter what cell type
% fprintf(fid,'POINT_DATA %8.0f\n',length(data));
fprintf(fid,['SCALARS ',name,' double 1\n'],length(data));
fprintf(fid,'LOOKUP_TABLE default\n');

for i = 1:length(data)
   fprintf(fid,'%23.15e\n',data(i));
end

