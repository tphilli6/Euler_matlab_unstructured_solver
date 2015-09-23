function [ti, vertex] = read_mesh_file(file)
%read the mesh file given by 'file'
% format:
%
% ? ? ? nnodes
% vertex 1
% vertex 2
% ...
% vertex nnodes
% file = '0012-1397.mesh';

fid = fopen(file,'r');
file_data_str = fgetl(fid);
file_data = sscanf(file_data_str, '%i %i %i %i')';

nvertex = file_data(4);
ncell = file_data(1);

for i = 1:file_data(4)
    node_str = fgetl(fid);
    vertex(i,:) = [sscanf(node_str, '%f %f')',0];
end

for i = 1:file_data(2)
    face_str = fgetl(fid);
    face(i,:) = sscanf(face_str, '%i %i %i %i')';
end
% face = sortrows(face,1);

face_nodes = face(:,3:4)+1; % nodes that make up a face
face_cells = face(:,1:2); % cell neighbor list for each face

% shift to matlab indexing
I=find(face_cells>=0);
face_cells(I) = face_cells(I) + 1;


for i = 1:ncell
    I = [find(face_cells(:,1)==i); find(face_cells(:,2)==i)];
    cell_faces(i,:) = sort(I)';
    ti(i,:) = unique(face_nodes(I,:))';
% ti/cell_nodes
end


% [face_nodes, cell_faces, face_cells] = find_neighbors(ti);


