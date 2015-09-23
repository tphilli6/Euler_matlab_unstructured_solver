function [vertex, face, cell, vertex_cv, face_cv, cell_cv, cell_all, all_to_cv] = read_unstructured_mesh(grid_in, vertex_centered, scale)

[ti, vertex] = read_mesh_file(grid_in);

vertex = vertex*scale;

[face, cell] =  compute_face_data_from_triangulation(ti,vertex);
if vertex_centered
  [cell_cv, face_cv, vertex_cv, cell_all, all_to_cv] = compute_vertex_centered_nonsense(cell, face, vertex);
else
    vertex_cv = vertex; cell_cv = cell; face_cv = face; cell_all = cell; all_to_cv = [1:length(cell.volume)]';
end
cell_cv = find_neighbor_cells( cell_cv, face_cv );