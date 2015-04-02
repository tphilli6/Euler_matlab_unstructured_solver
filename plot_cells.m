function plot_cells(ti,x)
% Plot the triangulation ti

[face_nodes, cell_faces] = find_neighbors(ti);
for j = 1:size(face_nodes,1)
    if (j==1); hold off; end
    plot( x(face_nodes(j,:),1), x(face_nodes(j,:),2) ,'k-o','LineWidth',2);
    hold on
end