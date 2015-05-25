function plot_cells(ti,x,h,style)
% Plot the triangulation ti

if nargin<3; h = 'off'; end
if nargin<4; style = 'k-o'; end
    

[face_nodes, cell_faces] = find_neighbors(ti);
for j = 1:size(face_nodes,1)
    if (j==1); eval(['hold ', h]); end
    plot( x(face_nodes(j,:),1), x(face_nodes(j,:),2) ,style,'LineWidth',2);
    hold on
end