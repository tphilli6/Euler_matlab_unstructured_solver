function [ti,xvec] = delaunay_triangulation(x,y)

rx = max(max(x))-min(min(x));
ry = max(max(y))-min(min(y));


% Generate bounding box
xbox = [min(min(x))-rx/2, min(min(y))-ry/2
        max(max(x))+rx/2, min(min(y))-ry/2
        max(max(x))+rx/2, max(max(y))+ry/2
        min(min(x))-rx/2, max(max(y))+ry/2];

xvec = [reshape(x,[numel(x),1]), reshape(y,[numel(y),1])];
    
% and concatinate to current list of nodes
nodes = [xbox
         xvec];

% Current cells
ti = [1, 2, 4
      4, 2, 3];

for i = 1:size(nodes,1)-4
    [ti, ~] = delaunay_add_point(ti,xvec(i,:),nodes(1:i+4-1,:));
end

ti = remove_corners(ti);