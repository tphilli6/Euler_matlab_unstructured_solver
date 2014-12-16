function exact = analytic_solution(vertex, cell, order, ...
                                   analytic_soln)
% This function evaluates an analytic function over the cells
% Inputs: 
%         vertex : list of nodes (nnodes x 2)
%           cell : cell structure
%  analytic_soln : anonymous function which evaluates the exact solution

%HARDWIRE : curtis-clenshaw quadrature points

dim = 2;
method = 'CC'; % the only one available currently
quad_range_option_for_tri_transform = 2; % 1 for [-1,1] range

[xcc_quad, wcc] = sparse_grid(dim, dim*order, method);
for i = 1:length(wcc)
  [xcc_tri(i,:)] = transform_quad_to_triangle(xcc_quad(i,:), ...
                                  quad_range_option_for_tri_transform);
end

test = analytic_soln([0,0]);
neq = size(test,2);
exact = zeros(cell.ncells,neq);

for n = 1:cell.ncells

  % Use simple geometric average of nodes to find the cell center
  nc = cell.nodes(n,1);
  if (nc==3) % triangles
    xcc = xcc_tri;
  elseif (nc==4) % quadrilateral
    xcc = xcc_quad;
  end

  for i = 1:length(wcc)
    x(1) = cell.map(n).x(xcc(i,:));
    x(2) = cell.map(n).y(xcc(i,:));

    exact(n, :) = exact(n, :) + analytic_soln(x).*wcc(i);
  end

end


