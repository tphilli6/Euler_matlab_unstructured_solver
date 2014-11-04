function exact = analytic_solution(vertex, cell, func)
% This function evaluates an analytic function over the cells
% Inputs: 
%         vertex : list of nodes (nnodes x 2)
%           cell : cell structure
%           func : func structure for primitive variables
%                  (func.rho, func.u, func.v, func.p)

%HARDWIRE : quadrature points, 1st order quadrature


for n = 1:cell.ncells

  % Use simple geometric average of nodes to find the cell center
  nc = cell.nodes(n,1);
  xc(1) = sum( vertex(cell.nodes(n,2:nc+1),1) )/nc;
  xc(2) = sum( vertex(cell.nodes(n,2:nc+1),2) )/nc;


  exact(n, 1) = func.rho(xc);
  exact(n, 2) = func.u(xc);
  exact(n, 3) = func.v(xc);
  exact(n, 4) = func.p(xc);

end


