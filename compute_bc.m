function [face] = compute_bc(vertex, cell, face_in)
% Computes the boundary conditions at the faces
% Inputs :
%          vertex : list of nodes array of size [nnodes x 2]
%            cell : cell structure
%            face : face structure
%
% Outputs :
%            face : output face structure

% HARDWIRE : quadrature points, 1st order quadrature

face = face_in;

for n = 1:numel(face)
  if (face(n).cell_neg == -1) % Dirichlet bc

    % Use simple geometric average of nodes to find the cell center
    nc = cell.nodes(n,1);
    xc(1) = sum( vertex(cell.nodes(n,2:nc+1),1) )/nc;
    xc(2) = sum( vertex(cell.nodes(n,2:nc+1),2) )/nc;

% HACK 
   

    face(n).ul(1,1) = face(n).func.rho(xc);
    face(n).ul(1,2) = face(n).func.u(xc);
    face(n).ul(1,3) = face(n).func.v(xc);
    face(n).ul(1,4) = face(n).func.p(xc);

  else
    fprintf('Error! Boundary condition not recognized!\n');
  end
end
