function [face] = compute_left_and_right_state(vertex, cell, face_in)
% Computes the left and right solution states for each face
% Inputs: 
%         vertex : list of nodes array of size [nnodes x 2]
%           cell : cell structure
%           face : face structure
%
% Outputs: 
%          face : face structure with face(n).ul and face(n).ur modified

% HARDWIRE : 1st order accurate only
% HARDWIRE : 1st order quadrature points for dirichlet bc 

face = face_in;


for n = 1 : length(face)

  r_cell = face(n).cell_plus;

  face(n).ur = cell.soln(r_cell, :);

  l_cell = face(n).cell_neg;
  % If not a boundary
  if (l_cell > 0)
    face(n).ul = cell.soln(l_cell, :);
  elseif (l_cell == -1)
    % DIRICHLET BC
    % Use simple geometric average of nodes to find the cell center
    xc(1) = sum( vertex(face(n).nodes(1:2),1) )/2;
    xc(2) = sum( vertex(face(n).nodes(1:2),2) )/2;

    rho = face(n).func.rho(xc);
    u = face(n).func.u(xc);
    v = face(n).func.v(xc);
    p = face(n).func.p(xc);
    face(n).ul = [rho, u, v, p];
  end

end
