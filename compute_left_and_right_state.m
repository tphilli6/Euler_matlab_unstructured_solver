function [face] = compute_left_and_right_state(vertex, cell, face_in,...
                                               korder, analytic_soln)
% Computes the left and right solution states for each face
% Inputs: 
%         vertex : list of nodes array of size [nnodes x 2]
%           cell : cell structure
%           face : face structure
%
% Outputs: 
%          face : face structure with face(n).ul and face(n).ur modified

% HARDWIRE : 1st order accurate only
% HARDWIRE : curtis-clenshaw quadrature for dirichlet bc 
if (korder >= 1)
  [xquad, wquad] = curtis_clenshaw( korder );
else
  xquad=[0.5];
  wquad=[1];
end

face = face_in;

% Get the number of equations
test = analytic_soln([0,0]);


for n = 1 : length(face)

  l_cell = face(n).cell_plus;

  face(n).ul = cell.soln(l_cell, :);

  r_cell = face(n).cell_neg;
  % If not a boundary
  if (r_cell > 0)
    face(n).ur = cell.soln(r_cell, :);
  elseif (r_cell == -1)
    % DIRICHLET BC
    soln = zeros(size(test));
    x1(1,:) = vertex(face(n).nodes(1:2),1);
    x2(2,:) = vertex(face(n).nodes(1:2),2);
    for i = 1:length(wquad)

% HARDWIRE : linear face mapping
      % x mapping
      x(1) = (x2(1) - x1(1)) .* xquad(i) + x1(1);
      x(2) = (x2(2) - x1(2)) .* xquad(i) + x1(2);
 
      soln = soln + analytic_soln(x)*wquad(i);
    end

    face(n).ur = soln;
  end

%    xc(1) = sum( vertex(face(n).nodes(1:2),1) )/2;
%    xc(2) = sum( vertex(face(n).nodes(1:2),2) )/2;
%
%    rho = face(n).func.rho(xc);
%    u = face(n).func.u(xc);
%    v = face(n).func.v(xc);
%    p = face(n).func.p(xc);
%    face(n).ur = [rho, u, v, p];
%    face(n).ul = [rho, u, v, p];

end
