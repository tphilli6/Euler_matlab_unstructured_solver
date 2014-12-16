function [face] = compute_left_and_right_state(vertex, cell, face_in,...
                                               analytic_soln)
% Computes the left and right solution states for each face
% Inputs: 
%         vertex : list of nodes array of size [nnodes x 2]
%           cell : cell structure
%           face : face structure
%
% Outputs: 
%          face : face structure with face(n).ul and face(n).ur modified

xquad = cell.reconstruction_param.xquad;
wquad = cell.reconstruction_param.wquad;

face = face_in;

% Get the number of equations
test = analytic_soln([0,0]);

% debug print
%cell.reconstruction_param.px
%cell.reconstruction_param.py
%xquad
%wquad


for n = 1 : length(face)

  l_cell = face(n).cell_plus;

  coef = cell.reconstruction(l_cell).coef;
  px = cell.reconstruction_param.px;
  py = cell.reconstruction_param.py;
  
% HARDWIRE : linear face mapping
  x1 = vertex(face(n).nodes(1),1:2);
  x2 = vertex(face(n).nodes(2),1:2);

   for j = 1:length(wquad)
      % x mapping
      x(1) = (x2(1) - x1(1)) .* xquad(j) + x1(1);
      x(2) = (x2(2) - x1(2)) .* xquad(j) + x1(2);
 
      % loop over the number of equations
      for i = 1:size(coef,2)
        % evaluate the reconstruction at the quadrature points
        face(n).ul(j,i) = sum( coef(:,i)'.*x(1).^px.*x(2).^py );

      end


    end
  % First order
  %face(n).ul = cell.soln(l_cell, :);



  r_cell = face(n).cell_neg;
  % If not a boundary
  if (r_cell > 0)
    face(n).ur = cell.soln(r_cell, :);
    coef = cell.reconstruction(r_cell).coef;

   for j = 1:length(wquad)
      % x mapping
      x(1) = (x2(1) - x1(1)) .* xquad(j) + x1(1);
      x(2) = (x2(2) - x1(2)) .* xquad(j) + x1(2);
 
      % loop over the number of equations
      for i = 1:size(coef,2)
        % evaluate the reconstruction at the quadrature points
        face(n).ur(j,i) = sum( coef(:,i)'.*x(1).^px.*x(2).^py );
      end

    end




  elseif (r_cell == -1)
    % DIRICHLET BC
    soln = zeros(size(test));
    x1 = vertex(face(n).nodes(1),1:2);
    x2 = vertex(face(n).nodes(2),1:2);


    for i = 1:length(wquad)

% HARDWIRE : linear face mapping
      % x mapping
      x(1) = (x2(1) - x1(1)) .* xquad(i) + x1(1);
      x(2) = (x2(2) - x1(2)) .* xquad(i) + x1(2);
 
      soln(i,:) = analytic_soln(x);
    end

    face(n).ur = soln;
  end

end
