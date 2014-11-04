function source = analytic_flux( vertex, cell, face, func)
% This function loops over the faces and computes the analytic source
% term.
% Inputs:
%         vertex : list of nodes (nnodes x 2)
%         cell   : cell structure for sizing
%         face   : face structure with relevant face data
%         func   : func structure for primitive variables 
%                  (func.rho, func.u, func.v, func.p)

% HARDWIRE : quadrature points, 1st order quadrature
xquad = 0.5;
wquad = 1.0;

source = zeros(size(cell.soln));
cnt = zeros(size(cell.soln(:,1)));

for n = 1:length(face)

  % store off face to prevent redudant memory access
  f = face(n);

  x1 = vertex( f.nodes(1), : );
  x2 = vertex( f.nodes(2), : );

  % Evaluate function at quadrature points
  flux = 0;
  for i = 1:length(wquad)

    % HARDWIRE : linear face mapping
    % x mapping
    x(1) = (x2(1) - x1(1)) .* xquad(i) + x1(1);
    x(2) = (x2(2) - x1(2)) .* xquad(i) + x1(2);

    flux = flux + euler_mms_flux(x, func, face(n).normal)*wquad(i);

  end

  flux = flux * f.area;

  source(f.cell_plus,:) = source(f.cell_plus,:) + flux;

  cnt(f.cell_plus) = cnt(f.cell_plus) + 1;
%  disp([f.cell_plus,f.cell_neg, flux])
  if (f.cell_neg>0) 
    source(f.cell_neg,:) = source(f.cell_neg,:) - flux;
    cnt(f.cell_neg) = cnt(f.cell_neg) + 1;
  end


end
