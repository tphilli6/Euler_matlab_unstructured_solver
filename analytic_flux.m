function source = analytic_flux( vertex, cell, face, flux_func, neq,...
                                 order, ncell)
% This function loops over the faces and computes the analytic source
% term.
% Inputs:
%         vertex : list of nodes (nnodes x 2)
%         cell   : cell structure for sizing
%         face   : face structure with relevant face data
%      flux_func : the flux function used to evaluate whatever
%         neq    : number of equations
%         order  : quadrature order for source terms

if (nargin == 6)
    face_loop = 1:numel(face);
elseif (nargin == 7)
    face_loop = cell.faces(ncell,:);
end

% HARDWIRE : curtis-clenshaw quadrature points
[xquad, wquad] = curtis_clenshaw( order );

source = zeros([numel(cell.volume),neq]);
    %temp
cnt = zeros( numel(cell.volume) );
    %temp
flux_check = zeros([numel(cell.volume),neq]);

for n = face_loop

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

%  if (face(n).cell_plus==1)
%  fprintf('anal flux for cell 1, face %8.4f\n',n)
%    fprintf('%8.4f %8.4f %23.15e %23.15e %23.15e %23.15e %8.4f \n',x(1), x(2), flux_func(x, face(n).normal),wquad(i))
%  end
    flux = flux + flux_func(x, face(n).normal)*wquad(i);

  end

  flux = flux * f.area;

  source(f.cell_plus,:) = source(f.cell_plus,:) + flux;

    %temp
  cnt(f.cell_plus) = cnt(f.cell_plus) + 1;
    %temp
  flux_check(f.cell_plus,cnt(f.cell_plus)) = flux(4);
%  disp([f.cell_plus,f.cell_neg, flux])
  if (f.cell_neg>0) 
    source(f.cell_neg,:) = source(f.cell_neg,:) - flux;
    %temp
    cnt(f.cell_neg) = cnt(f.cell_neg) + 1;
    %temp
    flux_check(f.cell_neg,cnt(f.cell_neg)) = -flux(4);
  end


end


for n = 1:size(source,2)
  source_new(:,n) = source(:,n)./cell.volume'; 
end
source = source_new;
