%setup reconstruction for each cell
% Requires the definition of 
% kexact_order : order of reconstruction
% kexact_type  : 'kexact', 'kexact-extended'

% This is a script that will setup the row in the left hand side which is related to the discrete integral a polynomial of order k 
% cell.recon(n).Arow

% For quads, the quadrature used is curtis clenshaw. For each control point a function evaluation can be computed for a given polynomial
% p(x) = a+bx+cy+dxy which in matrix form is 
% Ai*coef = ui
%and
% Ai = [1, map(i).x(xcc), map(i).y(xcc), map(i).x(xcc)*map(i).y(xcc),...]

%some basic definitions
dim = 2;
method = 'CC'; % the only one available currently
quad_range_option_for_tri_transform = 2; % 1 for [-1,1] range
                                         % 2 for [0,1] range

%Setup polynomial matrix based on kexact type
% if the original kexact method, then the sum of the exponents for each term is less than or equal to the kexact order
cell.nunknowns=1;

% if the kexact order is zero (order polynomial) then kexact is the only possible setup so the option is changed here because it is the earliest convienience that it is called.
if kexact_order == 0
  kexact_type = 'kexact';
end

if strcmp(kexact_type,'kexact')
  cnt = 1;
  for j = 0:kexact_order
    for i = 0:kexact_order-j
      p(1,cnt) = i;
      p(2,cnt) = j;
      cnt = cnt + 1;
    end
  end
  cell.nunknowns = size(p,2);

elseif strcmp(kexact_type,'kexact_extended')
  cnt = 1;
  for j = 0:kexact_order
    for i = 0:kexact_order
      p(1,cnt) = i;
      p(2,cnt) = j;
      cnt = cnt + 1;
    end
  end
  cell.nunknowns = size(p,2);
else
  error('Incorrect kexact_type!')
end


% Creates a row vector Ai for a given x
Ai_row = @(x) [x(1).^p(1,:).*x(2).^p(2,:)];

if (kexact_order == 0)
  xcc_quad = [0.5,0.5];
  wcc = [1];
else
  [xcc_quad, wcc] = sparse_grid(dim, dim*kexact_order, method);
end

for n = 1:cell.ncells
  nnodes = cell.nodes(n,1);
  if (nnodes == 3)
    % triangle quadrature
    % Transform the quadrature from a quad to a triangle 
    for i = 1:length(wcc)
      xcc(i,:) = transform_quad_to_triangle(xcc_quad(i,:),...
                                     quad_range_option_for_tri_transform);
    end

  elseif (nnodes == 4)
    % quad quadrature
    xcc = xcc_quad;
  end

  for nq = 1:length(wcc)
    xq(1,1) = cell.map(n).x(xcc(nq,:));
    xq(1,2) = cell.map(n).y(xcc(nq,:));

    A(nq,:) = Ai_row(xq)*wcc(nq);
  end
  Ai = sum(A).*cell.volume(n);
  cell.recon(n).Ai = Ai;
end
