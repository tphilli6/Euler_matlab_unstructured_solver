%setup reconstruction for each cell
% Requires the definition of 
% kexact_order : order of reconstruction
% kexact_type  : 'kexact', 'kexact-extended'

% This is a script that will setup the row in the left hand side which is related to the discrete integral a polynomial of order k 
% cell.reconstruction(n).Arow

% For quads, the quadrature used is curtis clenshaw. For each control point a function evaluation can be computed for a given polynomial
% p(x) = a+bx+cy+dxy which in matrix form is 
% Ai*coef = ui
%and
% Ai = [1, map(i).x(xcc), map(i).y(xcc), map(i).x(xcc)*map(i).y(xcc),...]
clear p A


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
  fit_type    = 'kexact';
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
elseif strcmp(kexact_type,'ts')
      cnt = 1;
      for n = 0:kexact_order
        for m = 0:n
          p(1,cnt) = m;
          p(2,cnt) = n-m;
          ts_coef(1,cnt) = 1/factorial(n)*factorial(n)/(factorial(m)*factorial(n-m));
          cnt = cnt + 1;
        end
      end
      cell.nunknowns = size(p,2);
    
    
else
  error('Incorrect kexact_type!')
end

% HARDWIRE : curtis-clenshaw quadrature for dirichlet bc 
% Stores the 1D quadrature for flux integration
%if (kexact_order == 0)
%  xquad=[0.5];
%  wquad=[1];
%else
  [xquad, wquad] = curtis_clenshaw( flux_integral_order );
%end
cell.reconstruction_param.xquad = xquad;
cell.reconstruction_param.wquad = wquad;

% Creates a row vector Ai for a given x
Ai_row = @(x) [x(1).^p(1,:).*x(2).^p(2,:)];

if (kexact_order == 0)
  xcc_quad = [0.5,0.5];
  wcc = [1];
else
  [xcc_quad, wcc] = sparse_grid(dim, dim*kexact_order, method);
end

for nn = 1:cell.ncells
  nnodes = cell.nodes(nn,1);
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
  else
    xcc = vertex_grid(nn,:); % Use the node of the regular grid.
      
  end

  if (nnodes ==3 || nnodes == 4);
      
      for nq = 1:length(wcc)
        xq(1,1) = cell.map(nn).x(xcc(nq,:));
        xq(1,2) = cell.map(nn).y(xcc(nq,:));

        A(nq,:) = Ai_row(xq)*wcc(nq);
      end
  else
      A(1,:) = Ai_row(xcc);
      
  end

  Ai = sum(A);%.*cell.volume(n);
  cell.reconstruction(nn).Ai = Ai;

  cell.reconstruction(nn).Axc = Ai_row(cell.xc(nn,1:2));
end

% This is a logical. If Ai is changed as it is in this script then when reconstruct_solution is called, the lhs is rebuild
cell.lhs_set = 0;


cell.reconstruction_param.px = p(1,:);
cell.reconstruction_param.py = p(2,:);
