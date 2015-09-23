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
%   error('Incorrect kexact_type!')
end

% HARDWIRE : curtis-clenshaw quadrature for dirichlet bc 
% Stores the 1D quadrature for flux integration
[xquad, wquad] = gauss_patterson( flux_integral_order );

cell.reconstruction_param.xquad = xquad;
cell.reconstruction_param.wquad = wquad;

[cell.reconstruction_param.moment] = compute_reconstruction_moments(vertex, cell, face, p );

% This is a logical. If Ai is changed as it is in this script then when reconstruct_solution is called, the lhs is rebuild
cell.lhs_set = 0;


cell.reconstruction_param.px = p(1,:);
cell.reconstruction_param.py = p(2,:);
% cell.reconstruction.p = p;
