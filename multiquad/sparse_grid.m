function [ x, w ] = sparse_grid( dim, order, method )
%% -----------------------------------------------------------------------------
%SPARSE_GRID Generates quadrature points and weights for Smolyak rule.
%   
%   Given a univariate quadrature rule (defined by order and method),
%   this routine produces the sparse grid quadrature points using the Smolyak
%   construction.  If the univariate rule generates nested points, then the
%   resulting multivariate quadrature rule will have sparsely defined points.
%
%   Usage:  [ qx, qw ] = sparse_grid( dim, order, method ); 
%
%   Variables:
%       qx       -  quadrature points   ( output )
%       qw       -  quadrature weights  ( output )
%    
%       dim      -  dimension of integration domain  ( input )
%       order    -  order of accuracy, the Smolyak rule will integrate monomials
%                   exactly up to and including this order  ( input )
%       method   -  quadrature method   ( optional input )
%                   default is to use Clenshaw-Curtis nested rules
%                   'CC' - Clenshaw-Curtis
%                   'GL' - Gauss-Legendre
%                   'GP' - Gauss-Patterson
%                   'method.x, method.w' a user-defined quadrature rule
%
%% -----------------------------------------------------------------------------

  if ( nargin<3 )
    method = 'CC';
  end

  %% Set up the univariate rule
  %
  if ( isstruct(method) )
    x1 = method.x;
    w1 = method.w;
    
  elseif ( strcmp(method,'CC') )
    for k=1:20
      [ x1{k}, w1{k} ] = curtis_clenshaw( k );
    end
    
  elseif ( strcmp(method,'GL') )
    for k=1:25
      [ x1{k}, w1{k} ] = gauss_legendre( k ); 
    end
    
  elseif ( strcmp(method,'GP') )
    for k=1:24  % loop over "rules" the order is 2*k+1
      [ x1{k}, w1{k} ] = gauss_patterson( k );
    end
    
  else
    error('SPARSE_GRID: unknown or improperly defined method')
    
  end
  clear method
  
  %% Create blank list of Smolyak quadrature points and weights, then collect
  %  points and weights for Smolyak rule
  x = [];
  w = [];

  minq = max(0,order-dim);
  maxq = order-1;
  for q = minq:maxq
    bq = (-1)^(maxq-q) * nchoosek(dim-1,order-1-q);
    
    %  Calculate the N_q^D range for 1d tensor products
    ii = NqD(q,dim);
    for j=1:size(ii,1)   % add the jth term in the sum to the list
      for k=1:dim
        rules.x{k} = x1{ii(j,k)};
        rules.w{k} = w1{ii(j,k)};
      end
      [xt, wt] = tensor_product( rules );
      
      x = [x;    xt];
      w = [w; bq*wt];
    end
  end

  %% Group together duplicated quadrature points (consolidate their weights)
  %  to eliminate unnecessary function evaluations
  [x index] = sortrows(x);
  w         = w(index);
  uniq      = 1; 
  lastuniq  = 1;

  %  sum weights associated with duplicated quadrature points
  for j=2:size(x,1)
    if ( norm(x(j,:)-x(j-1,:))<1e-7 )  % duplicate point
      w(lastuniq) = w(lastuniq) + w(j);
    else
      lastuniq = j;
      uniq     = [uniq ; j ];
    end
  end
  
  %  only keep nodes with collected weights
  x = x(uniq,:);
  w = w(uniq  );

end


function [ x, w ] = tensor_product( rules )
%% -----------------------------------------------------------------------------
%  TENSOR_PRODUCT - For a set of univariate rules, compute tensor product rule.
%
%  Usage:  [ x, w ] = tensor_product( rules )
%
%  Variables:  rules  - a cell array containing the univariate rule for each
%                       dimension
%
%              x      - a list of quadrature points (# points, dim)
%              w      - corresponding quadrature weights
%
%              dim    - the dimension of the domain to be integrated
%
%% -----------------------------------------------------------------------------
  x = rules.x{1} ; w = rules.w{1};
  for j=2:length(rules.x)
    jdim_nodes = rules.x{j};
    x = [kron(x,ones(size(jdim_nodes,1),1)) kron(ones(size(x,1),1),jdim_nodes)];
    w = kron(w,rules.w{j});
  end
end




