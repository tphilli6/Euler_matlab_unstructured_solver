function [l2norms, converged, l2_to_normalize] = check_convergence( resid, iter, l2_to_normalize, toler )
% checks the l2 norms of the iterative residual to see if they are below the convergence tolerance
% Inputs :
%          resid : residual arrau of size soln
%          iter  : current iteration
%l2_to_normalize : l2 norms used to normalize the residual
%          toler : convergence tolerance
%
% Outputs :
%           l2norms : the norms of the residual
%         converged : convergence flux (1=converged, 0=unconverged)
%   l2_to_normalize : l2 norm to normalize residual

l2norms = sqrt( sum( resid.^2 ) / size(resid,1) );

if (iter == 5)
  l2_to_normalize = l2norms;
end

% Normalizes residual
l2norms = l2norms./l2_to_normalize;

% Sets converged flag
if (l2norms < toler)
  converged = 1;
else
  converged = 0;
end


