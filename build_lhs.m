function [A] = build_lhs(stencil, reconstruction, fit_type)
% build_lhs constructs the left hand side in the linear system A*x=b
% Requires:
%           reconstruction(n).Axc  : contains a reference point
%           reconstruction(n).Ai   : contains the integral over the cell for a given quadrature
%           stencil                : vector of cells included in the stencil
%           fit_type               : special case if fit_type = 'lsq'


for i = 1:length(stencil);
  ncell = stencil(i);
  % Combine the rows associated with each cell in the stencil related 
  % to the integral of the given polynomial over the cell.
  A(i,:) = reconstruction(ncell).Ai;

  % A least squares fit with (n+1) cells where n is the number of unknowns
  % This is a constrained lsq fit where the fit is constrained only for the 
  % current cell. Must also do this for the right hand side b-b_1. 
  % The first coefficient is equal to u_i
  if strcmp(fit_type,'lsq')
    A(i,:) = A(i,:) - reconstruction(stencil(1)).Axc;
    %A(i,:) = A(i,:) - A(1,:);    
  end

end


