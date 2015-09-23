function [reconstruction, lhs_set]=reconstruct_solution(cell, fit_type, ge_terms, nc)
% This function computes the reconstruction for a given cell. The polynomial exponents are set in the setup_reconstruction routine and depend on the variables kexact_type = 'kexact' as the original kexact reconstruction or 'kexact-extended' which includes more terms. There is a fit type variable that changes how the reconstruction is computed. This variable is set by fit_type = 'lsq' which is the tranditional lsq reconstruction for kexact and satisfies the conservation of the mean for the current cell, 'kexact' maintains conservation of the mean over all stencil cells and 'lsq-true' does not satisfy conservation of the mean over any stencil cells. (The latter cannot be used for iterative solutions on truncation error estimation)

% Inputs : cell = cell structure variables required 
%            px, py, ncells, recon.Ai, stencil.cells
%        : fit_type = determines the type of reconstruction 
%            'kexact', 'lsq', 'lsq-true'

if (nargin == 3)
    loop_cells = 1:cell.ncells;
elseif (nargin == 4)
    loop_cells = nc;
end


reconstruction = cell.reconstruction;

% Iterate over each cell
if (cell.lhs_set==0)
  fprintf('Setting up LHS for reconstruction.\n')
  for n = loop_cells

    % extract cells over the stencil
    stencil = cell.stencil(n).cells;

    % Build the left hand side
    [Ainv, cond, A, wij, Aeval] = reconstruction_lhs(stencil, n, cell.reconstruction_param, cell.xc, fit_type, 0, ge_terms);
%     [Ainv2, ~, A2] = ts_2d_reconstruction(cell, n, 2);
    
    reconstruction(n).Ainv = Ainv;
    reconstruction(n).A = A;
    reconstruction(n).wij = wij;
    reconstruction(n).Aeval = Aeval;
  end

  cell.reconstruction = reconstruction;
end
lhs_set = 1;

test_max = 0;
neq = length(cell.soln(1,:));
% Setup rhs and solve for the coefficients
for n = loop_cells

    % extract cells over the stencil
    stencil = cell.stencil(n).cells;

%     for i = 1:length(stencil);
%       ncell = stencil(i);
      % Combine the rows associated with each cell in the stencil related
      % to the integral of the given polynomial over the cell.
%       b(i,:) = cell.soln(ncell,:);
      b = cell.soln(stencil,:);
      isten = find(n==stencil);

%     end



    % A least squares fit with (n+1) cells where n is the number of unknowns
    % This is a constrained lsq fit where the fit is constrained only for the 
    % current cell. Must also do this for the right hand side b-b_1. 
    % The first coefficient is equal to u_i
  clear rhs
  for j = 1:neq
    if strcmp(fit_type,'lsq')
        
      rhs = (b(:,j)-b(isten,j)).*reconstruction(n).wij(:,1);
      rhs(isten) = b(isten,j);
      
      %HARDCODE - Removes mean constraint.
%       rhs = b(:,j).*reconstruction(n).wij(:,1);
    else 
        
      rhs = b(:,j);
    end

    coef_temp = reconstruction(n).Ainv*rhs;

    coef(:,j) = coef_temp;

  end
      reconstruction(n).coef = coef;

end
