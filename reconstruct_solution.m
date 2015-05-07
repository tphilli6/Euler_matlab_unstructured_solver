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
    [Ainv, cond, A, wij] = reconstruction_lhs(stencil, n, cell.reconstruction_param, cell.xc, fit_type, 0, ge_terms);
%     [Ainv2, ~, A2] = ts_2d_reconstruction(cell, n, 2);
    
    reconstruction(n).Ainv = Ainv;
    reconstruction(n).A = A;
    reconstruction(n).wij = wij;
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
        
%       rhs = b(:,j);
%       for jj = 1:ge_terms
%         for ii = j+1:size(rhs,1);
%           rhs(ii,1) = rhs(ii,1) - reconstruction(n).wij(ii,jj)*rhs(jj,1);
%         end
%       end
      
      rhs = (b(:,j)-b(isten,j)).*reconstruction(n).wij;
      rhs(isten) = b(isten,j);
    else 
      rhs = b(:,j);
    end

    coef_temp = reconstruction(n).Ainv*rhs;

    %       face(n).ul(j,:) = ( (x(1)-cell.xc(l_cell,1)).^px.*(x(2)-cell.xc(l_cell,2)).^py )*coef;
%      ( (x(1)).^px.*(x(2)).^py )*coef;

%     if strcmp(fit_type,'lsq')
% % %      coef(1) = coef(1) + cell.soln(n,j);
% % %       coef_adjust = cell.reconstruction(n).Ai(2:end)*coef_temp;
%       coef(:,j) = [b(isten,j);coef_temp];
% % %       coef(:,j) = [cell.soln(n,j)-coef_adjust; coef_temp];
% % %    disp(reconstruction(n).coef(:,j)')
%     else

      coef(:,j) = coef_temp;
%     end



  end
      reconstruction(n).coef = coef;
% reconstruction(n).coef(:,1)
end
