function [resid] = compute_residual(cell, face, flux)
% Compute the residual
% Inputs :
%          cell : cell structure required for sizing
%          face : face structure with left and right states
%          flux : anonymous function used to compute the flux
%
% Outputs :
%           resid : residual array of size cell.soln


% Initialize residual
for i = 1:size(cell.mms_source,2)
  resid(:,i) = -cell.mms_source(:,i).*cell.volume';
end


for n = 1:numel(face)

  f = flux( face(n).ul, face(n).ur, face(n).normal )*face(n).area;

  i = face(n).cell_plus;
  resid(i,:) = resid(i,:) + f;

  % If there is a neighboring cell
  i = face(n).cell_neg;
  if (i~=-1)
    resid(i,:) = resid(i,:) - f;
  end


end


