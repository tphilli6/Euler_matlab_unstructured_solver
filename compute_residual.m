function [resid] = compute_residual(cell, face, flux, ncell)
% Compute the residual
% Inputs :
%          cell : cell structure required for sizing
%          face : face structure with left and right states
%          flux : anonymous function used to compute the flux
%
% Outputs :
%           resid : residual array of size cell.soln

if (nargin == 3)
    face_loop = 1:length(face);
elseif (nargin == 4)
    face_loop = cell.faces(ncell,:);
end

wquad = cell.reconstruction_param.wquad;

% Initialize residual
for i = 1:size(cell.mms_source,2)
  resid(:,i) = -cell.mms_source(:,i).*cell.volume';
end


for n = face_loop

  f = 0;
  for j = 1:length(wquad)
    f = f ...
      +flux( face(n).ul(j,:), face(n).ur(j,:), face(n).normal )*wquad(j);
  end
  f = f*face(n).area;

  i = face(n).cell_plus;
  resid(i,:) = resid(i,:) + f;

  % If there is a neighboring cell
  i = face(n).cell_neg;
  if (i>=1)
    resid(i,:) = resid(i,:) - f;
  end


end

%for n = 1:size(resid,2)
%  resid(:,n) = resid(:,n)./cell.volume'; 
%end

