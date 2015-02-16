function [resid] = compute_continuous_residual(cell, face, flux)
% Compute the residual
% Inputs :
%          cell : cell structure required for sizing
%          face : face structure with left and right states
%          flux : anonymous function used to compute the flux
%
% Outputs :
%           resid : residual array of size cell.soln

wquad = cell.reconstruction_param.wquad;

% Initialize residual
for i = 1:size(cell.mms_source,2)
  resid(:,i) = -cell.mms_source(:,i).*cell.volume';
end


for n = 1:numel(face)

  fl = 0;
  for j = 1:length(wquad)
    fl = fl ...
      +flux( face(n).ul(j,:), face(n).normal )*wquad(j);
  end
  fl = fl*face(n).area;

  i = face(n).cell_plus;
  resid(i,:) = resid(i,:) + fl;

  % If there is a neighboring cell compute the flux through the right face
  i = face(n).cell_neg;
  if (i>=1)

    fr = 0;
    for j = 1:length(wquad)
      fr = fr ...
        +flux( face(n).ul(j,:), face(n).normal )*wquad(j);
    end
    fr = fr*face(n).area;
      resid(i,:) = resid(i,:) - fr;

  end


end

%for n = 1:size(resid,2)
%  resid(:,n) = resid(:,n)./cell.volume'; 
%end

