function [dt] = euler_time_step(vertex, cell, face, CFL, glb)
% Computes the local time step for each cell
% Inputs:
%         vertex : list of cell verticies (do I need this?)
%         cell   : cell data structure
%         face   : face data structure
%         CFL    : You should know what this
%         glb    : logical (glb=1 uses global time step)

% pseudo code
% loop over cells
% compute speed of sound for the cell
% compute average face normals
% compute average face areas
% compute time step
%   dt = cfl * volume/ [ ( |vel . normal| + a)* face area  ]



for n = 1:cell.ncells
  a = sound_speed( cell.soln(n,4), cell.soln(n,1) );

  faces = cell.faces(n,:);
  % Cell boundary area
  % Hino, T., Martinelli, L., Jameson, A., "A FINITE-VOLUME METHOD WITH UNSTRUCTURED GRID FOR FREE SURFACE FLOW SIMULATIONS," Presented at the sixth international conference on numerical ship hydrodyamics, 1993.
  S = 0;
  for nn = 1:length(faces)
    S = S + face(nn).area;
  end

  % CFL * S / ( |V|+a )
  dt(n) = CFL*S / ( sqrt(cell.soln(n,2)^2 + cell.soln(n,3)^2) + a );

end

if (glb) 
  dt = min(dt)*ones(size(dt));
end
