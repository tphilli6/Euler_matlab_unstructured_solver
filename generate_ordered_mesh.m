function tip = generate_ordered_mesh(x,grid_type)
% This function generates a mesh and connectivity based on the ordered
% input of x which is sized [imax,jmax]


imax = size(x,1);
jmax = size(x,2);

%Define general grid type -----------------------------------------------------
%grid_type=0;%quads = 0, triangles = 1, mixed = 2

% cell type is an array to define a mixed cell type later if needed
if (grid_type == 0)
  ncells = (imax-1)*(jmax-1);
  cell.cell_type(1:ncells) = 0;

elseif (grid_type == 1)
  ncells = (imax-1)*(jmax-1)*2; %For triangles
  cell.cell_type(1:ncells) = 1;

elseif (grid_type == 2)
  nquads = (imax-1)*(jmax-1)/2;
  ntri   = (imax-1)*(jmax-1);
  ncells = nquads + ntri;
  cell.cell_type(1:nquads) = 0;
  cell.cell_type(nquads+1:ncells) = 1;
end


% Some parameters and useful functions
max_nsize = 5; % max size to define cell connectivity + 1
ij_to_vec = @(i,j) (j-1)*imax + i; %i,j location to cell vector location




% Write out connectivity
cnt = 1;
cell.nodes = zeros(ncells,max_nsize);
cell.size = 0;
for j = 1:jmax-1
  for i = 1:imax-1

   % Simple quads
   if (cell.cell_type(cnt) == 0)
     cell.nodes(cnt,1) = 4;
     cell.nodes(cnt,2:5) =[ij_to_vec(i,j),...
                            ij_to_vec(i+1,j),...
                            ij_to_vec(i+1,j+1),...
                            ij_to_vec(i,j+1) ];
     cell.type(cnt) = 9; %quad
     cnt = cnt + 1;

   elseif (cell.cell_type(cnt) == 1)
     % Triangle 1
     cell.nodes(cnt,1) = 3;
     cell.nodes(cnt,2:4) =[ij_to_vec(i,j+1),...
                                  ij_to_vec(i,j),...
                                  ij_to_vec(i+1,j)];
     cell.type(cnt) = 5; %triangle
     cnt = cnt + 1;

     % Triangle 2
     cell.nodes(cnt,1) = 3;
     cell.nodes(cnt,2:4) = [ij_to_vec(i,j+1), ...
                                   ij_to_vec(i+1,j),...
                                   ij_to_vec(i+1,j+1) ];
     cell.type(cnt) = 5; %triangle
     cnt = cnt + 1;

   end


  end
end


tip = cell.nodes(:,2:end);