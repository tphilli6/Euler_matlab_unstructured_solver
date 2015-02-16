function [vertex, cell, face] = compute_grid_derived_data(x,y, grid_type)

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



% Store cell vertex
cnt = 1;
for j = 1:jmax
  for i = 1:imax
    vertex(cnt,:) = [x(i,j),y(i,j),0];
    cnt = cnt + 1;
  end
end


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
cell.ncells = cnt-1;
cell.vtk_size = sum(cell.nodes(:,1)) + cell.ncells;



%% Loop through connectivity data to assemble faces
nface = 1;
face.nodes = []; % Initialize but not in memory
all_faces = zeros(cell.ncells,max_nsize);
face_cnt=1;
for nc = 1:cell.ncells % List out faces -simply adds initial face to list to close the polygon
  all_faces(nc,1:cell.nodes(nc,1)+1) = [cell.nodes(nc,2:cell.nodes(nc,1)+1),cell.nodes(nc,2)];
  nfaces(nc) = cell.nodes(nc,1);
  for nf = 1:nfaces(nc)
    % Create an array of faces [nfaces x 2]
    face_list(face_cnt,:) = [all_faces(nc,nf:nf+1)];
    % Store the cell that the face is attached to 
    cell_list(face_cnt) = nc;
   
    % Increment face counter
    face_cnt = face_cnt + 1;
  end 
end
nfaces_total = sum(nfaces);
face_list;


% Remove non-unique faces
% This method sorts all nodes for each face in ascending order, then the rows of each face are sorted. If there are two faces with the same two nodes then they will end up together. This algorithm removes the second for loop the previous method used by only checking one face ahead of the current face for repeats. This method scales by a factor of Nnodes.

% Initialize the counter for the number of faces per cell
face_list_unsorted = face_list;
face_list_sort = sort(face_list, 2, 'ascend'); % order node index per cell in ascending order
[face_list_sort, Isort] = sortrows(face_list_sort);
face_list = face_list_sort;

cell.nface = zeros(cell.ncells,1);
cell.nnbr = zeros(cell.ncells,1);
 nface = 1;
 for nf = 1:nfaces_total
   % If face has not been removed
   if face_list(nf,:) ~= [-1, -1]

     % Store the face and relevant nbr data
     %face(nface).nodes(1:2)  = face_list(nf,:);%Stores the node index that makes of the face
     face(nface).nodes(1:2)  = face_list_unsorted(Isort(nf),:);%Stores the node index that makes of the face
     face(nface).cell_plus = cell_list(Isort(nf)); %Stores the cell which the face is a member of
     face(nface).cell_neg  = -1; %blank nbr

     % Count the number of faces for each cell 
     cell.nface(face(nface).cell_plus) = cell.nface(face(nface).cell_plus) + 1;
     % Stores the face index for each cell : cell.faces( cell index, face counter index)
     cell.faces( face(nface).cell_plus, cell.nface(face(nface).cell_plus) )=nface;

     % Add face list to cells
     % Compute geometric data
      ds_vec = vertex(face(nface).nodes(2),:) - vertex(face(nface).nodes(1),:);
      face(nface).area = sqrt( sum(ds_vec.^2) );
      face(nface).normal(1:2) = [ ds_vec(2), -ds_vec(1) ]/face(nface).area;
 
 
 
       % If face matches one face ahead, remove the one face ahead from the check
       face_check = nf + 1;
       if (face_check<nfaces_total)
         if all( face_list(nf,:) == face_list(face_check,:) )
           % It's a match to another face!! So get rid of the match
           face_list(face_check,:) = -1;
     
           % Store the face and relevant cell connectivity
           face(nface).cell_neg  = cell_list(Isort(face_check)); %Stores the nbring cell
 
           % Count the number of faces for each cell 
           cell.nface(face(nface).cell_neg) = cell.nface(face(nface).cell_neg) + 1;
           % Stores the face index for each cell : cell.faces( cell index, face counter index)
           cell.faces( face(nface).cell_neg, cell.nface(face(nface).cell_neg))=nface;

           % get the current and nbr cell
           cur_cell = cell_list(Isort(nf));
           nbr_cell = cell_list(Isort(face_check));

           % store current nbr count and increment by 1 for current 
           % cell
           nnbr = cell.nnbr(cell_list(Isort(nf)));
           nnbr = 1 + nnbr;
           cell.nnbr(cur_cell) = nnbr;
           cell.nbrs(cur_cell,nnbr) = nbr_cell;

           % store current nbr count and increment by 1 for nbr
           % cell
           nnbr = cell.nnbr(cell_list(Isort(face_check)));
           nnbr = 1 + nnbr;
           cell.nnbr(nbr_cell) = nnbr;
           cell.nbrs(nbr_cell,nnbr) = cur_cell;

         end
       end
 
     nface = nface + 1;
   end 
 end
unique_faces = nface-1;


% Compute cell data -----------------------------------------------------------
for n = 1:cell.ncells
  if (cell.cell_type(n) == 0)
    % compute quad area using two half triangle areas
    cell.volume(n) = triangle_area(vertex(cell.nodes(n,2),:), ...
                                   vertex(cell.nodes(n,3),:), ...
                                   vertex(cell.nodes(n,4),:) );

    cell.volume(n) = triangle_area(vertex(cell.nodes(n,3),:), ...
                                   vertex(cell.nodes(n,4),:), ...
                                   vertex(cell.nodes(n,5),:) ) ...
                   + cell.volume(n);

  elseif (cell.cell_type(n) == 1)
    % compute tri area
    cell.volume(n) = triangle_area(vertex(cell.nodes(n,2),:), ...
                                   vertex(cell.nodes(n,3),:), ...
                                   vertex(cell.nodes(n,4),:) );

  end

  % Compute the cell center
  nnodes = cell.nodes(n,1);
  cell.xc(n,:) = zeros(1,3);
  for i = 1:nnodes
    cell.xc(n,:) = cell.xc(n,:) + vertex(cell.nodes(n,i+1),:);
  end
  cell.xc(n,:) = cell.xc(n,:)/nnodes;


end

