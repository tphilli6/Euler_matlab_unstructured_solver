function [vertex, cell, face] = compute_grid_derived_data(x,y, grid_type, ti)



% Store cell vertex
vertex = [ reshape(x,[numel(x),1]), reshape(y, [numel(y),1]), zeros(numel(x),1) ];


if nargin<4 % tip not set
    ti = generate_ordered_mesh(x,grid_type); 
end
    

for i = 1:size(ti,1)
    I = find(ti(i,:)~=0);
    cell.nodes(i,:) = [length(I),ti(i,:)];
    if length(I) == 3
        cell.type(i) = 5;
        cell.cell_type(i) = 1;
    elseif length(I) == 4
        cell.type(i) = 9;
        cell.cell_type(i) = 0;
    end
end
cell.ncells = size(ti,1);
cell.vtk_size = sum(cell.nodes(:,1)) + cell.ncells;
max_nsize = 5;


%% Loop through connectivity data to assemble faces
% I wrote a more efficient routine to extract the cell and face
% connectivity but I don't feel like making it work with the new routine...
% [face_list, cell_faces, face_cells] = find_neighbors(ti);

nface = 1;
face.nodes = []; % Initialize but not in memory
all_faces = zeros(size(ti));
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

