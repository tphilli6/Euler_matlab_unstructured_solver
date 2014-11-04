function [vertex, face, cell] = generate_mesh(imax, jmax, grid_type, neq)

x0 = 0;
xL = 1;
y0 = 0;
yL = 1;

%imax = 5;
%jmax = 5;


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

% Customize the grid ----------------------------------------------------------


% Generate grid points  -------------------------------------------------------
% Cartesian
x = linspace(x0,xL,imax);
y = linspace(y0,yL,jmax);

% Some parameters and useful functions
max_nsize = 5; % max size to define cell connectivity + 1
ij_to_vec = @(i,j) (j-1)*imax + i; %i,j location to cell vector location

%% Build unstructured grid ----------------------------------------------------

% Store cell vertex
cnt = 1;
for j = 1:jmax
  for i = 1:imax
    vertex(cnt,:) = [x(i),y(j),0];
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

%% Remove non-unique faces
% This method loops through each face and for each face loops through all faces to see if the current face matches another face. This method is terribly slow and scales at the rate of Nnodes^2. You shouldn't use this method

%% Initialize the counter for the number of faces per cell
%face_list_old = face_list;
%cell.nface = zeros(cell.ncells,1);
% nface = 1;
% for nf = 1:nfaces_total
%   % If face has not been removed
%   if face_list(nf,:) ~= [-1, -1]
%
%     % Store the face and relevant neighbor data
%     face(nface).nodes(1:2)  = face_list(nf,:);%Stores the node index that makes of the face
%     face(nface).cell_plus = cell_list(nf); %Stores the cell which the face is a member of
%     face(nface).cell_neg  = -1; %blank neighbor
%
%     % Count the number of faces for each cell 
%     cell.nface(face(nface).cell_plus) = cell.nface(face(nface).cell_plus) + 1;
%     % Stores the face index for each cell : cell.faces( cell index, face counter index)
%     cell.faces( face(nface).cell_plus, cell.nface(cell_list(nf)) ) = nf;
%
%% Add face list to cells
%     % Compute geometric data
%      ds_vec = vertex(face(nface).nodes(2),:) - vertex(face(nface).nodes(1),:);
%      face(nface).area = sqrt( sum(ds_vec.^2) );
%      face(nface).normal(1:2) = [ ds_vec(2), ds_vec(1) ]/face(nface).area;
% 
% 
%     for face_check = 1:nfaces_total 
% 
%       % If face matches another face 
%       %  -Flipped because counter clockwise ordering results from other faces
%       if ( face_list(nf,:) == fliplr(face_list(face_check,:)) )
%         % It's a match to another face!! So get rid of the match
%         face_list(face_check,:) = -1;
%     
%         % Store the face and relevant cell connectivity
%         face(nface).cell_neg  = cell_list(face_check); %Stores the neighboring cell
% 
%         % Count the number of faces for each cell 
%         cell.nface(face(nface).cell_neg) = cell.nface(face(nface).cell_neg) + 1;
%         % Stores the face index for each cell : cell.faces( cell index, face counter index)
%         cell.faces( face(nface).cell_neg, cell.nface(cell_list(nf)) ) = nf;
%       end
% 
%     end
%
%     nface = nface + 1;
%   end 
% end
%%for n = 1:unique_faces
%%fprintf('%8.0f  %8.0f\n',face(n).nodes)
%%end
%unique_faces = nface-1;


% Remove non-unique faces
% This method sorts all nodes for each face in ascending order, then the rows of each face are sorted. If there are two faces with the same two nodes then they will end up together. This algorithm removes the second for loop the previous method used by only checking one face ahead of the current face for repeats. This method scales by a factor of Nnodes.

% Initialize the counter for the number of faces per cell
face_list_sort = sort(face_list, 2, 'ascend'); % order node index per cell in ascending order
face_list_sort = sortrows(face_list_sort);
face_list = face_list_sort;

cell.nface = zeros(cell.ncells,1);
 nface = 1;
 for nf = 1:nfaces_total
   % If face has not been removed
   if face_list(nf,:) ~= [-1, -1]

     % Store the face and relevant neighbor data
     face(nface).nodes(1:2)  = face_list(nf,:);%Stores the node index that makes of the face
     face(nface).cell_plus = cell_list(nf); %Stores the cell which the face is a member of
     face(nface).cell_neg  = -1; %blank neighbor

     % Count the number of faces for each cell 
     cell.nface(face(nface).cell_plus) = cell.nface(face(nface).cell_plus) + 1;
     % Stores the face index for each cell : cell.faces( cell index, face counter index)
     cell.faces( face(nface).cell_plus, cell.nface(cell_list(nf)) ) = nf;

     % Add face list to cells
     % Compute geometric data
      ds_vec = vertex(face(nface).nodes(2),:) - vertex(face(nface).nodes(1),:);
      face(nface).area = sqrt( sum(ds_vec.^2) );
      face(nface).normal(1:2) = [ ds_vec(2), ds_vec(1) ]/face(nface).area;
 
 
 
       % If face matches one face ahead, remove the one face ahead from the check
       face_check = nf + 1;
       if (face_check<nfaces_total)
         if all( face_list(nf,:) == face_list(face_check,:) )
           % It's a match to another face!! So get rid of the match
           face_list(face_check,:) = -1;
     
           % Store the face and relevant cell connectivity
           face(nface).cell_neg  = cell_list(face_check); %Stores the neighboring cell
 
           % Count the number of faces for each cell 
           cell.nface(face(nface).cell_neg) = cell.nface(face(nface).cell_neg) + 1;
           % Stores the face index for each cell : cell.faces( cell index, face counter index)
           cell.faces( face(nface).cell_neg, cell.nface(cell_list(nf)) ) = nf;
         end
       end
 
     nface = nface + 1;
   end 
 end
unique_faces = nface-1



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
end
total_volume = sum(cell.volume);

size(cell.volume);
%% Write grid to file %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
write_vtk(vertex, cell);
fid = fopen('grid.vtk','a');
fprintf(fid,'CELL_DATA %8.0f\n',length(cell.volume));
fclose(fid);
write_vtk_data(cell.volume, 'volume');


for n = 1:length(face)
  face(n).ul = zeros(1,neq);
  face(n).ur = zeros(1,neq);
end

