function [face, cell] =  compute_face_data_from_triangulation(ti,vertex)

if (size(vertex,2)==2)
    Vtx = [vertex,zeros(size(vertex(:,1)))];
else
    Vtx = vertex;
end

[face_nodes, cell_faces, face_cells] = find_neighbors(ti);

cell.ncells = size(ti,1);
cell.vtk_size = cell.ncells + numel(ti);
cell.nnbr = zeros(cell.ncells,1);
cell.nbrs = zeros(cell.ncells,1);
for i = 1:size(ti,1)

    
    nc = length(ti(i,:));
    I = find(ti(i,:)~=0); nc = length(I);
    cell.nodes(i,1)      = nc;
    cell.nodes(i,2:nc+1) = ti(i,I);
    cell.xc(i,:) = mean(vertex(ti(i,I),:),1);
    
    I = find( cell_faces(i,:) ~= 0);
    cell.nface(i,1) = length(I);
    cell.faces(i,1:length(I)) = cell_faces(i,I);

    % Compute cell shape related parameters
    if (nc==3) %triangle
        cell.type(i) = 5;
        cell.cell_type(i) = 1;
        cell.volume(i) = triangle_area(Vtx(ti(i,1),:), ...
                                       Vtx(ti(i,2),:), ...
                                       Vtx(ti(i,3),:) );
    elseif (nc==4) %quad
        cell.type(i) = 9;
        cell.cell_type(i) = 0;
        
        % compute quad area using two half triangle areas
        cell.volume(i) = triangle_area(Vtx(ti(i,1),:), ...
                                       Vtx(ti(i,2),:), ...
                                       Vtx(ti(i,3),:) );

        cell.volume(i) = triangle_area(Vtx(ti(i,2),:), ...
                                       Vtx(ti(i,3),:), ...
                                       Vtx(ti(i,4),:) ) ...
                       + cell.volume(n);
    else %poly
        cell.type(i) = 15;
        cell.cell_type(i) = 2; 
        
        cell.volume(i) = 0;
    end
    
    
end



for i = 1:length(face_cells)
         face(i).nodes = face_nodes(i,:);
         fn = face_nodes(i,:);
         
         if (face_cells(i,2)==0);
             face_cells(i,2) = face_cells(i,1);
             face_cells(i,1) = -1;
            
             %Check to see if neg/pos nodes are correct
             cp = face_cells(i,2);
             cn = ti(cp,:);
             Icn = find(cn~=0);
             cn = [cn(Icn),cn(1)];
             In = find(fn(1) == cn,1);
             if (cn(In+1) ~= fn(2))
                 face_nodes(i,:) = fliplr(face_nodes(i,:));
             end
             
         else
             %Check to see if neg/pos nodes are correct
             cp = face_cells(i,2);
             cn = ti(cp,:);
             Icn = find(cn~=0);
             cn = [cn(Icn),cn(1)];
             In = find(fn(1) == cn,1);
             if (cn(In+1) ~= fn(2))
                 face_cells(i,:) = fliplr(face_cells(i,:));
             end
             
             % Add cell neighbors based on cp and cn
             cp = face_cells(i,2);
             cell.nbrs(cp,cell.nnbr(cp)+1) = face_cells(i,1);
             cell.nnbr(cp) = cell.nnbr(cp)+1;
             
             cn = face_cells(i,1);
             cell.nbrs(cn,cell.nnbr(cn)+1) = face_cells(i,2);
             cell.nnbr(cn) = cell.nnbr(cn)+1;
             
             
             
         end
         
         face(i).cell_plus = face_cells(i,2);
         face(i).cell_neg  = face_cells(i,1);
         
         ds_vec = vertex(face_nodes(i,2),:) - vertex(face_nodes(i,1),:);
         face(i).area = sqrt( sum(ds_vec.^2) );
         face(i).normal(1:2) = [ ds_vec(2), -ds_vec(1) ]/face(i).area;
     end