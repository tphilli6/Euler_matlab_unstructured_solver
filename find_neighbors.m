function [face_nodes, cell_faces, face_cells] = find_neighbors(cell)

j = [1:size(cell,2),1];
faces = cell(:,1:2);
for i = 2:size(cell,2)
    for ii = 1:size(cell,1)
        if cell(ii,j(i+1))~=0 && cell(ii,j(i))~=0
          faces = [faces; [cell(ii,j(i)), cell(ii,j(i+1))] ];
        end
    end
end

face_sort = sort(faces,2);

[face_nodes, IA, IC] = unique(face_sort,'rows');


cell_faces = reshape(IC,size(cell));%[IC(1:n:end),IC(2:n:end)]';

if nargout == 3
%     n = size(cell,1);
%     cell_neighbor = zeros(size(cell));
%     nneighbor=zeros(n,1);
%     for i = 1:length(IA)
%             I = find(IA(i)==cell_faces);
%             if (length(I) > 1)
%                 cell1 = rem(I(1),n);
%                 cell2 = rem(I(2),n);
%                 
%                 if cell1==0; cell1=n; end
%                 if cell2==0; cell2=n; end
%                 
%                 nneighbor(cell1) = nneighbor(cell1)+1;
%                 cell_neighbor(cell1,nneighbor(cell1)) = cell2;
%                 
%                 nneighbor(cell2) = nneighbor(cell2)+1;
%                 cell_neighbor(cell2,nneighbor(cell2)) = cell1;
%             end
% 
%     end 
    
    face_cells = zeros(length(IA),2);
    cell_cnt = zeros(length(IA),1);
    n = size(cell,1);
    for i = 1:length(IA)
        I = find(i==cell_faces);
          for j = 1:length(I)
              c = rem(I(j),n);
              if c==0; c=n; end
              
              face_cells(i,j) = c;
          end
        
    end
    
    
end