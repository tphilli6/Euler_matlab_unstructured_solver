function [face_nodes, cell_faces, face_cells] = find_neighbors(cell)

% j = [1:size(cell,2),1];
% faces= [];
% for i = 1:size(cell,1)
%     I = find(cell(i,:)~=0);
%     n1 = [cell(i,I)';cell(i,1)];
%     faces = [faces; [n1(1:end-1), n1(2:end)]];
% end



% j = [1:size(cell,2),1];
% faces2 = cell(:,1:2);
% for i = 2:size(cell,2)
%     for ii = 1:size(cell,1)
%         if cell(ii,j(i+1))~=0 && cell(ii,j(i))~=0
%           faces2 = [faces2; [cell(ii,j(i)), cell(ii,j(i+1))] ];
%         end
%     end
% end
% 
% face_sort2 = sort(faces2,2);
% 
% [face_nodes2, IA2, IC2] = unique(face_sort2,'rows');

%Cell mapping
Icell = reshape( repmat([1:size(cell,1)]',[1,size(cell,2)])',[numel(cell),1]);

cell2 = [cell, cell(:,1)];
faces = [reshape(cell2(:,1:end-1)',[numel(cell2(:,1:end-1)),1]),...
         reshape(cell2(:,2:end)',[numel(cell2(:,2:end)),1])];
I1 = find(faces(:,1)~=0);
I2 = find(faces(:,2)~=0);
faces = [faces(I1,1), faces(I2,2)];
Icell = Icell(I2);
face_sort = sort(faces,2);

[face_nodes, IA, IC] = unique(face_sort,'rows');

% if any(any(face_nodes2~=face_nodes))
%     faces2
% end

cell_faces = zeros(size(cell'));
Inonzero = find(cell'~=0);
cell_faces(Inonzero) = IC;
cell_faces=cell_faces';
% cell_faces2 = reshape(IC2,size(cell));%[IC(1:n:end),IC(2:n:end)]';

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
    for i = 1:max(max(cell_faces))
        I = find(i==cell_faces);
          for j = 1:length(I)
              c = rem( I(j), n);
%               c = rem(I(j),n);
              if c==0; c=n; end
              %Icell(IA);% Cell number
              face_cells(i,j) = c;
          end
        
    end
    
    
end