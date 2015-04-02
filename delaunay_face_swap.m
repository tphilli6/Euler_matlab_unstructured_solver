function [ti, face_nodes, face_cells] = delaunay_face_swap(cell1, cell2, face, ti, face_nodes, face_cells, x, override)

          % override option to force a face swap even if already delaunay
          if nargin <8; override = 0; end


          nodes1 = ti(cell1,:);
          nodes2 = ti(cell2,:);

          % Find the faces associated with the first new cell
         [~,~,face11] = intersect(nodes1([1,2]), face_nodes,'rows');
         [~,~,face12] = intersect(nodes1([2,1]), face_nodes,'rows');
         face1p(1) = [face11;face12];

         [~,~,face11] = intersect(nodes1([2,3]), face_nodes,'rows');
         [~,~,face12] = intersect(nodes1([3,2]), face_nodes,'rows');
         face1p(2) = [face11;face12];

         [~,~,face11] = intersect(nodes1([3,1]), face_nodes,'rows');
         [~,~,face12] = intersect(nodes1([1,3]), face_nodes,'rows');
         face1p(3) = [face11;face12];

         % Find the faces associated with the second new cell
         [~,~,face11] = intersect(nodes2([1,2]), face_nodes,'rows');
         [~,~,face12] = intersect(nodes2([2,1]), face_nodes,'rows');
         face2p(1) = [face11;face12];

         [~,~,face11] = intersect(nodes2([2,3]), face_nodes,'rows');
         [~,~,face12] = intersect(nodes2([3,2]), face_nodes,'rows');
         face2p(2) = [face11;face12];

         [~,~,face11] = intersect(nodes2([3,1]), face_nodes,'rows');
         [~,~,face12] = intersect(nodes2([1,3]), face_nodes,'rows');
         face2p(3) = [face11;face12];

          
          check_node1 = setdiff( nodes2, nodes1 );
          check_node2 = setdiff( nodes1, nodes2 );
          
          incircl = incircle_pred(x(nodes1(1),1), x(nodes1(1),2),... 
                                  x(nodes1(2),1), x(nodes1(2),2),...
                                  x(nodes1(3),1), x(nodes1(3),2),...
                                  x(check_node1,1), x(check_node1,2) );
         
          %Compare opposite angles of triangle
%           x0 = x(check_node2,:);
%           x1 = x(setdiff(nodes1,check_node2)',:);
%           V1 = x1(1,:) - x0;
%           V2 = x1(2,:) - x0;
%           L1 = sqrt( V1(1).^2 + V1(2).^2 );
%           L2 = sqrt( V2(1).^2 + V2(2).^2 );
%           th1 = acosd( (V1*V2')/(L1*L2) );
%               
%           x0 = x(check_node1,:);
%           x1 = x(setdiff(nodes2,check_node1)',:);
%           V1 = x1(1,:) - x0;
%           V2 = x1(2,:) - x0;
%           L1 = sqrt( V1(1).^2 + V1(2).^2 );
%           L2 = sqrt( V2(1).^2 + V2(2).^2 );
%           th2 = acosd( (V1*V2')/(L1*L2) );

         if incircl || override%|| th1+th2 >180
             shared_nodes = intersect(nodes1, nodes2);
             new_cell1 = [check_node2, shared_nodes(1), check_node1];
             new_cell2 = [check_node1, shared_nodes(2), check_node2];
             
             newface = [check_node1, check_node2];
             face_nodes(face,:) = newface;
             
             new_cell_nodes1 = order_nodes(new_cell1, x);
             new_cell_nodes2 = order_nodes(new_cell2, x);
             
             
             % Compare new cell nodes and face nodes to find the faces
             % associated with the new cells.
             
             % Find the faces associated with the first new cell
             [~,~,face11] = intersect(new_cell_nodes1([1,2]), face_nodes,'rows');
             [~,~,face12] = intersect(new_cell_nodes1([2,1]), face_nodes,'rows');
             face1(1) = [face11;face12];
             
             [~,~,face11] = intersect(new_cell_nodes1([2,3]), face_nodes,'rows');
             [~,~,face12] = intersect(new_cell_nodes1([3,2]), face_nodes,'rows');
             face1(2) = [face11;face12];
             
             [~,~,face11] = intersect(new_cell_nodes1([3,1]), face_nodes,'rows');
             [~,~,face12] = intersect(new_cell_nodes1([1,3]), face_nodes,'rows');
             face1(3) = [face11;face12];
             
             % Find the faces associated with the second new cell
             [~,~,face11] = intersect(new_cell_nodes2([1,2]), face_nodes,'rows');
             [~,~,face12] = intersect(new_cell_nodes2([2,1]), face_nodes,'rows');
             face2(1) = [face11;face12];
             
             [~,~,face11] = intersect(new_cell_nodes2([2,3]), face_nodes,'rows');
             [~,~,face12] = intersect(new_cell_nodes2([3,2]), face_nodes,'rows');
             face2(2) = [face11;face12];
             
             [~,~,face11] = intersect(new_cell_nodes2([3,1]), face_nodes,'rows');
             [~,~,face12] = intersect(new_cell_nodes2([1,3]), face_nodes,'rows');
             face2(3) = [face11;face12];
             
             face_swap1 = setdiff( face1p, face1 );
             face_swap2 = setdiff( face2p, face2 );
             
             I1 = find(cell1==face_cells(face_swap1,:));
             I2 = find(cell2==face_cells(face_swap2,:));
             
             face_cells(face_swap2,I2)=cell1;
             face_cells(face_swap1,I1)=cell2;
             
             ti(cell1,:) = new_cell1;
             ti(cell2,:) = new_cell2;
             
%              [face_nodes2, cell_faces2, face_cells2] = find_neighbors(ti);
             
             
%              figure(2)
%             [face_nodes1, cell_faces1] = find_neighbors(ti([cell1;cell2],:));
%             for j = 1:size(face_nodes1,1)
%                 if (j==1); hold off; end
%                 plot( x(face_nodes1(j,:),1), x(face_nodes1(j,:),2) ,'k-o','LineWidth',2);
%                 hold on
% 
%             end
%              
%              figure(3)
%             [face_nodes2, cell_faces2] = find_neighbors(tiold([cell1;cell2],:));
%             for j = 1:size(face_nodes2,1)
%                 if (j==1); hold off; end
%                 plot( x(face_nodes2(j,:),1), x(face_nodes2(j,:),2) ,'k-o','LineWidth',2);
%                 hold on
% 
%             end 
            
         end
