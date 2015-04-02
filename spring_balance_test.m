function x = spring_balance_test(ti,x,bfn)
%% Spring balance test
% clc
% clear all

% load('triangulation.mat')
[face_nodes, cell_faces, face_cells] = find_neighbors(ti);

figure(1)
plot_cells(ti,x)
axis equal

% x1 = [0,0];
% x2 = [0.5,1];
% x3 = [1,0];
% 
% d = @(x1,x2) sqrt( (x2-x1).^2*(x2'-x1').^2 );

[bndry, IA] = setdiff(face_nodes,bfn,'rows');
bndry_nodes = unique(bfn);

Fx = zeros(length(x),1);
Fy = zeros(length(x),1);
A = zeros(length(x),length(x));

for i = 1:size(x,1)
    if any(i==bndry_nodes)
        A(i,i) = 1;
        Fx(i,1) = x(i,1);
        Fy(i,1) = x(i,2);
        
    else
    
        % Find node neighbors
        clear node_nbor
        cnt = 1;
        for j = 1:size(face_nodes,1)
              I = find( any(face_nodes(j,:) == i ) );
              if length(I) > 0
                  J = find( face_nodes(j,:) ~= i ); 
                  node_nbor(cnt) = face_nodes(j,J);
                  cnt = cnt + 1;
              end
        end    

        A(i,i) = length(node_nbor);
        for j = 1:length(node_nbor)
            A(i,node_nbor(j)) = -1;
        end


    end
    
end



% i = 1;
% Ax = [x1(i)-x2(i), 0, x1(i)-x3(i)
%       x2(i)-x1(i), x2(i)-x3(i), 0
%       0, x3(i)-x2(i), x3(i)-x1(i)];
% i = 2;
% Ay = [x1(i)-x2(i), 0, x1(i)-x3(i)
%       x2(i)-x1(i), x2(i)-x3(i), 0
%       0, x3(i)-x2(i), x3(i)-x1(i)];
% k = [1, 1, 1]';

% Fx = A*x(:,1);
% Fy = A*x(:,2);
% for i = 1:size(face_cells,1)
%    I=find(face_cells(i,:)==0);
%    if length(I)>0
%        
%    end
% end
% for i=1:length(bfn_lbl)
%     if (bfn_lbl(i)>0)
%        A(i,:) = 0;
%        A(i,i) = 1;
%        Fx(i,1) = x(i,1);
%        Fy(i,1) = x(i,2);
%     end
% end

x(:,1) = A\Fx;
x(:,2) = A\Fy;

X = x(:,1);
Y = x(:,2);
figure(2)
plot_cells(ti,x)
title('quiver')
quiver(X,Y,Fx,Fy)
axis equal

