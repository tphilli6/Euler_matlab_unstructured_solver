function [tip1, nodes_out, bfn, bfun_lbl, cell_mod, cell_same] = delaunay_add_point(ti,xadd,nodes,bfn, bfun_lbl, bfun)


[face_nodes, cell_faces] = find_neighbors(ti);

fremove = 0;
if (nargin>=6)
    [xadd, bfn, bfun_lbl,fu_bndry] = chew_algorithm(xadd, bfn, nodes, bfun_lbl, bfun);
    fu_bndry = sort(fu_bndry);
    [~, fremove, ~] = intersect(face_nodes, fu_bndry, 'rows');
    if isempty(fremove); fremove = 0; end
end

  
%% Delaunay Kernel
incir = incircle_pred(nodes(ti(:,1),1), nodes(ti(:,1),2),...
                      nodes(ti(:,2),1), nodes(ti(:,2),2),...
                      nodes(ti(:,3),1), nodes(ti(:,3),2),...
                      xadd(1), xadd(2));

%% Plot circles for each triangle (red if xadd is not in circle, green otherwise) 
% figure(2)
% plot_cells(ti,nodes)
% theta = linspace(0,2*pi,100);
% plot(xadd(1),xadd(2),'ko')
% hold on
% for j = 1:size(ti,1)
%         [r,x0,y0] = find_circle(nodes(ti(j,1),1), nodes(ti(j,1),2),...
%                                 nodes(ti(j,2),1), nodes(ti(j,2),2),...
%                                 nodes(ti(j,3),1), nodes(ti(j,3),2));
%         if (incir(j)==1)
%             color='g--';
%         else
%             color='r--';
%         end
%         xcir = r*cos(theta)+x0;
%         ycir = r*sin(theta)+y0;
%         plot(xcir,ycir,color);
% %         axis([-2,2,-2,2])
% end               

I = find(incir==1);
Inot = find(incir~=1);

cf = cell_faces(I,:);
[fu, IA, IC] = unique(cf);

cnt = 1;
cnt2 = 1;
fdup=[];
uniq_faces=[];
for j = 1:length(fu)
  II = find(fu(j)==cf);
  if length(II)>1 || fu(j)==fremove
      fdup(cnt,1) = fu(j);
      cnt = cnt + 1;
  else
      uniq_faces(cnt2,1) = fu(j);
      cnt2 = cnt2 + 1;
  end
end

%% Plot interior faces
% for j = 1:length(fdup)
%     fn = face_nodes(fdup(j),:);
%     plot( nodes(fn,1), nodes(fn,2) ,'g-o','LineWidth',2);
% end
% axis([-2,2,-2,2])






nodes_out = [nodes;xadd];
n = size(nodes,1);

tip1=[face_nodes(uniq_faces,:),(n+1)*ones(length(uniq_faces),1)];


cell_mod = [1:size(tip1,1)] + length(Inot);
cell_same = Inot;

tip1 = order_nodes(tip1, nodes_out);
tip1 = [ti(Inot,:);tip1];












    