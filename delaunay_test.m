%% Delaunay mesh generation
clear all
clc

qual_max = 1;
r_max = .15;
nmax = 200;

% generate a circular boundary mesh
r = 1;
n = 41;
% th = linspace(0,2*pi,n)';
% x = r*[cos(th), sin(th)];

% xx = linspace(0,1,n);
% dx = xx(2)-xx(1);
% dy = dx*sind(60);
% yy = 0:dy:(n-1)*dy;
% [xx,yy]=meshgrid(xx,yy);
% xx=xx';
% yy=yy';
% for j = 1:size(xx,1)
%    if mod(j,2)==0
%        xx(:,j) = xx(:,j) + dx/2;
%    end
% end
% x(:,1) = reshape(xx,[numel(xx),1]);
% x(:,2) = reshape(yy,[numel(yy),1]);

xx = r*(rand(1,n)*2-1);
yy = r*(rand(1,n)*2-1);
x = [xx',yy'];

%Circular boundary
rad = @(x) sqrt( sum( x.^2,2 ) );
th = @(x) atan(x(:,2)./x(:,1)) + (x(:,1)<0)*pi + (x(:,1)>0 & x(:,2)<0)*2*pi;
dth = @(x1,x2) th(:,x2) - th(:,x1);
bfun(1).fun = @(x1,x2) [rad(x2).*cos( (th(x2) + th(x1))/2 ), rad(x2).*sin( (th(x2) + th(x1))/2  )];

%Airfoil
ynaca0012p = @(x) 0.594689181*(0.298222773*sqrt(x) - 0.127125232.*x - 0.357907906*x.^2 + 0.291984971*x.^3 - 0.105174606*x.^4);
ynaca0012n = @(x) - ynaca0012p(x);
xc = @(x1,x2) (x1+x2)/2;
yfun = @(x1,x2) ynaca0012p( xc(x1(:,1), x2(:,1)) ).*(xc(x1(:,2), x2(:,2))>=0) +ynaca0012n( xc(x1(:,1), x2(:,1)) ).*(xc(x1(:,2), x2(:,2))<0);
bfun(2).fun = @(x1,x2) [xc(x1(:,1), x2(:,1)), yfun(x1,x2) ];

        
x = [];
bfn = [];
bfun_lbl = [];
        
%Initial mesh outer boundary mesh
n1 = 40;
theta = linspace(0,2*pi,n1);
theta = theta(1:end-1);
xbndry = flipud([r*cos(theta'),r*sin(theta')]);
% x(:,2) = bfun(1).fun(x(1:end-1,:), x(2:end,:));
bfn1(:,1) = [1:size(xbndry,1)];
bfn1(:,2) = [2:size(xbndry,1),1];
bfun_lbl1 = ones(size(bfn1,1),1);
% bfun(1).fun = @(x1,x2) [xc(x1(:,1), x2(:,1)), xc(x1(:,2), x2(:,2)) ];
% bfn_faces = [xbndry];

x = [x;xbndry];
bfn = [bfn; bfn1];
bfun_lbl = [bfun_lbl; bfun_lbl1];


% Square Mesh
% clear bfn bfun_lbl
% x1 = linspace(0,r,n1);
% y1 = zeros(size(x1));
% x3 = fliplr(linspace(0,r/2,n1));
% y3 = x3*tand(60);
% x2 = x3+r/2;
% y2 = fliplr(y3);

% x = [x1,x2(2:end),x3(2:end-1)];
% y = [y1,y2(2:end),y3(2:end-1)];
% x = [x',y'];
% bfn(:,1) = [1:size(x,1)-1,1];
% bfn(:,2) = [2:size(x,1),size(x,1)];
% bfun_lbl = ones(size(bfn,1),1);




%Initial airfoil mesh
% n2 = 11;
% theta = linspace(0,pi,n2)';
% xairfoil = 1/2*cos(theta)+1/2;
% yairfoil = ynaca0012p(xairfoil);
% theta = linspace(pi,2*pi,n2)';
% xairfoil2 = 1/2*cos(theta)+1/2;
% xairfoil2 = xairfoil2(2:end);
% yairfoil2 = ynaca0012n(xairfoil2);
% xairfoil = [xairfoil, yairfoil
%             xairfoil2, yairfoil2];
%         
% bfn2(:,1) = [1:size(xairfoil,1)-1]+size(x,1);
% bfn2(:,2) = [2:size(xairfoil,1)]+size(x,1);
% bfun_lbl2 = 2*ones(size(bfn2,1),1);
% 
% x = [x;xairfoil];
% bfn = [bfn; bfn2];
% bfun_lbl = [bfun_lbl; bfun_lbl2];


% 
% xx = linspace(-r,r,n);
% [yy,xx] = meshgrid(xx,xx);
% x(:,1) = reshape(xx,[numel(xx),1]);
% x(:,2) = reshape(yy,[numel(xx),1]);

% load('mesh.mat')
% nodes=x;
% x = [0.8, -0.6
%      -0.75, 0.3
%      -0.3, -0.4
%      0.2,0.1];


plot(x(:,1),x(:,2),'k-o')

hold off
plot(0,0,'.')

% Generate bounding box
xbox = [-r*1.1, -r*1.5
         r*1.1, -r*1.5
         r*1.1,  r*1.5
         -r*1.1, r*1.5];

% and concatinate to current list of nodes
nodes = [xbox; x];

% Current cells
ti = [1, 2, 4
      4, 2, 3];

for i = 1:size(nodes,1)-4
    [ti, ~] = delaunay_add_point(ti,x(i,:),nodes(1:i+4-1,:));
    plot_cells(ti,nodes)
end

ti = remove_corners(ti);

%Remove exterior cells
[face, cell] =  compute_face_data_from_triangulation(ti,nodes);
for i = 1:length(face); fn(i,:) = face(i).nodes; end

cnt = 1;
bfn = bfn+4;
for i = 1:size(bfn,1)
    [~,Ia,Ib] = intersect( [bfn(i,:);fliplr(bfn(i,:))], fn , 'rows');
    if length(Ib>0)
        if (Ia==1)
        cn = face(i).cell_neg;
        if (cn>0)
        Inotkeep(cnt) = cn;
        cnt = cnt + 1;
        end
        
        elseif (Ia==2)
            cn = face(i).cell_plus;
            if (cn>0)
            Inotkeep(cnt) = cp;
            cnt = cnt + 1;
            end
        end
    end
end
Iall = 1:length(ti);
Ikeep = setdiff(Iall,Inotkeep);
% ti = ti(Ikeep,:);

for i = 1:size(ti,1)
    plot_cells(ti(i,:),x)
    hold on
    plot(x(:,1),x(:,2),'r*')
end

figure(1)
plot_cells(ti,x)
axis equal


% [face_nodes, cell_faces, face_cells] = find_neighbors(ti);
% tiold = ti;
% for i = 1:size(face_cells,1)
%       if all(face_cells(i,:)~=0)
%           cell1 = face_cells(i,1);
%           cell2 = face_cells(i,2);
% 
%           [ti, face_nodes, face_cells] = delaunay_face_swap(cell1, cell2, i, ti, face_nodes, face_cells, x);
%             figure(1)
%             plot_cells(ti,x)
%       end
% end




if ~is_delaunay(ti,x)
    fprintf('Mesh is a Delaunay triangulation! Woot!\n')
end

    %% Quality measure refinement
    bfn = bfn-4;
    quality=[];
    for j = 1:size(ti,1)
    xvec = x(ti(j,:),1)';
    yvec = x(ti(j,:),2)';
    [quality(j), r(j), x0(j,:)] = mesh_quality_r_to_d( xvec, yvec );
    end
    
    i = 1;
    max_qual = 2;
    max_r = 1;

    n = size(ti,1);
    while (max_qual>qual_max || max_r > r_max) && n < nmax;


        ismod = 0;
        [max_qual,I] = max(quality);
        [max_r, Ir]  = max(r);

        if (max_qual > qual_max*(1+eps))
            [ti, x, bfn, bfun_lbl, cell_mod, cell_same] = delaunay_add_point(ti,x0(I,:),x,bfn, bfun_lbl, bfun);
%             x = spring_balance_test(ti,x,bfn);
            ismod = 1;
        elseif (max_r > r_max)
            [ti, x, bfn, bfun_lbl, cell_mod, cell_same] = delaunay_add_point(ti,x0(Ir,:),x,bfn, bfun_lbl, bfun);
%             x = spring_balance_test(ti,x,bfn);
            ismod = 1;    
        end

        if ismod == 1;
            % Update quality and x0 for modified cells
            if (length(cell_same)>0)
            quality = quality(cell_same);
            x0 = x0(cell_same,:);
            r = r(cell_same);
            else
               quality = [];
               x0=[];
               r=[];
                
            end
            for j = cell_mod;
              xvec = x(ti(j,:),1)';
              yvec = x(ti(j,:),2)';
              [quality(j), r(j), x0(j,:)] = mesh_quality_r_to_d( xvec, yvec );
            end
           
        end

        [face_nodes, cell_faces, face_cells] = find_neighbors(ti);
        
        
% I = find(bfun_lbl==2);
% bfn2 = bfn(I,:);
% cell_remove = [];
% for j = 1:length(I)
%    bfntest = [];
%    bfntest(:,2) = I(j+2:end);
%    bfntest(:,1) = I(j);
%    [~,~,fremove] = intersect(bfntest,face_nodes,'rows');
% %    cell_remove = [cell_remove; face_cells(fremove,:)];
%    cell_remove = face_cells(fremove,:);
%    
% %    if ~isempty(fremove)
% %    cell_remove = unique(cell_remove);
% %    cell_keep =  setdiff([1:size(ti,1)]',cell_remove);
% %    ti = ti(cell_keep,:);
% %    [face_nodes, cell_faces, face_cells] = find_neighbors(ti);
% %    figure(1)
% %    plot_cells(ti,x)
% %    axis([-0.2,1.2,-0.5,0.5])
% %    end
% end
% cell_remove = unique(cell_remove);
% cell_keep =  setdiff([1:size(ti,1)]',cell_remove);
% ti = ti(cell_keep,:);
%         
        
        
        n = size(ti,1);

        if (mod(n,10)==0)
        figure(2)
        plot_cells(ti,x)
        axis equal
%         axis([-0.5,1.5,-0.5,0.5])
        end
        disp([max_qual , max_r , n ])
    
        
figure(2)
plot_cells(ti,x)
axis equal
disp([max_qual>qual_max , max_r > r_max, n < nmax])
    end

save('circle_mesh.mat','ti','x')
    
nodes = [xbox; x];
ti = [1, 2, 4
      4, 2, 3];
for i = 1:size(nodes,1)-4
    [ti, ~] = delaunay_add_point(ti,x(i,:),nodes(1:i+4-1,:));
end
ti = remove_corners(ti);
plot_cells(ti,x)
axis equal

% Remove the corners
% cnt = 1;
% for i = 1:size(ti,1)
%    if any(ti(i,:)==1) || ...
%       any(ti(i,:)==2) || ...
%       any(ti(i,:)==3) || ...
%       any(ti(i,:)==4) 
%   
%    else
%        ticlean(cnt,:) = ti(i,:);
%        cnt = cnt + 1;
%    end
% end
% ti = ticlean;
% 
%     figure(2)
%     [face_nodes, cell_faces] = find_neighbors(ti);
% %     figure(1)
%     for j = 1:size(face_nodes,1)
%         if (j==1); hold off; end
%     %     cell = [nodes(ti(i,[1;2;3;1]),1),nodes(ti(i,[1;2;3;1]),2)];
%         plot( nodes(face_nodes(j,:),1), nodes(face_nodes(j,:),2) ,'k-o','LineWidth',2);
%     %     plot(cell(:,1),cell(:,2),'k-');
%         hold on
% 
%     end
  
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    


%% Delaunay Kernel
% for i = 1:size(nodes,1)-4
%     incir = incircle_pred(nodes(ti(:,1),1), nodes(ti(:,1),2),...
%                           nodes(ti(:,2),1), nodes(ti(:,2),2),...
%                           nodes(ti(:,3),1), nodes(ti(:,3),2),...
%                           x(i,1), x(i,2));
%                       
% %     figure(1)
% %     hold off
% %     plot(nodes(i+4,1),nodes(i+4,2),'g*')
% %     hold on
% %     
% %     for j = 1:size(face_nodes,1)
% % %         if (i==1); hold off; end
% %     %     cell = [nodes(ti(i,[1;2;3;1]),1),nodes(ti(i,[1;2;3;1]),2)];
% %             fn = face_nodes(j,:);
% %             plot( nodes(fn,1), nodes(fn,2) ,'k-o','LineWidth',2);
% %             axis([-2*r,2*r,-2*r,2*r])
% %     %     plot(cell(:,1),cell(:,2),'k-');
% % 
% % 
% %     end
% %     
% %     for j = 1:size(ti,1)
% %         
% %         [r,x0,y0] = find_circle(nodes(ti(j,1),1), nodes(ti(j,1),2),...
% %                                 nodes(ti(j,2),1), nodes(ti(j,2),2),...
% %                                 nodes(ti(j,3),1), nodes(ti(j,3),2));
% %         if (incir(j)==1)
% %             color='g--';
% %         else
% %             color='r--';
% %         end
% %         xcir = r*cos(theta)+x0;
% %         ycir = r*sin(theta)+y0;
% %         plot(xcir,ycir,color);
% % %         axis equal
% %     end
% 
%                       
%                       
%                       
%                       
%                       
%                       
%     I = find(incir==1);
%     Inot = find(incir~=1);
%     cf = cell_faces(I,:);
%     [fu, IA, IC] = unique(cf);
%     nface_to_remove = length(IC)-length(IA);
%     cnt = 1;
%     cnt2 = 1;
% %     if (nface_to_remove==0)
% % figure(2)
% %     for j = 1:size(fu,1)
% %         if (i==1); hold off; end
% %     %     cell = [nodes(ti(i,[1;2;3;1]),1),nodes(ti(i,[1;2;3;1]),2)];
% % %         if any(j==fdup)
% %             fn = face_nodes(fu(j),:);
% %             plot( nodes(fn,1), nodes(fn,2) ,'g-o','LineWidth',2);
% %             hold on
% % %         end
% %     %     plot(cell(:,1),cell(:,2),'k-');
% % 
% % 
% %     end
%         
% %         figure(1)
% % %     else
%         fdup=[];
%         uniq_faces=[];
%         for j = 1:length(fu)
%           II = find(fu(j)==cf);
%           if length(II)>1
%               fdup(cnt,1) = fu(j);
%               cnt = cnt + 1;
%           else
%               uniq_faces(cnt2,1) = fu(j);
%               cnt2 = cnt2 + 1;
%           end
%         end
%     
% %     end
%     
%     
%     for j = 1:length(fdup)
% %         if (i==1); hold off; end
%     %     cell = [nodes(ti(i,[1;2;3;1]),1),nodes(ti(i,[1;2;3;1]),2)];
% %         if any(fu(j)==fdup)
%             fn = face_nodes(fdup(j),:);
%             plot( nodes(fn,1), nodes(fn,2) ,'g-o','LineWidth',2);
% %             
% %         end
%     %     plot(cell(:,1),cell(:,2),'k-');
% 
% 
%     end
% %     axis([-1.2*r,1.2*r,-1.6*r,1.6*r])
% %     fdup = cf(IA(end-nface_to_remove+1:end));
% %     uniq_faces = IA(1:end-nface_to_remove);
%     
%     tip1=[face_nodes(uniq_faces,:),(i+4)*ones(length(uniq_faces),1)];
%     tip1 = [ti(Inot,:);tip1];
%     ti = order_nodes(tip1, nodes); 
%     
%     [face_nodes, cell_faces] = find_neighbors(ti);
% %     figure(1)
% %     for j = 1:size(face_nodes,1)
% % %         if (j==1); hold off; end
% %     %     cell = [nodes(ti(i,[1;2;3;1]),1),nodes(ti(i,[1;2;3;1]),2)];
% %         plot( nodes(face_nodes(j,:),1), nodes(face_nodes(j,:),2) ,'k-o','LineWidth',2);
% %     %     plot(cell(:,1),cell(:,2),'k-');
% % %         hold on
% % 
% %     end
% %     % Plot the current triangulation
% 
% % plot(nodes(:,1),nodes(:,2),'ko')
% 
% 
% 
% 
%  end


    
    
    
    
    
    
%% FACE SWAP ALGORITHM TESTS
% 
% [face_nodes, cell_faces, face_cells] = find_neighbors(ti);
% 
% tiold = ti;
% for i = 1:size(face_cells,1)
%       if all(face_cells(i,:)~=0)
%           cell1 = face_cells(i,1);
%           cell2 = face_cells(i,2);
% %           common_face = intersect(cell1,cell2);
%           nodes1 = ti(cell1,:);
%           nodes2 = ti(cell2,:);
% 
%           % Find the faces associated with the first new cell
%          [~,~,face11] = intersect(nodes1([1,2]), face_nodes,'rows');
%          [~,~,face12] = intersect(nodes1([2,1]), face_nodes,'rows');
%          face1p(1) = [face11;face12];
% 
%          [~,~,face11] = intersect(nodes1([2,3]), face_nodes,'rows');
%          [~,~,face12] = intersect(nodes1([3,2]), face_nodes,'rows');
%          face1p(2) = [face11;face12];
% 
%          [~,~,face11] = intersect(nodes1([3,1]), face_nodes,'rows');
%          [~,~,face12] = intersect(nodes1([1,3]), face_nodes,'rows');
%          face1p(3) = [face11;face12];
% 
%          % Find the faces associated with the second new cell
%          [~,~,face11] = intersect(nodes2([1,2]), face_nodes,'rows');
%          [~,~,face12] = intersect(nodes2([2,1]), face_nodes,'rows');
%          face2p(1) = [face11;face12];
% 
%          [~,~,face11] = intersect(nodes2([2,3]), face_nodes,'rows');
%          [~,~,face12] = intersect(nodes2([3,2]), face_nodes,'rows');
%          face2p(2) = [face11;face12];
% 
%          [~,~,face11] = intersect(nodes2([3,1]), face_nodes,'rows');
%          [~,~,face12] = intersect(nodes2([1,3]), face_nodes,'rows');
%          face2p(3) = [face11;face12];
%           
% %           face1p = cell_faces(cell1,:);
% %           face2p = cell_faces(cell2,:);
%           
%           check_node1 = setdiff( nodes2, nodes1 );
%           check_node2 = setdiff( nodes1, nodes2 );
%           
%           incircl = incircle_pred(x(nodes1(1),1), x(nodes1(1),2),... 
%                                   x(nodes1(2),1), x(nodes1(2),2),...
%                                   x(nodes1(3),1), x(nodes1(3),2),...
%                                   x(check_node1,1), x(check_node1,2) );
%          
%           %Compare opposite angles of triangle
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
%  disp([incircl,th1, th2, th1+th2])
%  
%          if incircl || th1+th2 >180
%              shared_nodes = intersect(nodes1, nodes2);
%              new_cell1 = [check_node2, shared_nodes(1), check_node1];
%              new_cell2 = [check_node1, shared_nodes(2), check_node2];
%              
%              newface = [check_node1, check_node2];
%              face_nodes(i,:) = newface;
%              
%              new_cell_nodes1 = order_nodes(new_cell1, x);
%              new_cell_nodes2 = order_nodes(new_cell2, x);
%              
%              
%              % Compare new cell nodes and face nodes to find the faces
%              % associated with the new cells.
%              
%              % Find the faces associated with the first new cell
%              [~,~,face11] = intersect(new_cell_nodes1([1,2]), face_nodes,'rows');
%              [~,~,face12] = intersect(new_cell_nodes1([2,1]), face_nodes,'rows');
%              face1(1) = [face11;face12];
%              
%              [~,~,face11] = intersect(new_cell_nodes1([2,3]), face_nodes,'rows');
%              [~,~,face12] = intersect(new_cell_nodes1([3,2]), face_nodes,'rows');
%              face1(2) = [face11;face12];
%              
%              [~,~,face11] = intersect(new_cell_nodes1([3,1]), face_nodes,'rows');
%              [~,~,face12] = intersect(new_cell_nodes1([1,3]), face_nodes,'rows');
%              face1(3) = [face11;face12];
%              
%              % Find the faces associated with the second new cell
%              [~,~,face11] = intersect(new_cell_nodes2([1,2]), face_nodes,'rows');
%              [~,~,face12] = intersect(new_cell_nodes2([2,1]), face_nodes,'rows');
%              face2(1) = [face11;face12];
%              
%              [~,~,face11] = intersect(new_cell_nodes2([2,3]), face_nodes,'rows');
%              [~,~,face12] = intersect(new_cell_nodes2([3,2]), face_nodes,'rows');
%              face2(2) = [face11;face12];
%              
%              [~,~,face11] = intersect(new_cell_nodes2([3,1]), face_nodes,'rows');
%              [~,~,face12] = intersect(new_cell_nodes2([1,3]), face_nodes,'rows');
%              face2(3) = [face11;face12];
%              
%              face_swap1 = setdiff( face1p, face1 );
%              face_swap2 = setdiff( face2p, face2 );
%              
%              I1 = find(cell1==face_cells(face_swap1,:));
%              I2 = find(cell2==face_cells(face_swap2,:));
%              
%              face_cells(face_swap2,I2)=cell1;
%              face_cells(face_swap1,I1)=cell2;
%              
%              ti(cell1,:) = new_cell1;
%              ti(cell2,:) = new_cell2;
%              
% %              [face_nodes2, cell_faces2, face_cells2] = find_neighbors(ti);
%              
%              
%              figure(2)
%             [face_nodes1, cell_faces1] = find_neighbors(ti([cell1;cell2],:));
%         %     figure(1)
%             for j = 1:size(face_nodes1,1)
%                 if (j==1); hold off; end
%             %     cell = [nodes(ti(i,[1;2;3;1]),1),nodes(ti(i,[1;2;3;1]),2)];
%                 plot( x(face_nodes1(j,:),1), x(face_nodes1(j,:),2) ,'k-o','LineWidth',2);
%             %     plot(cell(:,1),cell(:,2),'k-');
%                 hold on
% 
%             end
%              
%              figure(3)
%             [face_nodes2, cell_faces2] = find_neighbors(tiold([cell1;cell2],:));
%         %     figure(1)
%             for j = 1:size(face_nodes2,1)
%                 if (j==1); hold off; end
%             %     cell = [nodes(ti(i,[1;2;3;1]),1),nodes(ti(i,[1;2;3;1]),2)];
%                 plot( x(face_nodes2(j,:),1), x(face_nodes2(j,:),2) ,'k-o','LineWidth',2);
%             %     plot(cell(:,1),cell(:,2),'k-');
%                 hold on
% 
%             end 
%             
%          end
%          
%       end
% end

    
    
    
    
    
    
    

