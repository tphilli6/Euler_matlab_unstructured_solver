function [moment] = compute_reconstruction_moments(vertex, cell, face, p )

% % Test input data
% nodes=6;
% kexact_order = 2;
% % x=linspace(0,1,nodes);
% [~,~,~,vertex, face, cell] = generate_mesh(5, 5, 1, 1, 1);
% % [y,x] = meshgrid(x);
% % [vertex, cell, face] = compute_grid_derived_data(x,y,1);
% % xc_m(:,:,1) = (x(1:end-1,1:end-1)+x(1:end-1,2:end)+x(2:end,1:end-1)+x(2:end,2:end))/4;
% % xc_m(:,:,2) = (y(1:end-1,1:end-1)+y(1:end-1,2:end)+y(2:end,1:end-1)+y(2:end,2:end))/4;
% % xc = reshape(xc_m,[numel(xc_m(:,:,1)),2]);
% for i = 1:length(cell.volume)
%     nn = cell.nodes(i,:);
% xc(i,:) = mean( vertex(cell.nodes(i,2:nn+1),:) );
% end
% 
% % Generate taylor series polynomials
% cnt = 1;
% for n = 0:kexact_order
% for m = 0:n
%   p(2,cnt) = m;
%   p(1,cnt) = n-m;
%   ts_coef(1,cnt) = 1/factorial(n)*factorial(n)/(factorial(m)*factorial(n-m));
%   cnt = cnt + 1;
% end
% end
nterms = size(p,2);
cell.reconstruction_param.px = p(1,:);
cell.reconstruction_param.py = p(2,:);


cell.nunknowns = size(p,2);
[xq, wq] = gauss_patterson( max(max(p')) );
% [xq, wq] = gauss_patterson( 2 );

for nn = 1:cell.ncells
  nnodes = cell.nodes(nn,1);
  xcc(1) = mean( vertex(cell.nodes(nn,2:nnodes+1),1) );
  xcc(2) = mean( vertex(cell.nodes(nn,2:nnodes+1),2) );
%   xcc = [0,0];
  
  I = find(cell.faces(nn,:)~=0);
  nf = length(I);
  f = cell.faces(nn,1:nf);
  
%     close all
    for j = 1:nterms
     ix = 0;
     iy = 0;

     int_x = 0;
     int_y = 0;

      for i = 1:nf
          fnormal = face(f(i)).normal;
          fx(:,1) = vertex( face(f(i)).nodes', 1);
          fx(:,2) = vertex( face(f(i)).nodes', 2);
          if (face(f(i)).cell_neg == nn); fnormal = -fnormal; end
          
         for n = 1:length(wq)
            x = fx(1,1) + (fx(2,1)-fx(1,1))*xq(n);
            y = fx(1,2) + (fx(2,2)-fx(1,2))*xq(n);
            
            int_x = int_x + ((x-xcc(1)).^(p(1,j)+1).*(y-xcc(2)).^p(2,j))*wq(n).*fnormal(1)*face(f(i)).area;
            int_y = int_y + ((x-xcc(1)).^(p(1,j)).*(y-xcc(2)).^(p(2,j)+1))*wq(n).*fnormal(2)*face(f(i)).area;

         end
      
%     plot(fx(:,1),fx(:,2),'k-o')
%     hold on
%     fxc = mean(fx,1);
%     quiver(fxc(1),fxc(2),fnormal(1),fnormal(2));
%     axis equal
    
         
      end
      
      int_x = int_x/( (p(1,j)+1)*cell.volume(nn) );
      int_y = int_y/( (p(2,j)+1)*cell.volume(nn) );
      cv_moment(nn,j) = int_x;% + int_y;
      cv_momenty(nn,j) = int_y;
      moment(nn).cv(  p(1,j)+1, p(2,j)+1 ) = cv_moment(nn,j);

      
    end
    
%     cv_moment(nn,:)*[1;0;0;1;0;1]
%     ftest = @(x,y) 1 + (x-xcc(1)).^2 + (y-xcc(2)).^2;
%     ftest(cell.xc(nn,1),cell.xc(nn,2))
 
  end

cell.reconstruction_param.moment = moment;


% % Test routines
% fit_type='lsq';
% build_kexact_stencil
% % cell.stencil(1,:) = [1:length(cell.volume)];
% % imax = floor((nodes-1)/2);
% for n = 1:length(cell.volume);
% cell.lhs(n).A = compute_reconstruction_lhs(cell.stencil(n).cells,n, moment, p, xc);
% % for j = 1:size(cell.stencil,2)
% %     xc(j,:) = cell.xc(cell.stencil(1,j),:);
% % end
% 
% xcc = cell.xc(n,:);
% ftest = @(x,y) 1 + (x-xcc(1)).^2 + 2*(y-xcc(2)).^2;% + (x-xcc(1)).^2 +(x-xcc(1)).*(y-xcc(2)) +  (y-xcc(2)).^2 + (x-xcc(1)).^3 + (y-xcc(2)).^3;
% 
% xtest = xc( cell.stencil(n).cells,:);
% u = ftest(xtest(:,1), xtest(:,2));
% % ucell = u(n);
% % u = u([1:imax*(nodes-1)+imax,imax*(nodes-1)+imax+2:end])-ucell;
% 
% % crow = cv_moment(5,2:end);
% % urow = u(5);
% % u = u([1:4,6:end]);
% % cv_moment = cv_moment([1:4,6:end],[1:4,6:end]);
% % for i = 1:size(cv_moment,1)
% %     cv_moment(i,:) = cv_moment(i,:) - crow;
% %     u(i) = u(i) - urow;
% % end
% %
% a = cell.lhs(n).A\u;
% 
% % coef = [ucell;a]'; 
% err(n) = max(cell.lhs(n).A*a - u);
% end
% 
% plot_cells(cell.nodes(:,2:end),vertex);
% hold on
% plot3(xc(:,1),xc(:,2),err,'ro');
% hold off
% 
