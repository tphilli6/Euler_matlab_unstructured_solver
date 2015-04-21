function [moment] = compute_reconstruction_moments(vertex, cell, face, p )
% kexact_order = 3;
% x=linspace(0,1,4);
% [y,x] = meshgrid(x);
% [vertex, cell, face] = compute_grid_derived_data(x,y,0);



% Generate taylor series polynomials
% cnt = 1;
% for n = 0:kexact_order
% for m = 0:n
%   p(1,cnt) = m;
%   p(2,cnt) = n-m;
%   ts_coef(1,cnt) = 1/factorial(n)*factorial(n)/(factorial(m)*factorial(n-m));
%   cnt = cnt + 1;
% end
% end
nterms = size(p,2);

cell.nunknowns = size(p,2);
[xq, wq] = curtis_clenshaw( max(max(p)) );



for nn = 1:cell.ncells
  nnodes = cell.nodes(nn,1);
  xcc(1) = mean( vertex(cell.nodes(nn,2:nnodes+1),1) );
  xcc(2) = mean( vertex(cell.nodes(nn,2:nnodes+1),2) );

  nf = length(cell.faces(nn,:));
  f = cell.faces(nn,1:nf);
  

      
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

         end




      end
      
      int_x = int_x/( (p(1,j)+1)*cell.volume(nn) );
      int_y = int_y/( (p(2,j)+1)*cell.volume(nn) );
      cv_moment(nn,j) = int_x;% + int_y;
      moment(nn).cv(  p(1,j)+1, p(2,j)+1 ) = cv_moment(nn,j);
%       cv_momenty(nn,j) = int_y;
      
    end
  
  
%     w(nn,1) = 1/sqrt(  (cell.xc(nn,1)-xcc(1)).^2 + (cell.xc(nn,2)-xcc(2)).^2 );
%     if sqrt(  (cell.xc(nn,1)-xcc(1)).^2 + (cell.xc(nn,2)-xcc(2)).^2 )==0; w(nn,1) = 1; end
%     cv_moment(nn,j) = cv_moment(nn,j)*w(nn,1);
    
end

% cell.moment = cv_moment;

% cell.stencil(1,:) = [1:9];


for i = 1:size(cell.stencil,1)
    xi = cell.xc(5,:);
    
    for j = 1:size(cell.stencil,2)
        xj = cell.xc(cell.stencil(i,j),:);

        for nn = 1:size(p,2)
            m = p(2,nn);
            n = p(1,nn);

            s = 0;
            for l = 0:m
                for k = 0:n
                    cv = moment( cell.stencil(i,j) ).cv(  n-k+1, m-l+1 );

                    s = s + factorial(m)/(factorial(l)*factorial(m-l))...
                    *factorial(n)/(factorial(k)*factorial(n-k))...
                    *(xj(1)-xi(1)).^k.*(xj(2)-xi(2)).^l.*cv;
                end
            end
            
        cell.lhs(i).A(j,nn)  = s;
        
        end
    end
end







% for j = 1:size(cell.stencil,2)
%     xc(j,:) = cell.xc(cell.stencil(1,j),:);
% end
%
% xcc = cell.xc(5,:);
% ftest = @(x,y) 1 + (x-cell.xc(5,1)) + (y-cell.xc(5,2)) + (x-xcc(1)).^2 +(x-xcc(1)).*(y-xcc(2)) +  (y-xcc(2)).^2 + (x-xcc(1)).^3 + (y-xcc(2)).^3;
% 
% u = ftest(xc(:,1), xc(:,2));%.*w;
% % crow = cv_moment(5,2:end);
% % urow = u(5);
% % u = u([1:4,6:end]);
% % cv_moment = cv_moment([1:4,6:end],[1:4,6:end]);
% % for i = 1:size(cv_moment,1)
% %     cv_moment(i,:) = cv_moment(i,:) - crow;
% %     u(i) = u(i) - urow;
% % end
% 
% a = cell.lhs(1).A\u