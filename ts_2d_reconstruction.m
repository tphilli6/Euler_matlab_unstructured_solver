function [Ainv, cell, A] = ts_2d_reconstruction(cell, n, kexact_order)
%% Notes
% Taylor series reconstruction

%F(x+dx,y+dy) = F(i) 
%             + dx*dfdx + dy*dfdy 
%             + dx^2/2*d2fdx2 + dy^2/2*d2fdy2 + dx*dy/2*d2fdxdy
%             + dx^3/6*d3fdx3 + dx^2*dy/6*d3f/dx^2dy + dx*dy^2/6*d3f/dxdy^2 + dy^3/6*d3f/dy^3
%
% Fbar = 1/(dx*dy) int(dx,dy) F(i)
%      = 1/vol*[F(i)*dx*dy
%                   + dx^2*dy/2*dfdx + dx*dy^2/2*dfdy
%                   + dx^3*dy/6*d2fdx2 + dy^3*dx/6*d2fdy2 + dx^2*dy^2/8*d2fdxdy
%                   + dx^4/24*dy*d3fdx3 + dx^3*dy^2/36*d3f/dx^2dy + dx^2*dy^3/36*d3f/dxdy^2 + dy^4/24*d3f/dy^3
%                   + ... ]
dim = 2;
% kexact_order = 1;
quad_range_option_for_tri_transform = 2;

cnt = 1;
for j = 0:kexact_order
    for i = 0:kexact_order-j
      p(1,cnt) = i;
      p(2,cnt) = j;
      cnt = cnt + 1;
    end
end
% cell.nunknowns = size(p,2);


% Creates a row vector Ai for a given x
Ai_row = @(x) [x(1).^p(1,:).*x(2).^p(2,:)];


if (kexact_order == 0)
  xcc_quad = [0.5,0.5];
  wcc = [1];
else
  [xcc_quad, wcc] = sparse_grid(dim, dim*kexact_order, 'CC');
end

% for n = 1:cell.ncells
    
    %Loop over all the cells in the stencil
    cnt = 1;
    for nn = cell.stencil(n).cells
        
      % Transform the quadrature points to either a quad or a triangle
      nnodes = cell.nodes(nn,1);
      if (nnodes == 3)
        % triangle quadrature
        % Transform the quadrature from a quad to a triangle 
        for i = 1:length(wcc)
          xcc(i,:) = transform_quad_to_triangle(xcc_quad(i,:),...
                                         quad_range_option_for_tri_transform);
        end

      elseif (nnodes == 4)
        % quad quadrature
        xcc = xcc_quad;
      end

      % Evalutate the quadrature points using the current cell map to find
      % physical location
      for nq = 1:length(wcc)
        xq(1,1) = cell.map(nn).x(xcc(nq,:));
        xq(1,2) = cell.map(nn).y(xcc(nq,:));

        % Evaluate the function at the quadrature point with the current
        % cell center xc(n) as the reference
        A(nq,:) = Ai_row( xq - cell.xc(n,1:2) )*wcc(nq);
      end

      Ai(cnt,:) = sum(A);%.*cell.volume(n);
      cnt = cnt + 1;
%       cell.reconstruction(nn).Ainv = Ai;

%       cell.reconstruction(nn).Axc = Ai_row(cell.xc(nn,1:2));
    end
    
         % Store the evaluation of the averaging operator for the current
         % cell
         cell.reconstruction(n).Ai = Ai(1,:);
         cell.reconstruction(n).A  = Ai;
         A = Ai;
          % Constrain the reconstruction
          Ai = Ai(2:end,2:end);
          Ainv = (Ai'*Ai).^(-1)*Ai';
%           cell.reconstruction(n).Ainv = Ainv;
    
% end
    