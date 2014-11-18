function [funcx, funcy, Ainv] = grid_mapping(x, Ainv)

%HARDWIRE: linear reconstruction


% The reconstruction is computed for a quad cell ranging from [0,1]


  if (size(x,1) == 4)
    % coef for p = a + bx + cy + dxy
    coef=@(x) [1, x(1), x(2), x(1).*x(2) ];
    xi = [0, 0
          0, 1
          1, 0
          1, 1];
  
    nunknowns = 4;
  
  elseif (size(x,1) == 3)
    coef=@(x) [1, x(1), x(2)];
    xi = [0, 0
          0, 1
          1, 0];
          
    nunknowns = 3;
  
  end
    
    for i = 1:size(x,1)
      A(i,:) = coef(xi(i,:));
    end

if (nargin == 1);
    Ainv = A^-1;
end


  Ainv;
  bx = x(:,1);
  coefx = Ainv*bx;

  by = x(:,2);
  coefy = Ainv*by;

  funcx = @(xi) coef(xi)*coefx;
  funcy = @(xi) coef(xi)*coefy;
