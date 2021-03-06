%setup mapping
% This is a script to add a mapping function for each cell similar to
%cell.map.funcx(xi)
%cell.map.funcy(xi)

%HARDWIRE: linear mapping for each cell

[~, ~ , Ainv_quad] = grid_mapping([0,0; 1,0; 1,1; 0,1]);
[~, ~ , Ainv_tri] = grid_mapping([0,0; 1,0; 0,1]);

for nn = 1:cell.ncells
  nnodes = cell.nodes(nn,1);
 
  if (nnodes == 3 || nnodes == 4)
      x = vertex( cell.nodes(nn,2:nnodes+1),: );
      if (nnodes == 3)
        Ainv = Ainv_tri;
      elseif (nnodes == 4)
        Ainv = Ainv_quad;
      end

      [cell.map(nn).x, cell.map(nn).y, ~] = grid_mapping(x, Ainv);

  end
end
