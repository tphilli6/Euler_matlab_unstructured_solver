function data_out = vertex_to_cell_average(data,cell,cell_cv)


if nargin == 3
    
    for i = 1:length(cell.volume)
        nn = cell.nodes(i,2:end);
        I = find(nn~=0);
        nn = nn(I);
        vol = cell_cv.volume(nn);
        for j = 1:size(data,2)
            data_out(i,j) = vol*data(nn,j)/sum(vol);
        end
    end
    
    
else

  for i = 1:length(cell.volume)
      nn = cell.nodes(i,2:end);
      I = find(nn~=0);
      nn = nn(I);
%       vol = cell.volume(nn);
      data_out(i,:) = mean( data(nn,:) );
  end
  
  
end