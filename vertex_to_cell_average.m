function data_out = vertex_to_cell_average(data,cell)

  for i = 1:length(cell.volume)
      data_out(i,:) = mean( data( cell.nodes(i,2:end)',:) );
  end