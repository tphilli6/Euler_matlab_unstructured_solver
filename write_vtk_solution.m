function write_vtk_solution(vertex, cell, data, soln_name, file_name, attribute, label)


fid = fopen(file_name,attribute);

if (attribute=='w')
  write_vtk(vertex, cell, fid);
  fprintf(fid,'CELL_DATA %8.0f\n',size(data(:,1),1) );
end


for i = 1:size(data,2)
  write_vtk_data(data(:,i), [soln_name, '-', label{i}], fid);
end 

fclose(fid);
