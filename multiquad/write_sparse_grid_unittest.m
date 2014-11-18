%% write NqD results to file for unit testing fortran version

order = 20;
d = 3;
method = 'CC';

[ x, w ] = sparse_grid( d, order, method );

fid = fopen(['sparse_grid_order',num2str(order),'_d',num2str(d),'_m',method,'.dat'],'w');
for i = 1:size(x,1)
    for j = 1:d
        fprintf(fid,'%23.15f', x(i,j) );
    end
    fprintf(fid,'\n');
end
for i = 1:size(w,1)
        fprintf(fid,'%23.15f\n', w(i) );
end
fclose(fid);