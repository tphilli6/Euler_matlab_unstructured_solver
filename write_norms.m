function write_norms(data, file_base)
% Write the L1, L2, and Linf norms of the columns in array data
% Inputs:
%         data : array of data sorted by columns
%    file_base : base_file name

l1norms = sum( abs(data) ) / size(data,1);
l2norms = sqrt( sum( data.^2 ) / size(data,1) );
linfnorms = max( abs(data) );

fid = fopen([file_base,'-l1.dat'], 'w');
fprintf(fid,'%23.15e',l1norms); fprintf(fid,'\n'); fclose(fid);

fid = fopen([file_base,'-l2.dat'], 'w');
fprintf(fid,'%23.15e',l2norms); fprintf(fid,'\n'); fclose(fid);

fid = fopen([file_base,'-linf.dat'], 'w');
fprintf(fid,'%23.15e',linfnorms); fprintf(fid,'\n'); fclose(fid);
