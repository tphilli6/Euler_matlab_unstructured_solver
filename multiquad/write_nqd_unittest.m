%% write NqD results to file for unit testing fortran version

q = 20;
d = 4;

ii = NqD(q,d);

fid = fopen(['Nqd_q',num2str(q),'_d',num2str(d),'.dat'],'w');
for i = 1:size(ii,1)
    for j = 1:d
        fprintf(fid,'%18.0f', ii(i,j) );
    end
    fprintf(fid,'\n');
end
fclose(fid);