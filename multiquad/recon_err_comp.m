clc
% clear all

scale = logspace(-5,5,11);
% scale=1;
order = [1,2,3,4,5,6];

for o=1:length(order)
for i = 1:length(scale)
   fprintf('%1.0f %12.4e',order(o),scale(i));
   error(o,i)=reconst_err_test(order(o),scale(i));
   fprintf('%12.4e\n',error(o,i));
end
fprintf('\n');
end
   

loglog(scale,error)

fid = fopen('fit_error_fitn_onlyroundoff.dat','w');
fprintf(fid,'Variables="scale""error"\n');

for i = 1:length(order)
   fprintf(fid,['zone T="order=',num2str(order(i)),'"\n']);
   fprintf(fid,'i=%2.0f\n',length(error));
   
   for j=1:length(error)
      fprintf(fid,'%23.15e %23.15e \n',scale(j),error(i,j));       
   end
       
end