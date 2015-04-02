function [coef] = pascals_triangle(row)


c = zeros(row);

c(1,1) = 1;
c(2,1:2) = 1;

for i = 3:row
   c(i,1) = 1;
   c(i,i) = 1;
   c(i,2:i-1) = c(i-1,1:i-2) + c(i-1,2:i-1);
end

coef = c(row,1:row);