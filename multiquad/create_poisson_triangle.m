function [P,exponents] = create_poisson_triangle(max_order)
max_order=5;

P = zeros(max_order+1);
exponents = zeros(max_order+1,max_order+1,2);

P(1,1) = 1;
exponents(1,1,:) = 0;

P(2,1:2) = 1;
exponents(2,1:2,1) = [1,0];
exponents(2,1:2,2) = [0,1];
for i = 3:max_order+1
    P(i,1) = 1;
    P(i,2:i-1) = P(i-1,1:i-2)+P(i-1,2:i-1);
    P(i,i) = 1;
    
    exponents(i,1:i,1) = i-1:-1:0;
    exponents(i,1:i,2) = 0:i-1;
end