function f = poly_int(x0, x1, coefs, exponents)
% Usage:
% f = poly_int(x0, x1, coefs, exponents)
%
% Inputs:
%         x0 = lower bound for the integral
%         x1 = upper bound for the integral
%      coefs = matrix of coefficients with increasing order. row =>
%              equation, column => term in equation
%  exponents = matrix of exponents corresponding to coefs with increasing
%              order
%
% =========================================================================
% This function integrates the product of several polynomials using
% integration by parts. There is no upper limit to the number of equations
% that can be integrated. However, using this function is not advised
% because the computational time increases exponentially with the number of
% equations.
%
% This is the first of the integration by parts series. The function is
% called recursively. It is the simplest algorithm to implement; however,
% it's painfully slow.
%
% Timing comparisons.
% Product of 3 4th order polynomials. 
% poly_int:  0.054  seconds
% poly_int2:  0.028 seconds
% poly_int3:  0.013 seconds
%
% Product of 4 4th order polynomials. 
% poly_int:  13.95  seconds
% poly_int2:  0.074 seconds
% poly_int3:  0.015 seconds
%
% Product of 5 4th order polynomials. 
% poly_int: >7376  seconds (Stopped because it was taking too long)
% poly_int2: 0.15  seconds
% poly_int3: 0.018 seconds



f0 = 1;
f1 = 1;
for n = 1:size(coefs,1)-1
   f0 = f0*poly_eval(x0, coefs(n,:), exponents(n,:) );
   f1 = f1*poly_eval(x1, coefs(n,:), exponents(n,:) );
end
f0 = f0*poly_eval(x0, coefs(end,:)./(exponents(end,:)+1), exponents(n,:) + 1);
f1 = f1*poly_eval(x1, coefs(end,:)./(exponents(end,:)+1), exponents(n,:) + 1);

f = f1-f0;

for n = 1:size(coefs,1)-1
    if exponents(n,end)~=0
    f = f - poly_int(x0, x1, [coefs(1:n-1,:); coefs(n,:).*exponents(n,:);coefs(n+1:end-1,:);coefs(end,:)./(exponents(end,:)+1)],...
            max([exponents(1:n-1,:); exponents(n,:)-1;exponents(n+1:end-1,:);exponents(end,:)+1],0) );
    end
end