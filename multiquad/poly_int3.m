function f = poly_int3(x0, x1, coefs, exponents)
% Usage:
% f = poly_int3(x0, x1, coefs, exponents)
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
% that can be integrated.
%
% This is the most optimized version of the integration by parts series. 
% The derivatives are pre computed to reduce computational expense and the 
% terms are combined using Pascal's triangle.
%
% Timing comparisons.
% Product of 4 4th order polynomials. 
% poly_int:  13.95  seconds
% poly_int2:  0.074 seconds
% poly_int3:  0.015 seconds
%
% Product of 5 4th order polynomials. 
% poly_int:  ----  seconds
% poly_int2: 0.15  seconds
% poly_int3: 0.018 seconds


%% If only one equation integrate and move on
if size(coefs,1)==1
    f = intf(x1, coefs(1,:), exponents(1,:),1)-intf(x0, coefs(1,:), exponents(1,:),1);
    return
end

%% Generate Pascals triangle
pmax = sum(exponents(1:end-1,end))+1;
P(1,1)   = 1;
P(2,1:2) = 1;
for i = 3:pmax
    P(i,1) = 1;
    P(i,2:i-1) = P(i-1,1:i-2)+P(i-1,2:i-1);
    P(i,i) = 1;
end         
         
%% Precompute derivatives
df1 = zeros(size(coefs,1),pmax);          
df0 = zeros(size(coefs,1),pmax);
for eq = 1:size(coefs,1)
    for term = 1:max(exponents(:,end))+2
        df1(eq,term) = df(x1, coefs(eq,:), exponents(eq,:),term-1);
        df0(eq,term) = df(x0, coefs(eq,:), exponents(eq,:),term-1);
    end
end

%% Setup first term
t1 = df1(1,:);
t0 = df0(1,:);


%% Product rule for the third or more functions
for eq = 2:size(coefs,1)-1
   for p = 1:pmax
      s1(p) = 0; s0(p) = 0;
      for n = 1:p
        s1(p) = s1(p) + P(p,n)*df1(eq,n)*t1(p-n+1);
        s0(p) = s0(p) + P(p,n)*df0(eq,n)*t0(p-n+1);
      end
   end
   t1 = s1;
   t0 = s0;
end

%% Compute integral terms
f1 = 0;
f0 = 0;
for p = 1:pmax
   f1 = f1 + (-1)^(p-1)*t1(p)*intf(x1, coefs(end,:), exponents(end,:),p);
   f0 = f0 + (-1)^(p-1)*t0(p)*intf(x0, coefs(end,:), exponents(end,:),p);
end

%% Compute final integration by parts integral
f = f1-f0;


