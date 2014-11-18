function f = df(x, coef, exponent,n)

I = find(exponent-n>=0);
f = coef.*(factorial(exponent))./factorial( max(exponent-n,0)).*x.^max(exponent-n,0);
f = sum(f(I));