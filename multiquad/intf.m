function f = intf(x, coef, exponent,n)

f = sum(coef.*(factorial(exponent))./factorial( exponent+n).*x.^(exponent+n));
