%% quadratic_conserved_fit


syms xim1 xi xip1 xip2 xip3 xip4 xip5 uim1 ui uip1 uip2 uip3 uip4 sp sm dx v1 v2 v3 v4 v5 reals


vol = 1/32;
xi=-v1-v2;
xip1=-v2;
xip2 = 0;
xip3=xip2+v3;
xip4=xip3+v4;
xip5=xip4+v5;

% Third order
A = [(xip1-xi), 1/2*(xip1^2-xi^2), 1/3*(xip1^3-xi^3), 1/4*(xip1^4-xi^4), 1/5*(xip1^5-xi^5)
     (xip2-xip1), 1/2*(xip2^2-xip1^2), 1/3*(xip2^3-xip1^3), 1/4*(xip2^4-xip1^4), 1/5*(xip2^5-xip1^5)
     (xip3-xip2), 1/2*(xip3^2-xip2^2), 1/3*(xip3^3-xip2^3), 1/4*(xip3^4-xip2^4), 1/5*(xip3^5-xip2^5)
     (xip4-xip3), 1/2*(xip4^2-xip3^2), 1/3*(xip4^3-xip3^3), 1/4*(xip4^4-xip3^4), 1/5*(xip4^5-xip3^5)
     (xip5-xip4), 1/2*(xip5^2-xip4^2), 1/3*(xip5^3-xip4^3), 1/4*(xip5^4-xip4^4), 1/5*(xip5^5-xip4^5)];
 
b = [ui*(xip1-xi);uip1*(xip2-xip1);uip2*(xip3-xip2);uip3*(xip4-xip3);uip4*(xip5-xip4)];

a=A\b;

Anum = subs(A,{v1,v2,v3,v4,v5},{vol,vol,vol,vol,vol})
bnum = subs(b,{v1,v2,v3,v4,v5,ui,uip1,uip2,uip3,uip4},{vol,vol,vol,vol,vol,2,1.99921152114,1.999153882007,1.99909228555,1.99902646022})

anum = Anum\bnum;

intp = @(x) anum(1)*x+anum(2)*x.^2/2+anum(3)*x.^3/3+anum(4)*x.^4/4+anum(5)*x.^5/5;
(intp(vol)-intp(0))/vol