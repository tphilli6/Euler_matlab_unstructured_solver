function [f, dpdx, dpdy] = peval(x, coef, px, py)
% Evaluates a function 
 f = sum( coef'.*x(1).^px.*x(2).^py );

dpx = px-1;
I = (dpx>=0);
coefdx = coef'.*I; %Remove the coef if derivative disappears
dpx = max(dpx,0);
dpdx =  sum( coefdx.*x(1).^dpx.*x(2).^py );


dpy = py-1;
I = (dpy>=0);
coefdy = coef'.*I; %Remove the coef if derivative disappears
dpy = max(dpy,0);
dpdy =  sum( coefdy.*x(1).^px.*x(2).^dpy );


