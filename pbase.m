function [f] = pbase(x,coef,px,py)

f = sum( coef'.*x(1).^px.*x(2).^py );
