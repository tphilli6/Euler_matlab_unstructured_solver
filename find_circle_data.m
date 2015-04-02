%% Find circle
syms x1 y1 x2 y2 x3 y3

A=[ x2-x1, y2-y1
    x3-x1, y3-y1 ];

b=1/2*[ x2^2-x1^2 + y2^2 - y1^2
        x3^2-x1^2 + y3^2 - y1^2 ];
    
f0 = A\b;
x0 = simplify(expand(f0(1)));
y0 = simplify(expand(f0(2)));