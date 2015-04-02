function [J] = compute_triangle_jacobian(x,y)

dx = 1;
dy = 0.5*tand(60);

xi = [0, dx, dx/2];
eta = [0, 0, dy];

xi = xi - mean(xi);
eta = eta - mean(eta);

% A = [1, xi(1), eta(1)
%      1, xi(2), eta(2)
%      1, xi(3), eta(3)];

% Precomputed Ainv of A;
x1 = xi(1); x2 = xi(2); x3 = xi(3); y1 = eta(1); y2 = eta(2); y3 = eta(3);
Ainv = [ (x2*y3 - x3*y2)/(x1*y2 - x2*y1 - x1*y3 + x3*y1 + x2*y3 - x3*y2), -(x1*y3 - x3*y1)/(x1*y2 - x2*y1 - x1*y3 + x3*y1 + x2*y3 - x3*y2), (x1*y2 - x2*y1)/(x1*y2 - x2*y1 - x1*y3 + x3*y1 + x2*y3 - x3*y2)
               (y2 - y3)/(x1*y2 - x2*y1 - x1*y3 + x3*y1 + x2*y3 - x3*y2),       -(y1 - y3)/(x1*y2 - x2*y1 - x1*y3 + x3*y1 + x2*y3 - x3*y2),       (y1 - y2)/(x1*y2 - x2*y1 - x1*y3 + x3*y1 + x2*y3 - x3*y2)
              -(x2 - x3)/(x1*y2 - x2*y1 - x1*y3 + x3*y1 + x2*y3 - x3*y2),        (x1 - x3)/(x1*y2 - x2*y1 - x1*y3 + x3*y1 + x2*y3 - x3*y2),      -(x1 - x2)/(x1*y2 - x2*y1 - x1*y3 + x3*y1 + x2*y3 - x3*y2)];

x = reshape(x,[numel(x),1]);
y = reshape(y,[numel(y),1]);

ax = Ainv*x;
ay = Ainv*y;

J = [ax(2), ax(3)
     ay(2), ay(3)];

