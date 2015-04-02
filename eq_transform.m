% syms x1 y1 x2 y2 x3 y3 c xc yc reals
clear all
clc

[ti, xi, c_cntr] = generate_equilateral_mesh(4);

% for i = 1:length(ti(:,1));
figure(1)
% i=64;
i=c_cntr;
hold off
plot_cells(ti,xi)
hold on
plot([xi(ti(i,:),1);xi(ti(i,1),1)],[xi(ti(i,:),2);xi(ti(i,1),2)],'g-o')
ximesh = xi;




x = [0, 1, 0.5];
y = [0, 0, sqrt(3)/2];

J = compute_triangle_jacobian(x,y);

xfun = @(xi,eta) mean(x) + J(1,1)*xi + J(1,2)*eta;
yfun = @(xi,eta) mean(y) + J(2,1)*xi + J(2,2)*eta;

xtrans(:,1) = xfun(ximesh(:,1), ximesh(:,2));
xtrans(:,2) = yfun(ximesh(:,1), ximesh(:,2));
 
xi = xtrans;
figure(2)
i=c_cntr;
hold off
plot_cells(ti,xi)
hold on
plot([xi(ti(i,:),1);xi(ti(i,1),1)],[xi(ti(i,:),2);xi(ti(i,1),2)],'g-o')
ximesh = xi;
axis equal

