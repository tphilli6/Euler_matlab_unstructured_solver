function [quality, r, X0] = mesh_quality_r_to_d(x,y)


X = [x, x(:,1)];
Y = [y, y(:,1)];

dx = X(:,2:end) - X(:,1:end-1);
dy = Y(:,2:end) - Y(:,1:end-1);

ds = sqrt(dx.^2+dy.^2);

d = min(ds')';

[r,x0,y0] = find_circle( x(:,1), y(:,1) ,...
                         x(:,2), y(:,2) ,...
                         x(:,3), y(:,3)  );
                    
quality = r./d;
X0 = [x0,y0];