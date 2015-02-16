function func = setup_mms_crossterm(mms)

% Supersonic manufactured solution

l = 1;

switch mms
case (1)
rho = [1, ...
     -0.15, 2, 1/3, ...
     0.1, 1, -1/5, ...
     -0.05, 1, 1/4];

u = [800, ...
     50, 2/3, 1/2, ...
     30, 1.5, -1/8, ...
     -60, 1.5, 1/8];

v = [800, ...
     -75, 5/3, -1/2, ...
     40, 3, 0, ...
     60, 1.5, 1/8];

p = [100000, ...
     20000, 1/2, 1/2, ...
     50000, 1, -1/2, ...
     -10000, 3/4, -1/4];

case (2)
  rho = [1, ...
       0.15, 2, 1/3, ...
       -0.1, 1, -1/5, ...
       -0.05, 1, 1/4];
  
  u = [80, ...
       5, 2/3, 1/2, ...
       -3, 1.5, -1/8, ...
       -6, 1.5, 1/8];
  
  v = [80, ...
       -7.5, 5/3, -1/2, ...
       4.0, 3, 0, ...
       6.0, 1.5, 1/8];
  
  p = [100000, ...
       -20000, 1/2, 1/2, ...
       50000, 1, -1/2, ...
       -10000, 3/4, -1/4];
end
func.rho = @(x) base_function(x, rho, l);
func.u = @(x) base_function(x, u, l);
func.v = @(x) base_function(x, v, l);
func.p = @(x) base_function(x, p, l);



function f = base_function(x,c,l)

% f = c(1) + c(2)*x(1).^2 + c(5)*x(2).^2;% + c(8)*x(1)*x(2);
% f = c(1) + c(2)*x(1)*x(2) + c(5)*x(1).^2*x(2).^2 + c(8)*x(1).^3*x(2);% + c(8)*x(1)*x(2);

f = c(1) + c(2)*sin(c(3)*pi*x(1)/l + c(4)*pi) ...
        + c(5)*sin(c(6)*pi*x(2)/l + c(7)*pi) ...
        + c(8)*sin(c(9)*pi*x(1)*x(2)/l^2 + c(10)*pi); 
