function [V] = euler_analytic_solution(x, func)
% evaluates the input analytic solution for euler equations

V(1, 1) = func.rho(x);
V(1, 2) = func.u(x);
V(1, 3) = func.v(x);
V(1, 4) = func.p(x);
