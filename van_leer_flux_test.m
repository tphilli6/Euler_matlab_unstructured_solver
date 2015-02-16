%% Test the van Leer flux function

rho = 1;
u = 800;
v = 800;
p = 100000;
normal = [sqrt(2)/2, sqrt(2)/2];
normal = [1,0];


ul = [rho,u,v,p];
ur = [rho,u,v,p]*0.95;

flux_ex = euler_flux(rho,u,v,p,normal);
flux = flux_vanleer(ul, ur, normal);
fprintf('Error in test 1: %12.4e %12.4e %12.4e %12.4e\n', flux-flux_ex)

flux_ex = euler_flux(rho*0.95,u*0.95,v*0.95,p*0.95,-normal);
flux = flux_vanleer(ul, ur, -normal);
fprintf('Error in test 2: %12.4e %12.4e %12.4e %12.4e\n', flux-flux_ex)


ul = [rho,-u,-v,p];
ur = [rho,-u,-v,p]*0.95;

flux_ex = euler_flux(rho*0.95, -u*0.95, -v*0.95, p*0.95, -normal)
flux = flux_vanleer(ur, ul, -normal)
fprintf('Error in test 3: %12.4e %12.4e %12.4e %12.4e\n', flux-flux_ex)

flux_ex = euler_flux(rho,-u,-v,p, normal)
flux = flux_vanleer(ul, ur, normal)
fprintf('Error in test 4: %12.4e %12.4e %12.4e %12.4e\n', flux-flux_ex)
