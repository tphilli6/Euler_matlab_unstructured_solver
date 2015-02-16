% Script to estimate the truncation error for all available options


write_vtk_solution( vertex, cell, cell.te_exact, 'exact-te', 'te.vtk','w',eq )

te_exact_l2 = sqrt( sum( cell.te_exact.^2 ) / size(cell.te_exact,1) );
% Generic setup, some variables are over written and some are not used
equation_setup


kexact_order_vec=[1,2,3,4]; 
resid_type_vec={'L','Lh'};
flux_integral_order_vec=[4];
kexact_type_vec={'kexact','kexact_extended'};
fit_type_vec={'kexact','lsq','lsq-true'};

load('output_solution.mat');

kexact_order = kexact_order_vec(4);
resid_type=resid_type_vec{1};
flux_integral_order=flux_integral_order_vec(1);
kexact_type = kexact_type_vec{1};
fit_type = fit_type_vec{2};


te = estimate_te(cell, face, vertex, kexact_order, resid_type, flux_integral_order, kexact_type,...
                 fit_type, eq_flux, flux, analytic_soln);

write_vtk_solution( vertex, cell, te, 'estimated-te', 'te.vtk','a',eq )
te_estimate_l2 = sqrt( sum( te.^2 ) / size(te,1) );


% Writes l2 norms of exact and estimate to file
fid = fopen('te_norms.dat','w');
fprintf(fid,'%12.8e  ', te_exact_l2); fprintf(fid,'\n');
fprintf(fid,'%12.8e  ', te_estimate_l2); fprintf(fid,'\n');

