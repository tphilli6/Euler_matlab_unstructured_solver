clc
clear all

equation_setup


refinement = 2;
kexact_order = 2;
grid_type = -1; %quad = 0, triangles = 1, mixed = 2, from file = -1 (predecided mix)
grid_in = 'equilateral.mesh';%'0012-1397.mesh';
[ti, xi] = generate_equilateral_mesh(kexact_order+1, 1); % generate and write 'equilateral.mesh'

source_term_order=kexact_order+1;
exact_order=kexact_order+1;
flux_integral_order=kexact_order+1;

for ref = 1:10
    scale = 1/(refinement.^(ref-1));

%     scale = 1;
    % generate_equilateral_mesh(4, 1); % generate and write 'equilateral.mesh'
    [~, ~, ~, vertex, face, cell_new,  cell_tri, cell_tri_to_cv] = read_unstructured_mesh(grid_in, vertex_centered, scale);

%     vertex(:,1) = vertex(:,1) + 1;
%     vertex(:,2) = vertex(:,2) + 1;
    
    
    % Computes the cell mapping required for quadratures
    cell = cell_new;
    setup_reconstruction;
    build_kexact_stencil;
    reconstruction=cell.reconstruction;

    fit = [1, cell.nbrs(1,1:cell.nnbr(1))];
    vol(ref) = cell.volume(1);


    for n = 1 : length(face); face(n).func = exact_fun; end % Assign exact function to faces
    cell.mms_source = analytic_flux(vertex, cell, face, exact_flux, neq, source_term_order); % Compute source term
    cell.exact = compute_generic_moments(vertex, cell, face, exact_fun );
    cell.exact2 = compute_analytic_exact(cell_tri, vertex, exact_order, analytic_soln, cell_tri_to_cv); % Compute exact solution
    cell.soln = cell.exact;

%     dx=0.01;
%     y=0.1;
%     for ij = 1:5
%         der_a=(exact_fun.intx([dx,y],exact_fun.cex(1,:)) -exact_fun.intx([-dx,y],exact_fun.cex(1,:) ))/(2*dx);
%         der = exact_fun.rho([0,y]);
%         der_e(ij) = der_a - der;
%         dx = dx/2;
%     end
    
    cell.reconstruction = reconstruct_solution(cell, fit_type, 1);

    face = compute_left_and_right_state(vertex, cell, face, analytic_soln);

%     for i=1:12
%         nap=cell.faces(1,i);
%         face(nap).ul(:,1)-face(nap).ur(:,1)
%     end
%     face(cell.faces(1,12)).nodes = cell.nodes(1,[12,1]);

    cell.te_exact = compute_residual( cell, face, flux );

    for jj = 1:size(cell.te_exact,2)
        cell.te_exact(:,jj) = cell.te_exact(:,jj)./cell.volume';
    end

    flux_int(ref,:) = abs(cell.te_exact(1,:));
%     flux_int(ref,:) = sqrt( sum( (cell.te_exact(1:7,:)).^2 )/7 );

end

h = vol.^(1/2);%(vol./vol(end)).^(1/2);

rp = (vol(1:end-1)./vol(2:end)).^(1/2);
p = log(  flux_int(1:end-1,:)./flux_int(2:end,:) )./log( repmat(rp',[1,size(flux_int,2)]) );

ref_h = ([h(1),h(end)]/h(end)).^(1/2);
ref_line = 0.1*max(max(flux_int)).*[ref_h.^(kexact_order+1)];


figure(2)
subplot(1,2,1)

loglog(h,flux_int(:,1), 'r-o',...
       h,flux_int(:,2), 'g-o',...
       h,flux_int(:,3), 'k-o',...
       h,flux_int(:,4), 'b-o',...
       ref_h, ref_line, 'k--');
title('Error')

subplot(1,2,2)
semilogx(h(2:end),p(:,1), 'r-o',...
       h(2:end),p(:,2), 'g-o',...
       h(2:end),p(:,3), 'k-o',...
       h(2:end),p(:,4), 'b-o',...
       [h(1),h(end)], kexact_order+1*[1,1], 'k--');
axis([h(end),h(1),0,kexact_order+2])
title('Order of Convergence')