function te_smooth_grid(cell, face, vertex, ...
                        kexact_order, ...
                        kexact_type, ...
                        fit_type, ...
                        te_order,...
                        flux_integral_order,...
                        analytic_soln,...
                        flux,...
                        vertex_centered,...
                        cell_grid,...
                        face_grid,...
                        vertex_grid,...
                        cell_cv_to_tri,...
                        exact_flux_orig)
 
                    
itermax = 1;
                    
eq = {'mass','xmtm', 'ymtm', 'nrgy'};
var = {'rho','u','v','p'};

base_order = kexact_order;
base_type  = kexact_type;
base_fit_type = fit_type;
soln_old = cell.soln;
cell_old = cell;
face_old = face;
vertex_old = vertex;


  % Step 1: compute a lsq reconstruction over a large stencil
  % - setup_reconstruction sets up the polynomial for the reconstruction
  kexact_order = te_order;
  fit_type = base_fit_type;

  setup_reconstruction;
  higher_order_recon_param = cell.reconstruction_param;
  
  extra_terms = cell.nunknowns*2;
  fit_type = 'extended';
  build_kexact_stencil;

  % ******* Important output ********* higher_order_recon
  % reconstruct the higher order solution
  % Sets up higher_order_recon(n).coef  % coefficients of higher order reconstruction
  cell.lhs_set = 0; % Reset stencil
  [higher_order_recon, ~] = reconstruct_solution(cell, fit_type, 1);
  
  % Relaxes LSQ constraint by overwritting the Ainv previously computed.
  cell.lhs_set = 1;
  for i = 1:cell.ncells
      A = higher_order_recon(i).Aeval;
      cell.reconstruction(i).Ainv = (A'*A)^(-1)*A';      
  end
[higher_order_recon, ~] = reconstruct_solution(cell, fit_type, 1);
cell.lhs_set = 0; % Reset stencil

% ***********************************************
% ***********************************************
% Compute the jacobian for the original triangle mesh
for i = 1:cell_grid.ncells
    x = vertex(cell_grid.nodes(i,2:4),1);
    y = vertex(cell_grid.nodes(i,2:4),2);
    J(:,:,i) = compute_triangle_jacobian(x,y);
end


[ti, xi, icell] = generate_equilateral_mesh(te_order);
% volume(x)/det(J) = volume(xi)
% scale jacobian so that the areas are the same even if the shape is not.
dxideta = 0.866025403784439; %6.4952; %det( [ 0.5, 0; 0.25, 0.5*tand(60)] ); % Volume of equilateral triangle from generate_equilateral_mesh()
icell = floor(size(xi,1)/2);

vertex = [xi,zeros(size(xi(:,1)))];

[face_tri, cell_tri] =  compute_face_data_from_triangulation(ti,vertex);
[cell, face, vertex, cell_all, all_to_cv] = compute_vertex_centered_nonsense(cell_tri, face_tri, vertex);
cell = find_neighbor_cells(cell, face);
cell.soln = zeros( cell.ncells, size(cell_old.soln,2) );

% cell_te, face_te, and vertex_te store information that's required for
% higher order reconstruction
kexact_order = te_order;
kexact_type = base_type;
setup_reconstruction;
% build_kexact_stencil;
cell_te = cell;

% cell_smooth stores the solver order information for the smooth grid
fit_type = base_fit_type;
kexact_order = base_order;
kexact_type = base_type;
setup_reconstruction;
build_kexact_stencil;

stencil=cell.stencil(icell).cells;
for i = cell.stencil(icell).cells
    stencil = [stencil, cell.stencil(i).cells];
end
stencil = unique(stencil);

I = find(stencil~=icell);
stencil = [icell, stencil(I)];
% stencil = 1:length(cell.volume);

source_order = max([higher_order_recon_param.px, higher_order_recon_param.py]);

pmap = [0,1,0; 0, 0, 1];




cell_smooth = cell;
face_smooth = face;
for n = 1:cell_old.ncells
    cell = cell_smooth;
    face = face_smooth;


    fprintf('\nTE estimate for cell %4.0f\n',n);
    % Euler hardcode for source term
    func.rho = @(x) pbase(x, higher_order_recon(n).coef(:,1), higher_order_recon_param.px, higher_order_recon_param.py);
    func.u = @(x)   pbase(x, higher_order_recon(n).coef(:,2), higher_order_recon_param.px, higher_order_recon_param.py);
    func.v = @(x)   pbase(x, higher_order_recon(n).coef(:,3), higher_order_recon_param.px, higher_order_recon_param.py);
    func.p = @(x)   pbase(x, higher_order_recon(n).coef(:,4), higher_order_recon_param.px, higher_order_recon_param.py);
    exact_flux = @(x, normal) euler_mms_flux(x, func, normal);
    neq = 4;

    % Curve fit volume
%     clear Avol
    nbrs = cell_old.nbrs(n,:);
    I = find(nbrs~=0);
    nbrs = [n,nbrs(I)];
    cell_volume = mean(cell_old.volume(nbrs));
%     V = cell_old.volume(nbrs)';
%     xc = cell_old.xc(nbrs,:);
%     for i = 1:length(V)
%         Avol(i,:) = [1, xc(i,1), xc(i,2)];        
%     end
%     avol = Avol\V;
%     cell_volume  = avol(1) + avol(2)*cell_old.xc(n,1) + avol(3)*cell_old.xc(n,2);    
    
   

    % Compute scaling
    scale = sqrt( cell_volume/dxideta );
    %   scale = 1;
    vertex_scaled_0(:,1) = (vertex(:,1)-cell_smooth.xc(icell,1) )*scale;
    vertex_scaled_0(:,2) = (vertex(:,2)-cell_smooth.xc(icell,2) )*scale;
    vertex_scaled(:,1) = vertex_scaled_0(:,1) + cell_old.xc(n,1);
    vertex_scaled(:,2) = vertex_scaled_0(:,2) + cell_old.xc(n,2);
    vertex_scaled(:,3) = 0;

    for i = 1:cell.ncells
        cell.volume(i) = cell.volume(i)*scale^2;
        I = find( cell.nodes(i,2:end)~=0 );
        cell.xc(i,:) = mean( vertex_scaled( cell.nodes(i,I+1), 1:2), 1);
    end

    for i = 1:length(face)
        face(i).area = face(i).area*scale;
    end

  
  [face, ~] =  compute_face_data_from_triangulation(cell.nodes(:,2:end),vertex_scaled);

                             
  % Need to shift vertex_scaled => vertex_scaled - xc_{cur_cell}   
  cell.mms_source = analytic_flux(vertex_scaled_0, cell, face, exact_flux, neq,...
                                  source_order);
                  
% TEMP FOR TESTING
% cell.mms_source = analytic_flux(vertex_scaled, cell, face, exact_flux_orig, neq,...
%                                 flux_integral_order);
  
pte(1,:) = cell_te.reconstruction_param.px;
pte(2,:) = cell_te.reconstruction_param.py;

% [recon_soln_eval, ~] = reconstruct_solution(cell, fit_type, 1);
[moment] = compute_reconstruction_moments(vertex_scaled, cell, face, pte, stencil );
[~, ~, Aeval] = compute_reconstruction_lhs(stencil, icell, moment, pte, cell.xc, 1);

% Step 5: Evaluate the higher order reconstruction over the smooth cell stencils.
for i = 1:size(higher_order_recon(1).coef,2)
    cell_eval = Aeval*higher_order_recon(n).coef(:,i);
    cell.soln(stencil,i) = cell_eval;  
end

% TEMP FOR TESTING ========================================================
% cell.soln = compute_analytic_exact(cell_all, vertex_scaled, 4, analytic_soln, all_to_cv);


[cell.reconstruction_param.moment] = compute_reconstruction_moments(vertex_scaled, cell, face, p, stencil );
  [smooth_recon, ~] = reconstruct_solution(cell, fit_type, 1, [icell, cell.nbrs(icell,:)] );

  cell.reconstruction = smooth_recon;
  % Step 6: Evaluate the residual over the smooth cells
  ncell = cell.faces(icell,:);
  face_out = compute_left_and_right_state(vertex_scaled, cell, face,...
                                          analytic_soln, icell);
  face = face_out;

  % Step 7: Estimate the truncation error over the reconstruction
  % This is just a regular higher order estimate
  te_temp = compute_residual( cell, face, flux,icell);


%   ivec = (icell-1)*(imax_local) + icell;
  te_smooth(n,:) = te_temp(icell,:);
  
% data_write = vertex_to_cell_average(cell.soln,cell_tri);
% write_vtk_solution( vertex_scaled, cell_tri, data_write, 'soln', ['soln-te-',num2str(n),'.vtk'],'w',var ) 
% 
% data_write = vertex_to_cell_average(cell.mms_source,cell_tri);
% write_vtk_solution( vertex_scaled, cell_tri, data_write, 'source', ['soln-te-',num2str(n),'.vtk'],'a',eq ) 
% 
% data_write = vertex_to_cell_average(-te_temp,cell_tri);
% write_vtk_solution( vertex_scaled, cell_tri, data_write, 'resid', ['soln-te-',num2str(n),'.vtk'],'a',var ) 

end



% for nn = 1:length(te_smooth(1,:))
% nn=1;
%    tvolume = sum(cell_old.volume);
%    cbar = cell_old.volume*abs(te_smooth(:,nn))/tvolume;
%    
%    cdiff(iter,nn) = max( abs( abs(te_smooth(:,nn)) - cbar) );
%    
%    A = diag( abs(te_smooth(:,nn)) );
%    A = A - ones(size(A));
%    
%    b = (cbar - tvolume)*ones(length(te_smooth(:,1)),1);
%    
%    temp = A\b;
%    
%    temp_scaled = temp/sum(temp);
%    
%    new_volume(:,nn) = cell_old.volume' + 0.05*(temp_scaled - cell_old.volume');
   

% xx = reshape(cell_old.xc(:,1),[16,16]);
% yy = reshape(cell_old.xc(:,2),[16,16]);
% % v = reshape(cell_old.volume,[16,16]);
% v = reshape(new_volume(:,nn),[16,16]);
% figure(1)
% surf(xx,yy,v);
% view([0,0,90])


% figure(2)
% t = reshape(abs(te_smooth(:,nn)),[16,16]);
% surf(xx,yy,t/max(max(t)))
% view([0,0,90])
% end

% fprintf('%12.4e %12.4e %12.4e %12.4e\n',cdiff(iter,:));

% cell_old.volume = mean(new_volume,2)';






% xx = reshape(cell_old.xc(:,1),[16,16]);
% yy = reshape(cell_old.xc(:,2),[16,16]);
% v = reshape(cell_old.volume,[16,16]);
% surf(xx,yy,v);
% view([0,0,90])


% end

% fprintf('%12.4e\n',cdiff);

%   write_vtk_solution( vertex_old, cell_old, te_non_smooth, 'te_estimate_nonsmooth', 'grid.vtk', 'a', eq )

data_write = vertex_to_cell_average(te_smooth, cell_grid, cell_old);
  write_vtk_solution( vertex_grid, cell_grid, data_write, 'te_estimate_smooth', 'grid.vtk', 'a', eq )
primal_soln = cell_old.soln;
save('te_estimate.mat','te_smooth','primal_soln')
