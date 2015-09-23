function te_standard = te_kexact_estimate(cell, face, vertex, ..., ...
                                        kexact_type, ...
                                        fit_type, ...
                                        te_order,...
                                        flux_integral_order,...
                                        analytic_soln,...
                                        flux,...
                                        exact_flux,...
                                        vertex_centered,...
                                        cell_grid)
  
eq = {'mass','xmtm', 'ymtm', 'nrgy'};

cell_old = cell;


  % Step 1: compute a lsq reconstruction over a large stencil
  % - setup_reconstruction sets up the polynomial for the reconstruction
    kexact_order = te_order;
    setup_reconstruction;
    build_kexact_stencil;


  % ******* Important output ********* higher_order_recon
  % reconstruct the higher order solution
  % Sets up higher_order_recon(n).coef  % coefficients of higher order reconstruction
  cell.lhs_set = 0; % Reset stencil
  [higher_order_recon, ~] = reconstruct_solution(cell, fit_type, 3);

  cell.reconstruction = higher_order_recon;
  
  
  % compute discrete TE (i.e. Uses flux scheme)
  face_out = compute_left_and_right_state(vertex, cell, face, analytic_soln); face = face_out;


  % This is just a regular higher order estimate
  te_standard = -compute_residual( cell, face, flux );

te_smooth = te_standard;
primal_soln = cell_old.soln;
save('te_estimate_kexact.mat','te_smooth','primal_soln')
  
if vertex_centered;   data_write = vertex_to_cell_average(te_standard,cell_grid); else data_write=te_standard; end;
write_vtk_solution( vertex, cell, data_write, 'te_k_disc', 'grid.vtk', 'a', eq )
  



  % This makes an adjustment to the discrete residual to investigate the
  % differences in TE reconstruction and solution reconstruction
  
   source_order = max([cell.reconstruction_param.px, cell.reconstruction_param.py]);
   for n = 1:numel(cell.volume)
        func.rho = @(x) pbase(x, higher_order_recon(n).coef(:,1), cell.reconstruction_param.px, cell.reconstruction_param.py);
        func.u = @(x)   pbase(x, higher_order_recon(n).coef(:,2), cell.reconstruction_param.px, cell.reconstruction_param.py);
        func.v = @(x)   pbase(x, higher_order_recon(n).coef(:,3), cell.reconstruction_param.px, cell.reconstruction_param.py);
        func.p = @(x)   pbase(x, higher_order_recon(n).coef(:,4), cell.reconstruction_param.px, cell.reconstruction_param.py);
        exact_flux = @(x, normal) euler_mms_flux(x, func, normal);
        neq = 4;
        
        vertex_shift(:,1) = vertex(:,1) - cell.xc(n,1);
        vertex_shift(:,2) = vertex(:,2) - cell.xc(n,2);
        vertex_shift(:,3) = vertex(:,3);
        source = analytic_flux(vertex_shift, cell, face, exact_flux, neq,...
                                        source_order, n);
        mms_source_recon(n,:) = source(n,:);
   end
  
for n = 1:size(mms_source_recon,2);
    te_cont(:,n) = (cell.mms_source(:,n) - mms_source_recon(:,n)).*cell.volume'; 
end

if vertex_centered;   data_write = vertex_to_cell_average(te_cont,cell_grid); else data_write=te_cont; end;
  write_vtk_solution( vertex, cell, data_write, 'te_k_cont', 'grid.vtk', 'a', eq )

te_smooth = te_cont;
primal_soln = cell_old.soln;
save('te_estimate_kexact_cont.mat','te_smooth','primal_soln')


% Some test stuff
return
cell_te = cell;
% cell.mms_source = cell_te.mms_source;

  cell = cell_old;
  %cell.mms_source = mms_source_recon;
  nold = numel(cell.reconstruction(1).coef(:,1));
  % Sets the TE reconstruction coef to the solution reconstruction coef.
for i = 1:numel(cell.volume)
  cell.reconstruction(i).coef = cell_te.reconstruction(i).coef(1:nold,:);
end
% compute discrete residual (i.e. Uses flux scheme)
face_out = compute_left_and_right_state(vertex, cell, face, analytic_soln); face = face_out;


% This is just a regular higher order estimate
te_adjustment = compute_residual( cell, face, flux );
% te_adjustment = mms_source_recon;
   
if vertex_centered;   data_write = vertex_to_cell_average(te_adjustment,cell_grid); else data_write=te_standard; end;
  write_vtk_solution( vertex, cell, data_write, 'te_k_ad', 'grid.vtk', 'a', eq )
  
te_smooth = te_adjustment;
primal_soln = cell_old.soln;
save('te_estimate_kexact_disc_const.mat','te_smooth','primal_soln')
  
  
te_adjusted = te_standard-te_adjustment;
if vertex_centered;   data_write = vertex_to_cell_average(te_adjusted,cell_grid); else data_write=te_standard; end;
  write_vtk_solution( vertex, cell, data_write, 'te_k_adjed', 'grid.vtk', 'a', eq )
  
te_smooth = te_adjusted;
primal_soln = cell_old.soln;
save('te_estimate_kexact_disc_const_adjed.mat','te_smooth','primal_soln')
%   
%   px = cell.reconstruction_param.px;
%   py = cell.reconstruction_param.py;
%   for i = 1:size(cell.mms_source,2)
%     te_cont(:,i) = -cell.mms_source(:,i).*cell.volume';
%   end
%   
%   for i = 1:cell.ncells
%       xc = cell.xc(i,:);
% %       xc = [0,0];
%       coef = cell.reconstruction(i).coef;
%       
%       nodes = cell.nodes(i,2:end);
%       I = find(nodes~=0);
%       nodes = [nodes(I),nodes(1)];
%       
%      
%       for j = 1:length(nodes)-1
%             x1 = vertex(nodes(j),1:2);
%             x2 = vertex(nodes(j+1),1:2);
%             
%             ds_vec = x2 - x1;
%             area = sqrt( sum(ds_vec.^2) );
%             normal(1:2) = [ ds_vec(2), -ds_vec(1) ]/face(i).area;
%           
%             f_cell = 0;
%             f_cell2 = 0;
%             for ii = 1:length(wquad)
% 
%                 % HARDWIRE : linear face mapping
%                 % x mapping
%                 x(1) = (x2(1) - x1(1)) .* xquad(ii) + x1(1);
%                 x(2) = (x2(2) - x1(2)) .* xquad(ii) + x1(2);
% 
%                 V = ( (x(1)-xc(1)).^px.*(x(2)-xc(2)).^py )*coef;
%                 
%                 f_cell = f_cell + flux(V, V, normal)*wquad(ii);
%                 f_cell2 = f_cell2 + exact_flux(V, normal)*wquad(ii);
%             end
%           
%             te_cont(i,:) = te_cont(i,:) + f_cell*area;
%           
%           
%       end
%   end
%   
%   te_cont = -te_cont;
%   
%   if vertex_centered;   data_write = vertex_to_cell_average(te_cont,cell_grid); else data_write=te_cont; end;
%   write_vtk_solution( vertex, cell, data_write, 'te_kexact_cont', 'grid.vtk', 'a', eq )
%   
  