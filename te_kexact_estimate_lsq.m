function te_standard = te_kexact_estimate_lsq(cell, face, vertex, ..., ...
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
%   [higher_order_recon, ~] = reconstruct_solution(cell, fit_type, 3);
%   cell.reconstruction = higher_order_recon;
  neq = size(cell.soln,2);
  nunknowns = length(cell.reconstruction_param.px);
  for i = 1:cell.ncells
      sten = cell.stencil(i).cells;

    
      
      
      [~, ~, ~, ~, Aeval] = reconstruction_lhs(sten, i, cell.reconstruction_param, cell.xc, fit_type, 0, 1);
      higher_order_recon(i).Aeval = Aeval;
      
      %LSQ fit
      A = Aeval;
      b = cell.soln(sten,:);
      Ainv = (A'*A)^(-1)*A';

      % Polynomial fit
%       r = sqrt( (cell.xc(:,1)-cell.xc(i,1)).^2  + (cell.xc(:,2)-cell.xc(i,2)).^2  );
%       [rsort,Isort] = sort(r,1,'ascend');
% 
%       sten_poly = Isort(1:nunknowns);
%       cnt = 1;
%       for j = 1:length(sten_poly)
%           I = find(sten_poly(j)==sten);
%           A(cnt,:) = Aeval(I,:);
%           b(cnt,:) = cell.soln(sten(I),:);
%           cnt = cnt + 1;
%       end
%       Ainv = A^-1;

      
    higher_order_recon(i).coef = Ainv*b;

  end


  cell.reconstruction = higher_order_recon;
  

  for i = 1:cell.ncells
      I = find(cell.stencil(i).cells==i);
      cell_old.soln(i,:) = cell.reconstruction(i).Aeval(I,:)*cell.reconstruction(i).coef;
  end
  cell.soln = cell_old.soln;

  
  
  
  
  
  
  
  
  
  
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
  

cell_old.mms_source = mms_source_recon;

cell = cell_old;
[reconstruction, cell.lhs_set] = reconstruct_solution(cell, fit_type, 1);
cell.reconstruction = reconstruction;

%Compute left and right face states
% Not effecient but can only do this in matlab
face_out = compute_left_and_right_state(vertex, cell, face,...
                                      analytic_soln);
face = face_out;

%Compute residual
te_standard = compute_residual( cell, face, flux );


te_adjusted = te_standard;
if vertex_centered;   data_write = vertex_to_cell_average(te_adjusted,cell_grid); else data_write=te_standard; end;
  write_vtk_solution( vertex, cell, data_write, 'te_lsq', 'grid.vtk', 'a', eq )
  
te_smooth = te_adjusted;
primal_soln = cell_old.soln;
save('te_estimate_lsq.mat','te_smooth','primal_soln')
