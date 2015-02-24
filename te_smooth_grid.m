function te_smooth_grid(cell, face, vertex, ...
                        kexact_order, ...
                        kexact_type, ...
                        fit_type, ...
                        te_order,...
                        flux_integral_order,...
                        analytic_soln,...
                        flux)

eq = {'mass','xmtm', 'ymtm', 'nrgy'};
base_order = kexact_order;
base_type  = kexact_type;
base_fit_type = fit_type;
soln_old = cell.soln;
cell_old = cell;
face_old = face;
vertex_old = vertex;

% xmat = reshape(cell_old.xc(:,1),[8,8]);
% ymat = reshape(cell_old.xc(:,2),[8,8]);
% s = reshape(cell_old.soln(:,4), [8,8]);
% v = reshape(cell_old.volume, [8,8]);

  % Step 1: compute a lsq reconstruction over a large stencil
  % - setup_reconstruction sets up the polynomial for the reconstruction
  kexact_order = te_order;
  fit_type = 'lsq';%base_fit_type;
  extra_terms = 12; % manually add more cells to the stencil. = 0 will be a kexact fit. This can be variable to play with the smoothness of the fit

  setup_reconstruction;
  higher_order_recon_param = cell.reconstruction_param;
  build_kexact_stencil;

  % ******* Important output ********* higher_order_recon
  % reconstruct the higher order solution
  % Sets up higher_order_recon(n).coef  % coefficients of higher order reconstruction
  cell.lhs_set = 0; % Reset stencil
  [higher_order_recon, ~] = reconstruct_solution(cell, fit_type);
   cell.reconstruction = higher_order_recon;

  face_out = compute_left_and_right_state(vertex, cell, face,...
                                          analytic_soln);
  face = face_out;

  % This is just a regular higher order estimate
  te_standard = compute_residual( cell, face, flux );

  write_vtk_solution( vertex_old, cell_old, te_standard, 'te_estimate_k4', 'grid.vtk', 'a', eq )



  % Step 1: compute a lsq reconstruction over a large stencil
  % - setup_reconstruction sets up the polynomial for the reconstruction
  kexact_order = te_order;
  kexact_type  = base_type;
  fit_type = 'lsq';
  extra_terms = 12; % manually add more cells to the stencil. = 0 will be a kexact fit. This can be variable to play with the smoothness of the fit

  setup_reconstruction;
  higher_order_recon_param = cell.reconstruction_param;
  build_kexact_stencil;


  % ******* Important output ********* higher_order_recon
  % reconstruct the higher order solution
  % Sets up higher_order_recon(n).coef  % coefficients of higher order reconstruction
  cell.lhs_set = 0; % Reset stencil
  [higher_order_recon, ~] = reconstruct_solution(cell, fit_type);
   cell.reconstruction = higher_order_recon;
  
     % Copy the old cell setup for the same reconstruction stencil that was used for the solve
  cell = cell_old;
  
for n = 1:cell.ncells
fprintf('\nTE estimate for cell %4.0f\n',n);
  % Euler hardcode
  func.rho = @(x) pbase(x, higher_order_recon(n).coef(:,1), higher_order_recon_param.px, higher_order_recon_param.py);
  func.u = @(x)   pbase(x, higher_order_recon(n).coef(:,2), higher_order_recon_param.px, higher_order_recon_param.py);
  func.v = @(x)   pbase(x, higher_order_recon(n).coef(:,3), higher_order_recon_param.px, higher_order_recon_param.py);
  func.p = @(x)   pbase(x, higher_order_recon(n).coef(:,4), higher_order_recon_param.px, higher_order_recon_param.py);
  exact_flux = @(x, normal) euler_mms_flux(x, func, normal);
  neq = 4;

%   tic
%   higher_order_recon(n).coef;
% 
  cell = cell_old;
  vertex = vertex_old;
  face = face_old;
  
  kexact_order = base_order;
  kexact_type  =  base_type;
  fit_type = base_fit_type;
% 
% %   cell.mms_source = analytic_flux(vertex, cell, face, exact_flux, neq,...
% %                                   flux_integral_order, n);
% 
%                               
%                               
%   % Evaluate the higher order reconstruction for cell n over the surrounding unsmooth stencils
%   % Loops over the all the stencil cells included in stencil cell n
%   % NOTE: Several duplicate evaluations are done
%   for i = 1:length(cell.stencil(n).cells)
%     icell = cell.stencil(n).cells(i);
%     for j = 1:length(cell.stencil(icell).cells)
%       jcell = cell.stencil(icell).cells(j);
%       for nneq = 1:size(higher_order_recon(n).coef,2)
%         cell.soln(jcell,nneq) = higher_order_recon(jcell).Ai*higher_order_recon(n).coef(:,nneq);
%        % cell.reconstruction(jcell).coef = higher_order_recon(n).coef;
%       end
%     end
%   end
% % cell.soln(jcell,nneq)
% 
%   [reconstruction, ~] = reconstruct_solution(cell, base_fit_type, n);
%    cell.reconstruction = reconstruction;
% % reconstruction(n).coef
%   % Computes the left and right states over all the cells
%   % NOTE: a lot of duplicate calculations here
%   face_out = compute_left_and_right_state(vertex, cell, face,...
%                                           analytic_soln);
%   face = face_out;
%  
%   % Step 3: Evaluate the residual over the non-smooth cells
%   te_temp = compute_residual( cell, face, flux, n);
% 
%   te_non_smooth(n,:) = te_temp(n,:);
%   te_non_smooth(n,:) = cell.mms_source(n,:);

%   time1 = toc;
%   fprintf('Evaluate higher order reconstruction over non-smooth cells %8.6f\n', time1)
  tic
  % Step 4: Build a smooth grid classcell.ncells
  % Call to construct_te_smooth_stencil
  % sets
  %      subcell.nunknowns                  % polynomial unknowns
  %      subcell.reconstruction_param.xquad % quadrature points
  %      subcell.reconstruction_param.wquad % quadrature weights
  %      subcell.reconstruction_param.px    % x coefficients
  %      subcell.reconstruction_param.py    % y coefficients
  %      subcell.reconstuction.Ai           % integral over structured cell (smooth) geometry
  %      subcell.reconstruction.Axc         % reference integral location at cell center
  [cell, face, vertex, icell, cell_phy] = construct_te_smooth_stencil(...
                                 cell_old, face_old, vertex_old,...
                                 te_order, ...
                                 kexact_type, ...
                                 fit_type, ...
                                 flux_integral_order, n );
  imax_local = 2*(te_order + 1)+1;
  row = (imax_local-1)/2;
  cntr = (imax_local)*row + row + 1;
  ivec = (icell-1)*(imax_local) + icell;                             
                             
                             
  cell.mms_source = analytic_flux(vertex, cell, face, exact_flux, neq,...
                                  flux_integral_order, ivec);


  
  % Step 5: Evaluate the higher order reconstruction over the smooth cell stencils. 
  for i = 1:cell.ncells
    for nneq = 1:size(higher_order_recon(n).coef,2)
      cell.soln(i,nneq) = cell_phy.reconstruction(i).Ai*higher_order_recon(n).coef(:,nneq);
    end
  end

%   
% xmat = reshape(cell.xc(:,1),[imax_local,imax_local]);
% ymat = reshape(cell.xc(:,2),[imax_local,imax_local]);
% s = reshape(cell.soln(:,4), [imax_local,imax_local]);
% v = reshape(cell.volume, [imax_local,imax_local]);
% 
% subplot(1,2,1)
% surf(xmat,ymat,s)
% view([0,0,90])
% 
% subplot(1,2,2)
% surf(xmat,ymat,s./v)
% view([0,0,90])

  % Step 6: Setup
  % Setup the grid classes for a second-order reconstruction
%  [subcell, subface, subvertex, icell2] = construct_te_smooth_stencil(...
%                                           cell, face, vertex,...
%                                           base_order, ...
%                                           base_type, ...
%                                           base_fit_type, ...
%                                           flux_integral_order, n );


  kexact_order = base_order;
%   kexact_type  = 'kexact_extended';%base_type;
%   fit_type = 'kexact';%base_fit_type;
  kexact_type  = base_type;
  fit_type = base_fit_type;
  
  setup_reconstruction;
%   build_kexact_stencil;

  sten = base_order+2;
  sten_width = floor(sten/2);
  cnt2=1;
  for jj = row:row+2
      for ii = row:row+2
%   for jj = 1:imax_local
%       for ii = 1:imax_local
          iivec = (jj-1)*(imax_local) + ii;
          
          ilow = ii - sten_width;
          ihigh = ilow + sten;
          if (ilow<1)
              ilow = 1;
              ihigh = ilow + sten;
          elseif (ihigh>imax_local)
              ihigh = imax_local;
              ilow = ihigh - sten;
          end
          
          jlow = jj - sten_width;
          jhigh = jlow + sten;
          if (jlow<1)
              jlow = 1;
              jhigh = jlow + sten;
          elseif (jhigh>imax_local)
              jhigh = imax_local;
              jlow = jhigh - sten;
          end
          
          cnt = 1;
          for jjj = jlow : jhigh
              for iii = ilow : ihigh              
                  iiivec = (jjj-1)*(imax_local) + iii;
                  sub_sten(cnt) = iiivec;
                  cnt = cnt + 1;
              end
          end  
            cell.stencil(iivec).cells = sub_sten;

            recon_sten(cnt2)=iivec;
            cnt2 = cnt2 + 1;
      end
  end


  time3 = toc;
%   fprintf('Generate the smooth cell stencil for the regular order evaluation %8.6f\n', time3)
  tic

%  imax2_local = 2*(base_order + 1)+1;
%  imax_diff = icell - icell2;
%  cnt = 1;
%
%  for j2 = 1:imax2_local
%    for i2 = 1:imax2_local
%      i4 = i2 + imax_diff;
%      j4 = j2 + imax_diff;
%    
%      ivec = (j4-1)*(imax2_local) + i4;
%      subcell.soln(cnt,:) = local_soln(ivec,:);
%      cnt = cnt + 1;
%    end
%  end
%  cell.mms_source = zeros(size(cell.soln));

  % Compute the source term for the higher order curve fit on the smooth grid
  % Compute cell center
%  for i = 1:cell.ncells
%    nvertex = cell.nodes(i,1);
%    cindex = cell.nodes(i,2:end);
%    xc = 0;
%
%    for j = 1:nvertex
%      xc = vertex(cindex,:) + xc;
%    end
%    xc = xc/nvertex;
%
%    % Evaluate source term at cell center
%    cell.mms_source(i,:) = euler_source(xc, higher_order_recon(n).coef, higher_order_recon_param.px, higher_order_recon_param.px);
%
%  end

%  % Copy the reconstruction for cell i to all cells in the subcell
%  for i = 1:length(subcell)
%    subcell.reconstruction(i).Ainv = higher_order_recon(n).Ainv;
%    subcell.reconstruction(i).Axc = higher_order_recon(n).Axc;
%    subcell.reconstruction(i).Ai = higher_order_recon(n).Ai;
%  end

  cell.lhs_set = 0; % Reset stencil
  [reconstruction, ~] = reconstruct_solution(cell, fit_type, recon_sten);
  cell.reconstruction = reconstruction;



  % Step 6: Evaluate the residual over the smooth cells
  ncell = cell.faces(ivec,:);
  face_out = compute_left_and_right_state(vertex, cell, face,...
                                          analytic_soln, ivec);
  face = face_out;

  % Step 7: Estimate the truncation error over the reconstruction
  % This is just a regular higher order estimate
  te_temp = compute_residual( cell, face, flux, ivec);
  
  time4 = toc;
%   fprintf('Compute the residual over the smaller stencil smooth stencil %8.6f\n', time4)

  ivec = (icell-1)*(imax_local) + icell;
  te_smooth(n,:) = -te_temp(ivec,:);
%   te_smooth(n,:) = cell.soln(ivec,:);
%    te_smooth(n,:) = cell.mms_source(ivec,:);
% R = reshape(te_temp(:,1), [imax_local,imax_local]);
% surf(xmat(3:end-2,3:end-2),ymat(3:end-2,3:end-2),R(3:end-2,3:end-2))
% view([0,0,90])
%   fprintf('Total time for cell %4.0f: %8.6f\n', n, time1 + time2 + time3 + time4)
end

%   write_vtk_solution( vertex_old, cell_old, te_non_smooth, 'te_estimate_nonsmooth', 'grid.vtk', 'a', eq )
  write_vtk_solution( vertex_old, cell_old, te_smooth, 'te_estimate_smooth', 'grid.vtk', 'a', eq )
primal_soln = cell_old.soln;
save('te_estimate.mat','te_smooth','primal_soln')
