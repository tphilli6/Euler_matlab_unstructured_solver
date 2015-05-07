function te_smooth_grid(cell, face, vertex, ...
                        kexact_order, ...
                        kexact_type, ...
                        fit_type, ...
                        te_order,...
                        flux_integral_order,...
                        analytic_soln,...
                        flux)
 
itermax = 1;
                    
eq = {'mass','xmtm', 'ymtm', 'nrgy'};


% base_order = kexact_order;
% base_type  = kexact_type;
% base_fit_type = fit_type;
% soln_old = cell.soln;
% cell_old = cell;
% face_old = face;
% vertex_old = vertex;

% xmat = reshape(cell_old.xc(:,1),[8,8]);
% ymat = reshape(cell_old.xc(:,2),[8,8]);
% s = reshape(cell_old.soln(:,4), [8,8]);
% v = reshape(cell_old.volume, [8,8]);

  % Step 1: compute a lsq reconstruction over a large stencil
  % - setup_reconstruction sets up the polynomial for the reconstruction
  kexact_order = te_order;
  fit_type = base_fit_type;

  setup_reconstruction;
  higher_order_recon_param = cell.reconstruction_param;
  build_kexact_stencil;

  % ******* Important output ********* higher_order_recon
  % reconstruct the higher order solution
  % Sets up higher_order_recon(n).coef  % coefficients of higher order reconstruction
  cell.lhs_set = 0; % Reset stencil
  [higher_order_recon, ~] = reconstruct_solution(cell, fit_type);
%   for i = 1:length(higher_order_recon)
% %       higher_order_recon(i).coef(:,1) = [1, 0, -0.15, 0, 0, 0.1]';
% %       higher_order_recon(i).coef(:,2) = [800, 0, 50, 0, 0, 30]';
% %       higher_order_recon(i).coef(:,3) = [800, 0, -75, 0, 0, 40]';
% %       higher_order_recon(i).coef(:,4) = [100000, 0, 20000, 0, 0, 50000]';
%       higher_order_recon(i).coef(:,1) = [1, -0.15, 0, 0.1, 0, 0]';
%       higher_order_recon(i).coef(:,2) = [800, 50, 0, 30, 0, 0]';
%       higher_order_recon(i).coef(:,3) = [800, -75, 0, 40, 0, 0]';
%       higher_order_recon(i).coef(:,4) = [100000, 20000, 0, 50000, 0, 0]';
%   end
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
  fit_type = 'lsq-true';
 
  setup_reconstruction;
  higher_order_recon_param = cell.reconstruction_param;
  build_kexact_stencil;


  % ******* Important output ********* higher_order_recon
  % reconstruct the higher order solution
  % Sets up higher_order_recon(n).coef  % coefficients of higher order reconstruction
  cell.lhs_set = 0; % Reset stencil
  [higher_order_recon, ~] = reconstruct_solution(cell, fit_type);
%   for i = 1:length(higher_order_recon)
%       higher_order_recon(i).coef(:,1) = [1, 0, -0.15, 0, 0, 0.1]';
%       higher_order_recon(i).coef(:,2) = [800, 0, 50, 0, 0, 30]';
%       higher_order_recon(i).coef(:,3) = [800, 0, -75, 0, 0, 40]';
%       higher_order_recon(i).coef(:,4) = [100000, 0, 20000, 0, 0, 50000]';
%       higher_order_recon(i).coef(:,1) = [1, -0.15, 0, 0.1, 0, 0]';
%       higher_order_recon(i).coef(:,2) = [800, 50, 0, 30, 0, 0]';
%       higher_order_recon(i).coef(:,3) = [800, -75, 0, 40, 0, 0]';
%       higher_order_recon(i).coef(:,4) = [100000, 20000, 0, 50000, 0, 0]';
      
%   end
   cell.reconstruction = higher_order_recon;
  
%      % Copy the old cell setup for the same reconstruction stencil that was used for the solve
%   cell = cell_old;
  
% for iter = 1:itermax
  % Copy the old cell setup for the same reconstruction stencil that was used for the solve
%   cell = cell_old;
  
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

%   cell.mms_source = analytic_flux(vertex, cell, face, exact_flux, neq,...
%                                   flux_integral_order, n);

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
                                 flux_integral_order, base_order, n );
%   imax_local = 2*(te_order + 1)+1;
%   row = (imax_local-1)/2;
%   cntr = (imax_local)*row + row + 1;
%   ivec = (icell-1)*(imax_local) + icell;
                             
                             
  cell.mms_source = analytic_flux(vertex, cell, face, exact_flux, neq,...
                                  flux_integral_order, icell);

% mms_source_test = analytic_flux(vertex_old, cell_old, face_old, exact_flux, neq,...
%                                   flux_integral_order,n);
  
  % Step 5: Evaluate the higher order reconstruction over the smooth cell stencils.
  for i = 1:cell.ncells
    for nneq = 1:size(higher_order_recon(n).coef,2)
      cell.soln(i,nneq) = cell_phy.reconstruction(i).Ai*higher_order_recon(n).coef(:,nneq);
    end
  end


  kexact_order = base_order;
%   kexact_type  = 'kexact_extended';%base_type;
%   fit_type = 'kexact';%base_fit_type;
  kexact_type  = base_type;
  fit_type = base_fit_type;
  
  setup_reconstruction;
%   build_kexact_stencil;

if cell_old.cell_type(n) == 0
  sten = base_order+2;
  sten_width = floor(sten/2);
  cnt2=1;
  for jj = row:row+2
      for ii = row:row+2
% %   for jj = 1:imax_local
% %       for ii = 1:imax_local
          iivec = (jj-1)*(imax_local) + ii;
%           
%           ilow = ii - sten_width;
%           ihigh = ilow + sten;
%           if (ilow<1)
%               ilow = 1;
%               ihigh = ilow + sten;
%           elseif (ihigh>imax_local)
%               ihigh = imax_local;
%               ilow = ihigh - sten;
%           end
%           
%           jlow = jj - sten_width;
%           jhigh = jlow + sten;
%           if (jlow<1)
%               jlow = 1;
%               jhigh = jlow + sten;
%           elseif (jhigh>imax_local)
%               jhigh = imax_local;
%               jlow = jhigh - sten;
%           end
%           
%           cnt = 1;
%           for jjj = jlow : jhigh
%               for iii = ilow : ihigh              
%                   iiivec = (jjj-1)*(imax_local) + iii;
%                   sub_sten(cnt) = iiivec;
%                   cnt = cnt + 1;
%               end
%           end  
%             cell.stencil(iivec).cells = sub_sten;
% 
            recon_sten(cnt2)=iivec;
            cnt2 = cnt2 + 1;
      end
  end
  
elseif cell_old.cell_type(n) == 1

    cell_faces = cell.faces(icell,:);
    recon_sten(1) = icell;
    for jj = 1:length(cell_faces)
        recon_sten(jj+1) = setdiff([face(cell_faces(jj)).cell_neg, face(cell_faces(jj)).cell_plus],icell);
    end
end

  cell.lhs_set = 0; % Reset stencil
  [reconstruction, ~] = reconstruct_solution(cell, fit_type, recon_sten);
  cell.reconstruction = reconstruction;



  % Step 6: Evaluate the residual over the smooth cells
  ncell = cell.faces(icell,:);
  face_out = compute_left_and_right_state(vertex, cell, face,...
                                          analytic_soln, icell);
  face = face_out;

  % Step 7: Estimate the truncation error over the reconstruction
  % This is just a regular higher order estimate
  te_temp = compute_residual( cell, face, flux,icell);
  

%   ivec = (icell-1)*(imax_local) + icell;
  te_smooth(n,:) = -te_temp(icell,:);
%   te_smooth(n,:) = cell.soln(ivec,:);
%    te_smooth(n,:) = cell.mms_source(ivec,:);
% R = reshape(te_temp(:,1), [imax_local,imax_local]);
% surf(xmat(3:end-2,3:end-2),ymat(3:end-2,3:end-2),R(3:end-2,3:end-2))
% view([0,0,90])
%   fprintf('Total time for cell %4.0f: %8.6f\n', n, time1 + time2 + time3 + time4)

% fprintf('%23.15e %23.15e %23.15e %23.15e\n',te_smooth(n,:) - te_standard(n,:));
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
  write_vtk_solution( vertex_old, cell_old, te_smooth, 'te_estimate_smooth', 'grid.vtk', 'a', eq )
primal_soln = cell_old.soln;
save('te_estimate.mat','te_smooth','primal_soln')
