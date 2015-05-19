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
                        cell_cv_to_tri)
 
                    
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
  build_kexact_stencil;

  % ******* Important output ********* higher_order_recon
  % reconstruct the higher order solution
  % Sets up higher_order_recon(n).coef  % coefficients of higher order reconstruction
  cell.lhs_set = 0; % Reset stencil
  [higher_order_recon, ~] = reconstruct_solution(cell, fit_type, 1);

% ***********************************************
% ***********************************************
% Compute the jacobian for the original triangle mesh
for i = 1:cell_grid.ncells
    x = vertex(cell_grid.nodes(i,2:4),1);
    y = vertex(cell_grid.nodes(i,2:4),2);
    J(:,:,i) = compute_triangle_jacobian(x,y);
end



% cell_grid_to_cv = zeros(cell.ncells,1);
% for i = 1:cell.ncells
%     for j = 1:cell.nodes(i,1)
%       test_mat = repmat( vertex(cell.nodes(i,j+1),1:2), [cell_grid.ncells,1]);
%       I = find(  test_mat == cell_grid.xc(:,1:2) );
%       if length(I)>0
%           lcgc = length( find(cell_grid_to_cv(i,:)~=0) ) + 1;
%           cell_grid_to_cv(i,lcgc) = I(1);
%       end
%       
%     end
% end
% 
% I = find(cell_grid_to_cv(1,:)~=0);
% Jave = mean( J(:,:,cell_grid_to_cv(1,I)),3 );
% Jave = J(:,:,1);

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
build_kexact_stencil;
cell_te = cell;

% cell_smooth stores the solver order information for the smooth grid
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

% WAIT. WAIT. WAIT. I'm not quite sure if my process is correct.

cell_smooth = cell;
face_smooth = face;
for n = 1:cell.ncells
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

  % Compute scaling
  scale = sqrt( cell_old.volume(n)/dxideta );
%   scale = 1;
  vertex_scaled(:,1) = (vertex(:,1)-cell_smooth.xc(icell,1) )*scale + cell_old.xc(n,1);
  vertex_scaled(:,2) = (vertex(:,2)-cell_smooth.xc(icell,2) )*scale + cell_old.xc(n,2);
  vertex_scaled(:,3) = 0;

  for i = 1:cell.ncells
      cell.volume(i) = cell.volume(i)*scale^2;
      I = find( cell.nodes(i,:)~=0 );
      cell.xc(i,:) = mean( vertex_scaled( cell.nodes(i,I), 1:2), 1);
  end
  
  for i = 1:length(face)
      face(i).area = face(i).area*scale;
  end
%       nnodes = cell.nodes(i,1);
%       inodes = cell.nodes(i,2:nnodes+1);
%       cell.xc(i,1) = mean( vertex_scaled(inodes,1) );
%       cell.xc(i,2) = mean( vertex_scaled(inodes,2) );
%       cell.volume(i) = 0;
%       
%       I = find(cell.faces(n,:)~=0);
%       nf = length(I);
%       f = cell.faces(n,I);
%       for j = f
%           fn = face_smooth(f).nodes;
%           ax = vertex_scaled(fn(1),1);
%           ay = vertex_scaled(fn(1),2);
%           bx = vertex_scaled(fn(2),1);
%           by = vertex_scaled(fn(2),2);
%           cx = cell.xc(i,1);
%           cy = cell.xc(i,2);
%           
% %           a = compute_area(ax,ay,bx,by,cx,cy);
%           a = abs( triangle_area([ax,ay,0], [bx,by,0], [cx,cy,0]) );
%           cell.volume(i) = cell.volume(i) + a;
%       end
%       
%   end
  
  [face, ~] =  compute_face_data_from_triangulation(cell.nodes(:,2:end),vertex_scaled);
%   for i = 1:length(face_smooth)
%          fn = face_smooth(i).nodes;
%          ds_vec = vertex_scaled(fn(2),:) - vertex_scaled(fn(1),:);
%          face_smooth(i).area = sqrt( sum(ds_vec.^2) );
%          face_smooth(i).normal(1:2) = [ ds_vec(2), -ds_vec(1) ]/face(i).area;
%   end
  
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
%   [cell, face, vertex, icell, cell_phy] = construct_te_smooth_stencil(...
%                                  cell_old, face_old, vertex_old,...
%                                  te_order, ...
%                                  kexact_type, ...
%                                  fit_type, ...
%                                  flux_integral_order, base_order, n );
%   imax_local = 2*(te_order + 1)+1;
%   row = (imax_local-1)/2;
%   cntr = (imax_local)*row + row + 1;
%   ivec = (icell-1)*(imax_local) + icell;
                             
  % Need to shift vertex_scaled => vertex_scaled - xc_{cur_cell}   
  cell.mms_source = analytic_flux(vertex_scaled, cell, face, exact_flux, neq,...
                                  flux_integral_order, icell);

% mms_source_test = analytic_flux(vertex_old, cell_old, face_old, exact_flux, neq,...
%                                   flux_integral_order,n);


% SOMETHING IS WRONG HERE.
  
pte(1,:) = cell_te.reconstruction_param.px;
pte(2,:) = cell_te.reconstruction_param.py;

% [recon_soln_eval, ~] = reconstruct_solution(cell, fit_type, 1);
[moment] = compute_reconstruction_moments(vertex_scaled, cell, face, pte );
[~, ~, Aeval] = compute_reconstruction_lhs(stencil, icell, moment, pte, cell.xc, 1);

% for i = 1:size(A,1)
%     A(i,:) = A(i,:)./wij(i);
% end


% Step 5: Evaluate the higher order reconstruction over the smooth cell stencils.
for i = 1:size(higher_order_recon(1).coef,2)
    cell_eval = Aeval*higher_order_recon(n).coef(:,i);
    cell.soln(stencil,i) = cell_eval;  
end


  [smooth_recon, ~] = reconstruct_solution(cell, fit_type, 1, cell.stencil(icell).cells);


%   kexact_order = base_order;
%   kexact_type  = 'kexact_extended';%base_type;
%   fit_type = 'kexact';%base_fit_type;
%   kexact_type  = base_type;
%   fit_type = base_fit_type;
  
%   setup_reconstruction;
% %   build_kexact_stencil;
% 
% if cell_old.cell_type(n) == 0
%   sten = base_order+2;
%   sten_width = floor(sten/2);
%   cnt2=1;
%   for jj = row:row+2
%       for ii = row:row+2
% % %   for jj = 1:imax_local
% % %       for ii = 1:imax_local
%           iivec = (jj-1)*(imax_local) + ii;
% %           
% %           ilow = ii - sten_width;
% %           ihigh = ilow + sten;
% %           if (ilow<1)
% %               ilow = 1;
% %               ihigh = ilow + sten;
% %           elseif (ihigh>imax_local)
% %               ihigh = imax_local;
% %               ilow = ihigh - sten;
% %           end
% %           
% %           jlow = jj - sten_width;
% %           jhigh = jlow + sten;
% %           if (jlow<1)
% %               jlow = 1;
% %               jhigh = jlow + sten;
% %           elseif (jhigh>imax_local)
% %               jhigh = imax_local;
% %               jlow = jhigh - sten;
% %           end
% %           
% %           cnt = 1;
% %           for jjj = jlow : jhigh
% %               for iii = ilow : ihigh              
% %                   iiivec = (jjj-1)*(imax_local) + iii;
% %                   sub_sten(cnt) = iiivec;
% %                   cnt = cnt + 1;
% %               end
% %           end  
% %             cell.stencil(iivec).cells = sub_sten;
% % 
%             recon_sten(cnt2)=iivec;
%             cnt2 = cnt2 + 1;
%       end
%   end
%   
% elseif cell_old.cell_type(n) == 1
% 
%     cell_faces = cell.faces(icell,:);
%     recon_sten(1) = icell;
%     for jj = 1:length(cell_faces)
%         recon_sten(jj+1) = setdiff([face(cell_faces(jj)).cell_neg, face(cell_faces(jj)).cell_plus],icell);
%     end
% end
% 
%   cell.lhs_set = 0; % Reset stencil
%   [reconstruction, ~] = reconstruct_solution(cell, fit_type, recon_sten);
%   cell.reconstruction = reconstruction;
% 
% 

  cell.reconstruction = smooth_recon;
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
  
data_write = vertex_to_cell_average(cell.soln,cell_tri);
write_vtk_solution( vertex_scaled, cell_tri, data_write, 'soln', ['soln-te-',num2str(n),'.vtk'],'w',var ) 

data_write = vertex_to_cell_average(cell.mms_source,cell_tri);
write_vtk_solution( vertex_scaled, cell_tri, data_write, 'source', ['soln-te-',num2str(n),'.vtk'],'a',eq ) 


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
