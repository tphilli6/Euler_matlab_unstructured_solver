function te_standard = te_kexact_test(cell, face, vertex, ..., ...
                                        kexact_type, ...
                                        fit_type, ...
                                        te_order,...
                                        flux_integral_order,...
                                        analytic_soln,...
                                        flux,...
                                        exact_flux,...
                                        vertex_centered,...
                                        cell_grid)
                            
                            
% This function is trying to reconstruct the higher order terms. It assumes 
% that the solver reconstruction stencil and the te reconstruction stencil
% are the same size. This may not be necessary, but it's a first cut. The
% final reconstruction from the solver is evaluated over the stencil and
% subtracted from the solution values. The TE estimate reconstruction
% removes these solver terms from the reconstruction LHS and solves for the
% coefficients of only the higer order terms. This completes the Taylor
% series expansion so that there are shared terms. Otherwise, the solver
% reconstruction and the TE reconstruction have different initial terms.

cell_old = cell;

eq = {'mass','xmtm', 'ymtm', 'nrgy'};




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
  % The previous call to reconstruct_solution is the standard TE
  % reconstruction. It is called so that Aeval is setup for the TE
  % reconstruction. The actual coefficients are disgarded.
  % TEST KERNEL
  
  old_nunknowns = size(cell_old.reconstruction(1).coef(:,1));
  neq = size(cell_old.reconstruction(1).coef(1,:));
  for i = 1:cell.ncells
      ueval = cell_old.reconstruction(i).Aeval*cell_old.reconstruction(i).coef;
      urhs = cell_old.soln(cell_old.stencil(i).cells,:) - ueval;
      Aeval = cell.reconstruction(i).Aeval(:,old_nunknowns+1:end);
      A = Aeval;
      
      
    I = find(cell.stencil(i).cells==i);
    A(I,:) = A(I,:)./A(I,1);
    row = A(I,:);
    cnt=1;
    clear Aout
    for j = 1:size(A,1)
        if j~=I
            
            c = A(j,1)/row(1);
            
            Aout(cnt,:) = A(j,:)-c*row;
            cnt = cnt + 1;
        else
%             Aout(cnt,:) = zeros(size(row));
            Aout(cnt,:) = A(j,:);
            cnt = cnt + 1;
            
        end
    end
    A = Aout;
      
      
      
      
      
      for j = 1:size(A,2)
         A(:,j) = A(:,j).*cell.reconstruction(i).wij;
      end
      for j = 1:size(urhs,2)
         urhs(:,j) = urhs(:,j).*cell.reconstruction(i).wij;
      end
      Alhs = (A'*A)\A';

      
      for j = 1:neq
          coef_te(:,j) = Alhs*urhs(:,j);
          cell.reconstruction(i).coef(:,j) = [cell_old.reconstruction(i).coef(:,j); coef_te];
          
      end
      
  end
  
  % Now that the TE reconstruction is partly constrained to correspond with
  % the numerical solution, go on with computing a residual as normal.
  
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
   te = analytic_flux(vertex_shift, cell, face, exact_flux, neq,...
                                    source_order, n);
   te_standard(n,:) = -(te(n,:)-cell.mms_source(n,:)).*cell.volume(i);
   end
  
  % compute discrete TE (i.e. Uses flux scheme)
%   face_out = compute_left_and_right_state(vertex, cell, face, analytic_soln); face = face_out;


  % This is just a regular higher order estimate
%   te_standard = compute_residual( cell, face, flux );

  if vertex_centered;   data_write = vertex_to_cell_average(te_standard,cell_grid); else data_write=te_standard; end;
  write_vtk_solution( vertex, cell, data_write, 'te_k_disc_const', 'grid.vtk', 'a', eq )
  
  
te_smooth = te_standard;
primal_soln = cell_old.soln;
save('te_estimate_kexact_const.mat','te_smooth','primal_soln')

  
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
  
  