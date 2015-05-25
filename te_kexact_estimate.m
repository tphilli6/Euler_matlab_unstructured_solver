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

  if vertex_centered;   data_write = vertex_to_cell_average(te_standard,cell_grid); else data_write=te_standard; end;
  write_vtk_solution( vertex, cell, data_write, 'te_kexact_discrete', 'grid.vtk', 'a', eq )
  
  px = cell.reconstruction_param.px;
  py = cell.reconstruction_param.py;
  for i = 1:size(cell.mms_source,2)
    te_cont(:,i) = -cell.mms_source(:,i).*cell.volume';
  end
  
  for i = 1:cell.ncells
      xc = cell.xc(i,:);
%       xc = [0,0];
      coef = cell.reconstruction(i).coef;
      
      nodes = cell.nodes(i,2:end);
      I = find(nodes~=0);
      nodes = [nodes(I),nodes(1)];
      
     
      for j = 1:length(nodes)-1
            x1 = vertex(nodes(j),1:2);
            x2 = vertex(nodes(j+1),1:2);
            
            ds_vec = x2 - x1;
            area = sqrt( sum(ds_vec.^2) );
            normal(1:2) = [ ds_vec(2), -ds_vec(1) ]/face(i).area;
          
            f_cell = 0;
            f_cell2 = 0;
            for ii = 1:length(wquad)

                % HARDWIRE : linear face mapping
                % x mapping
                x(1) = (x2(1) - x1(1)) .* xquad(ii) + x1(1);
                x(2) = (x2(2) - x1(2)) .* xquad(ii) + x1(2);

                V = ( (x(1)-xc(1)).^px.*(x(2)-xc(2)).^py )*coef;
                
                f_cell = f_cell + flux(V, V, normal)*wquad(ii);
                f_cell2 = f_cell2 + exact_flux(V, normal)*wquad(ii);
            end
          
            te_cont(i,:) = te_cont(i,:) + f_cell*area;
          
          
      end
  end
  
  te_cont = -te_cont;
  
  if vertex_centered;   data_write = vertex_to_cell_average(te_cont,cell_grid); else data_write=te_cont; end;
  write_vtk_solution( vertex, cell, data_write, 'te_kexact_cont', 'grid.vtk', 'a', eq )
  
  