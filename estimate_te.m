function te = estimate_te(cell, face, vertex, kexact_order, resid_type, flux_integral_order, ...
                          kexact_type, fit_type, flux_func, reimann_solver, analytic_soln)
% Estimates the truncation error using any number of different options
% Input:  
%        cell : cell structure with relevant data
%      vertex : cell vertex physical points
%        face : face structure
% kexact_order: polynomial order for kexact reconstruction
%  resid_type : residual type ('Lh', 'L') 'Primary difference is Riemann solver
% flux_integral_order : quadrature order to integrate flux
% kexact_type : Controls number of terms in polynomial ('kexact','kexact_extended')
%    fit_type : Type of reconstruction ('kexact', 'lsq', 'lsq_true')


  % Sets up the reconstruction
  % Requires setup_mapping be called
  % Uses kexact_order, flux_integral_order
  % Sets up the cell integral, sets variables px, py, Ai, Axc, reconstruction_param.xquad[wquad]
  setup_reconstruction;
  % Sets up the stencil used for reconstruction
  % Requires setup_mapping and setup_reconstruction be called
  % Uses fit_type
  build_kexact_stencil;
 
  % reconstructs solution
  % returns the coefficients of polynomial set in setup_reconstruction
  [reconstruction, cell.lhs_set] = reconstruct_solution(cell, fit_type);
  cell.reconstruction = reconstruction;

  %Compute left and right face states
  % Not effecient but can only do this in matlab
  % Evaluates said polynomial at specified quadrature points set in call to setup_reconstruction
  face_out = compute_left_and_right_state(vertex, cell, face,...
                                          kexact_order, analytic_soln);
  face = face_out;
 
  %Compute residual
  if strcmp(resid_type,'Lh')
    % Discrete residual
    % Lh(.) = L(.) + te(.)
    % L(I_h u_h) ~ Lh_bar(I_h u_h)
    % te(.) ~ -Lh_bar(.)
    te = compute_residual( cell, face, reimann_solver );
%    te = -te;
  elseif strcmp(resid_type,'L')
    % Continuous residual
    % Lh(.) = L(.) = te(.)
    % te(.) = -L(I_h u_h)
    te = compute_continuous_residual( cell, face, flux_func );
    te = -te;
  end
