function [Alhs,cond_num,A, wij, Aeval] = reconstruction_lhs(stencil, n, reconstruction_param, xc, ...
                                           fit_type, test, ge_terms)
% This function is the companion of reconstruct_solution. The exact routine is needed to check the condition number of the stencil

Alhs = 0;

% A = build_lhs(stencil, reconstruction, fit_type);
p = [reconstruction_param.px; reconstruction_param.py];
[A, wij, Aeval] = compute_reconstruction_lhs(stencil,n, reconstruction_param.moment, p, xc, ge_terms);


% If a lsq fit then compute the left hand side for a non-square matrix
% Ax = b => x = (A'A)^(-1) A' b
if strcmp(fit_type(1:3),'lsq')

  if strcmp(fit_type,'lsq')
%     Alhs = A(2:end,2:end);
    cond_num = rcond(A'*A);
    if (test==0)
      Alhs = (A'*A)^-1*A';
    end
  elseif strcmp(fit_type,'lsq-true')
    cond_num = rcond(A'*A);
    if (test==0)
      Alhs = (A'*A)^-1*A';
    end
  end


%   cond(Alsq'*Alsq)
%   Alhs_temp = (Alsq'*Alsq)^-1*Alsq';
%   reconstruction(n).Ainv = Alhs_temp;

elseif strcmp(fit_type,'kexact')
  cond_num = rcond(A);
  if (test==0)
    Alhs = A^-1;
  end
%   cond(A)
%   reconstruction(n).Ainv = A^-1;

elseif strcmp(fit_type,'extended')
    cond_num = rcond(A'*A);
    if (test==0)
      Alhs = (A'*A)^-1*A';
    end

end

%if test==1
%A
%disp([stencil,cond_num])
%end
