function soln = karman_trefftz_airfoil(x, nout, rhoinf, pinf, vinf, alpha, mu, te_angle)

%     rhoinf    = argin(2)
%     pinf      = argin(3)
%     vinf      = argin(4)
%     alpha     = argin(5)
%     mu        = cmplx(argin(6), argin(7), cd)

    alpha_rad = alpha*pi/180;
    R = sqrt( (1-real(mu)).^2 + imag(mu).^2);
    gamma = 4*pi*vinf*R*sin( alpha_rad + asin( imag(mu)/R ) );

%     xp = x(:,1)*4-2;
%     yp = x(:,2)*4;
    z = complex(x(:,1), x(:,2));
    zeta = kt_transform(z, te_angle);
    
    wtilde = cylinder_solution(zeta, vinf, alpha_rad, gamma, R, mu);
    w = wtilde./( 1 - 1./zeta.^2 );
    I=find(isnan(w));
    w(I)=0;

    switch nout
        case(1)
            soln = rhoinf*ones(size(w));
        case(2)
            soln = real(w);
        case(3)
            soln = -imag(w);
        case(4)
            % Bernouilli's to compute pressure
            u = real(w);
            v = -imag(w);
            soln = pinf + (vinf.^2 - u.^2 - v.^2)/2;
        otherwise
            soln = 0;
    end
    
end
    
function zeta = kt_transform(z, te_angle)
        
    n = 2 - te_angle/pi;

    zeta = -( ((z-n)./(z+n)).^(1/n) + 1) ./ ( ((z-n)./(z+n)).^(1/n) - 1);
end

function z = kt_transform_reverse(zeta, te_angle)
        
    n = 2 - te_angle/pi;

    z = n*(  (1+1./zeta).^n + (1-1./zeta).^n )./(  (1+1./zeta).^n - (1-1./zeta).^n );
end

function w = cylinder_solution(zeta, vinf, alpha_rad, gamma, R, mu)

  alpha_im = complex(0, alpha_rad);
  gamma_im = complex(0, gamma);

  w = vinf*exp(-alpha_im) + gamma_im./(2*pi*(zeta-mu)) - vinf*R.^2*exp(alpha_im)./(zeta-mu).^2;

end