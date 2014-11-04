function flux = flux_vanleer( ul, ur, normal)
% Computes vanleer riemann solver flux
% Inputs :
%          ul : left primitive variable state
%          ur : right primitive variable state
%      normal : face normal

% HARDWIRE : gamma = 1.4
gamma = 1.4;

%Compute left state
U = ul;
vel = normal(1)*U(2) + normal(2)*U(3);
a = sound_speed( U(4), U(1));
M = vel/a;

FL = zeros(1,4);
if ( abs(M) < 1 )
  fmass = 0.25*U(1)*a*(M+1)^2;
  fenergy = 0.5*fmass*( 1/(gamma^2-1)*( (gamma-1)*vel + 2*a )^2 ...
            + U(2)^2 + U(3)^2 - vel^2 );

  FL(1) = fmass;
  FL(2) = fmass*( normal(1)*( -vel + 2*a)/gamma + U(2) );
  FL(3) = fmass*( normal(2)*( -vel + 2*a)/gamma + U(3) );
  FL(4) = fenergy;

elseif ( M >= 1 )
  FL(1) = U(1)*vel;
  FL(2) = U(1)*vel*U(2) + normal(1)*U(4);
  FL(3) = U(1)*vel*U(3) + normal(2)*U(4);
  FL(4) = U(1)*vel*( gamma/(gamma-1)*U(4)/U(1) + 0.5*(U(2)^2 + U(3)^2) );
  
end


U = ur;
vel = normal(1)*U(2) + normal(2)*U(3);
a = sound_speed( U(4), U(1));
M = vel/a;

FR = zeros(1,4);
if ( abs(M) < 1 )
  fmass = 0.25*U(1)*a*(M+1)^2;
  fenergy = 0.5*fmass*( 1/(gamma^2-1)*( (gamma-1)*vel + 2*a )^2 ...
            +U(2)^2 + U(3)^2 - vel^2 );

  FR(1) = fmass;
  FR(2) = fmass*( normal(1)*( -vel + 2*a)/gamma + U(2) );
  FR(3) = fmass*( normal(2)*( -vel + 2*a)/gamma + U(3) );
  FR(4) = fenergy;

elseif ( M >= 1 )
  FR(1) = U(1)*vel;
  FR(2) = U(1)*vel*U(2) + normal(1)*U(4);
  FR(3) = U(1)*vel*U(3) + normal(2)*U(4);
  FR(4) = U(1)*vel*( gamma/(gamma-1)*U(4)/U(1) + 0.5*(U(2)^2 + U(3)^2) );
  
end

flux = FL + FR;
