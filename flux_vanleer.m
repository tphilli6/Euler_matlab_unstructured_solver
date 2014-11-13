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

%disp([vel, M])

FL = zeros(1,4);
if ( abs(M) < 1 )

  fmass = 0.25*U(1)*a*(M+1)^2;
  %fenergy = 0.5*fmass*( 1/(gamma^2-1)*( (gamma-1)*vel + 2*a )^2 ...
  %          + U(2)^2 + U(3)^2 - vel^2 );
  fenergy = 0.5*fmass*( 1/(gamma^2-1)*( (gamma-1)*vel + 2*a )^2 );
  %FL(1) = fmass;
  %FL(2) = fmass*( normal(1)*( -vel + 2*a)/gamma + U(2) );
  %FL(3) = fmass*( normal(2)*( -vel + 2*a)/gamma + U(3) );
  %FL(4) = fenergy;

  FL(1) = fmass;
  FL(2) = fmass*( ( (gamma-1)*vel + 2*a)/gamma );
  FL(3) = fmass*( ( (gamma-1)*vel + 2*a)/gamma );
  FL(4) = fenergy;

%if (nn==1 || nn==2 || nn==4 || nn==18 )
%nn
%vel
%a
%U
%  disp([fmass, fenergy, FL])
%end

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

%disp([vel,M])

FR = zeros(1,4);
if ( abs(M) < 1 )

  fmass = -0.25*U(1)*a*(M-1)^2;
  %fenergy = 0.5*fmass*( 1/(gamma^2-1)*( (gamma-1)*vel - 2*a )^2 ...
  %          +U(2)^2 + U(3)^2 - vel^2 );
  fenergy = 0.5*fmass*( 1/(gamma^2-1)*( (gamma-1)*vel + 2*a )^2 );

  %FR(1) = fmass;
  %FR(2) = fmass*( normal(1)*( -vel - 2*a)/gamma + U(2) );
  %FR(3) = fmass*( normal(2)*( -vel - 2*a)/gamma + U(3) );
  %FR(4) = fenergy;

  FR(1) = fmass;
  FR(2) = fmass*( ( (gamma-1)*vel + 2*a)/gamma );
  FR(3) = fmass*( ( (gamma-1)*vel + 2*a)/gamma );
  FR(4) = fenergy;

%if (nn==1 || nn==2 || nn==4 || nn==18 )
%vel
%a
%U
%  disp([fmass, fenergy, FR])
%end

elseif ( M <= -1 )

%  disp([U, vel, gamma])
  FR(1) = U(1)*vel;
  FR(2) = U(1)*vel*U(2) + normal(1)*U(4);
  FR(3) = U(1)*vel*U(3) + normal(2)*U(4);
  FR(4) = U(1)*vel*( gamma/(gamma-1)*U(4)/U(1) + 0.5*(U(2)^2 + U(3)^2) );
  
end
%disp(FL)
%disp(FR)

%if (nn==2 || nn==4 || nn==16 || nn==17 )
%  disp([FL+FR])
%end
flux = FL + FR;
