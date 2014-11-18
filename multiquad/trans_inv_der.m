clc
clear all

syms dxdxi dydxi dzdxi dxdeta dydeta dzdeta dxdzeta dydzeta dzdzeta dxi deta dzeta dx dy dz  dxidx dxidy dxidz detadx detady detadz zdetadx zdetady zdetadz J Jinv reals

%3D
fprintf('3D ================================\n')

Tinv = [ dxdxi dxdeta dxdzeta
      dydxi dydeta dydzeta
      dzdxi dzdeta dzdzeta];

T = [ dxidx dxidy dxidz
      detadx detady detadz
      zdetadx zdetady zdetadz]; 
  
ddxi = [dxi; deta; dzeta];
ddx = [dx; dy; dz];



ddxideq = det([ddx,Tinv(:,2:end)])/Jinv;
ddxideq2 = dxidx*dx+dxidy*dy+dxidz*dz;

disp(subs(ddxideq2,{dx,dy,dz},{1,0,0}))
fprintf(' = ')
disp(subs(ddxideq,{dx,dy,dz},{1,0,0}))

disp(subs(ddxideq2,{dx,dy,dz},{0,1,0}))
fprintf(' = ')
disp(subs(ddxideq,{dx,dy,dz},{0,1,0}))

disp(subs(ddxideq2,{dx,dy,dz},{0,0,1}))
fprintf(' = ')
disp(subs(ddxideq,{dx,dy,dz},{0,0,1}))

fprintf('\n\n')
ddetadeq = det([Tinv(:,1),ddx,Tinv(:,3)])/Jinv;
ddetadeq2 = detadx*dx+detady*dy+detadz*dz;

disp(subs(ddetadeq2,{dx,dy,dz},{1,0,0}))
fprintf(' = ')
disp(subs(ddetadeq,{dx,dy,dz},{1,0,0}))

disp(subs(ddetadeq2,{dx,dy,dz},{0,1,0}))
fprintf(' = ')
disp(subs(ddetadeq,{dx,dy,dz},{0,1,0}))

disp(subs(ddetadeq2,{dx,dy,dz},{0,0,1}))
fprintf(' = ')
disp(subs(ddetadeq,{dx,dy,dz},{0,0,1}))


fprintf('\n\n')
dzdetadeq = det([Tinv(:,[1,2]),ddx])/Jinv;
dzdetadeq2 = zdetadx*dx+zdetady*dy+zdetadz*dz;

disp(subs(dzdetadeq2,{dx,dy,dz},{1,0,0}))
fprintf(' = ')
disp(subs(dzdetadeq,{dx,dy,dz},{1,0,0}))

disp(subs(dzdetadeq2,{dx,dy,dz},{0,1,0}))
fprintf(' = ')
disp(subs(dzdetadeq,{dx,dy,dz},{0,1,0}))

disp(subs(dzdetadeq2,{dx,dy,dz},{0,0,1}))
fprintf(' = ')
disp(subs(dzdetadeq,{dx,dy,dz},{0,0,1}))

fprintf('\n\nJ = ')
disp(det(T))
fprintf('\nJinv = ')
disp(det(Tinv))


%3D=>2D
fprintf('3D=>2D ================================\n')

ddxideq = det([ddx,Tinv(:,2:end)])/Jinv;
ddxideq2 = dxidx*dx+dxidy*dy+dxidz*dz;

disp(subs(ddxideq,{dx,dy,dz,dxdzeta, dydzeta, dzdxi, dzdeta, dzdzeta},{1,0,0,0,0,0,0,1}))
fprintf(' = ')
disp(subs(ddxideq2,{dx,dy,dz},{1,0,0}))

disp(subs(ddxideq,{dx,dy,dz,dxdzeta, dydzeta, dzdxi, dzdeta, dzdzeta},{0,1,0,0,0,0,0,1}))
fprintf(' = ')
disp(subs(ddxideq2,{dx,dy,dz},{0,1,0}))

% disp(subs(ddxideq,{dx,dy,dz,dxdzeta, dydzeta, dzdxi, dzdeta, dzdzeta},{0,0,1,0,0,0,0,1}))
% fprintf(' = ')
% disp(subs(ddxideq2,{dx,dy,dz},{0,0,1}))

fprintf('\n\n')
ddetadeq = det([Tinv(:,1),ddx,Tinv(:,3)])/Jinv;
ddetadeq2 = detadx*dx+detady*dy+detadz*dz;


disp(subs(ddetadeq,{dx,dy,dz,dxdzeta, dydzeta, dzdxi, dzdeta, dzdzeta},{1,0,0,0,0,0,0,1}))
fprintf(' = ')
disp(subs(ddetadeq2,{dx,dy,dz},{1,0,0}))

disp(subs(ddetadeq,{dx,dy,dz,dxdzeta, dydzeta, dzdxi, dzdeta, dzdzeta},{0,1,0,0,0,0,0,1}))
fprintf(' = ')
disp(subs(ddetadeq2,{dx,dy,dz},{0,1,0}))
% 
% disp(subs(ddetadeq,{dx,dy,dz,dxdzeta, dydzeta, dzdxi, dzdeta, dzdzeta},{0,0,1,0,0,0,0,1}))
% fprintf(' = ')
% disp(subs(ddetadeq2,{dx,dy,dz},{0,0,1}))


% fprintf('\n\n')
% dzdetadeq = det([Tinv(:,[1,2]),ddx])/Jinv;
% dzdetadeq2 = zdetadx*dx+zdetady*dy+zdetadz*dz;
% 
% disp(subs(dzdetadeq,{dx,dy,dz,dxdzeta, dydzeta, dzdxi, dzdeta, dzdzeta},{1,0,0,0,0,0,0,1}))
% fprintf(' = ')
% disp(subs(dzdetadeq2,{dx,dy,dz},{1,0,0}))
% 
% disp(subs(dzdetadeq,{dx,dy,dz,dxdzeta, dydzeta, dzdxi, dzdeta, dzdzeta},{0,1,0,0,0,0,0,1}))
% fprintf(' = ')
% disp(subs(dzdetadeq2,{dx,dy,dz},{0,1,0}))
% 
% disp(subs(dzdetadeq,{dx,dy,dz,dxdzeta, dydzeta, dzdxi, dzdeta, dzdzeta},{0,0,1,0,0,0,0,1}))
% fprintf(' = ')
% disp(subs(dzdetadeq2,{dx,dy,dz},{0,0,1}))

fprintf('\nJinv = ')
disp(subs(det(Tinv),{dxdzeta, dydzeta, dzdxi, dzdeta, dzdzeta},{0,0,0,0,1}))



%2D
fprintf('2D ================================\n')

Tinv = [ dxdxi dxdeta 
      dydxi dydeta];

T = [ dxidx dxidy
      detadx detady];  

ddxi = [dxi; deta];
ddx = [dx; dy];


ddxideq = det([ddx,Tinv(:,2)])/Jinv;
ddxideq2 = dxidx*dx+dxidy*dy+dxidz*dz;

disp(subs(ddxideq2,{dx,dy,dz},{1,0,0}))
fprintf(' = ')
disp(subs(ddxideq,{dx,dy,dz},{1,0,0}))

disp(subs(ddxideq2,{dx,dy,dz},{0,1,0}))
fprintf(' = ')
disp(subs(ddxideq,{dx,dy,dz},{0,1,0}))


fprintf('\n\n')
ddetadeq = det([Tinv(:,1),ddx])/Jinv;
ddetadeq2 = detadx*dx+detady*dy+detadz*dz;

disp(subs(ddetadeq2,{dx,dy,dz},{1,0,0}))
fprintf(' = ')
disp(subs(ddetadeq,{dx,dy,dz},{1,0,0}))

disp(subs(ddetadeq2,{dx,dy,dz},{0,1,0}))
fprintf(' = ')
disp(subs(ddetadeq,{dx,dy,dz},{0,1,0}))



fprintf('\n\nJ = ')
disp(det(T))
fprintf('\nJinv = ')
disp(det(Tinv))

% %1D
% fprintf('1D ================================\n')
% 
% Tinv = dxdxi;
%       
% 
% T = dxidx ;
%       
% 
% ddxi = [dxi];
% ddx = [dx];
% 
% 
% ddxideq = det([ddxi,T(:,2)])/J;
% disp(ddxideq)
% fprintf(' = ')
% disp(dxdxi*dxi+dxdeta*deta)
% 
% 
% ddxideq = det([ddx,Tinv(:,2)])/Jinv;
% ddxideq2 = dxidx*dx+dxidy*dy+dxidz*dz;
% 
% disp(subs(ddxideq2,{dx,dy,dz},{1,0,0}))
% fprintf(' = ')
% disp(subs(ddxideq,{dx,dy,dz},{1,0,0}))
% 
% disp(subs(ddxideq2,{dx,dy,dz},{0,1,0}))
% fprintf(' = ')
% disp(subs(ddxideq,{dx,dy,dz},{0,1,0}))
% 
% 
% fprintf('\n\n')
% ddetadeq = det([Tinv(:,1),ddx])/Jinv;
% ddetadeq2 = detadx*dx+detady*dy+detadz*dz;
% 
% disp(subs(ddetadeq2,{dx,dy,dz},{1,0,0}))
% fprintf(' = ')
% disp(subs(ddetadeq,{dx,dy,dz},{1,0,0}))
% 
% disp(subs(ddetadeq2,{dx,dy,dz},{0,1,0}))
% fprintf(' = ')
% disp(subs(ddetadeq,{dx,dy,dz},{0,1,0}))
% 
% 
% 
% fprintf('\n\nJ = ')
% disp(det(T))
% fprintf('\nJinv = ')
% disp(det(Tinv))
% 
% 
% 
