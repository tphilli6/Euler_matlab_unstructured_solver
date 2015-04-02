function [ti, xi, c_cntr] = generate_equilateral_mesh(nrows)
% This function generates an equilateral mesh
% The purpose of this function is to generate a mesh which has nrows
% surrounding a central cells. (i.e. a triangle with all node neighbors
% included for the central cell and cell face neighbors would be nrows = 1)
% This function is meant to be used for reconstructions to compute a
% residual for the central cell c_cntr which requires the face neighbors.
%
% The mesh is centered about the c_cntr geometric center and the size
% corresponds to the equilateral cell used for the mapping transformation
% in compute_triangle_jacobian

% nrows = 1; % Number of rows around the cell and face neighbors
nx = 2*nrows + 4;
ny = nx-1;
ic = (ny-1)/2 + 1;
dx = 1;
dy = 0.5*tand(60);

xivec = [];
etavec = [];
c_cntr = 0;
for i = 1:ny
   irow = nx - abs(i-ic);
   xit = [0:dx:(irow-1)*dx] + dx/2*abs(i-ic);
   etat = (i-1)*dy*ones(size(xit));
   
   xivec = [ xivec, xit];
   etavec = [ etavec, etat];
   if i-ic < 0
       c_cntr = c_cntr + (nx - abs(i-ic)-1) + (nx - abs(i-ic+1)-1);
   elseif i-ic == 1
       c_cntr = c_cntr + ((nx-1)+(nx-2)-1)/2;
   end
   
end


[ti,xi] = delaunay_triangulation(xivec,etavec);

xc = mean(xi(ti(c_cntr,:),1));
yc = mean(xi(ti(c_cntr,:),2));
xi(:,1) = xi(:,1)-xc;
xi(:,2) = xi(:,2)-yc;

