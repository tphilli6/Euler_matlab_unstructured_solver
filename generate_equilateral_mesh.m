function [ti, xi, c_cntr] = generate_equilateral_mesh(nrows, write_mesh_file)
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

if (nargin==1); write_mesh_file=0; end
% 
%     
% % nrows = 2; % Number of rows around the cell and face neighbors
% nx = 2*nrows + 4;
% ny = nx-1;
% ic = (ny-1)/2 + 1;
% dx = 1;
% dy = 0.5*tand(60);
% 
% xivec = [];
% etavec = [];
% c_cntr = 0;
% for i = 1:ny
%    irow = nx - abs(i-ic);
%    xit = [0:dx:(irow-1)*dx] + dx/2*abs(i-ic);
%    etat = (i-1)*dy*ones(size(xit));
%    
%    xivec = [ xivec, xit];
%    etavec = [ etavec, etat];
%    if i-ic < 0
%        c_cntr = c_cntr + (nx - abs(i-ic)-1) + (nx - abs(i-ic+1)-1);
%    elseif i-ic == 1
%        c_cntr = c_cntr + ((nx-1)+(nx-2)-1)/2;
%    end
%    
% end
% 
% 
% [ti,xi] = delaunay_triangulation(xivec,etavec);
% 
% xc = mean(xi(ti(c_cntr,:),1));
% yc = mean(xi(ti(c_cntr,:),2));
% xi(:,1) = xi(:,1)-xc;
% xi(:,2) = xi(:,2)-yc;
% 
% if which_eq_mesh==1; return; end

ndivisions = 6;

neigh = [0:nrows+1]*ndivisions;

x(1,:) = [0,0];
xp(1) = 0;
cnt = size(x,1);
ti = [];
faces = [];
low = 2;
high = low + ndivisions-1;
for i = low:high;
   x(i,:) = [cos( (i-2)*(2*pi/ndivisions) ), sin( (i-2)*(2*pi/ndivisions) )];
   xp(i) = 1;
   
   faces = [faces; xp(i), i];
   if i==high
    ti = [ti; xp(i), i, low];
    faces = [faces; i,low];
   else
    ti = [ti; xp(i), i, i+1];
    faces = [faces; i, i+1];
   end
end
length(faces)

cnt = size(x,1);
cell_cnt = 1;
for n = 2:length(neigh)-1
    low = length(xp)-neigh(n)+1;
    high = length(xp);
    skip = 0;
    face_bndry = [];
    
    for i = low:high
        dx = x(i,1)-x(xp(i),1);
        dy = x(i,2)-x(xp(i),2);
        th = atan( dy/dx );
        if dx<0;
            th=th+pi;
        end

        if skip == 0
            for j = 1:2;
              x(cnt+1,:) = [ cos( th + (j-1)*(2*pi/ndivisions) ), sin( th + (j-1)*(2*pi/ndivisions) ) ] + x(i,:);
              xp(cnt+1) = i;
              cnt = cnt + 1;
            end
            
            if i>low
                ti = [ti; cnt-2, cnt-1, xp(cnt-1)];
                faces = [faces; cnt-2, cnt-1];
                face_bndry = [face_bndry; cnt-2, cnt-1];
            end
            ti = [ti; xp(cnt-1), cnt-1, cnt];
            faces = [faces; xp(cnt-1), cnt-1];
            faces = [faces; cnt-1, cnt];
            face_bndry = [face_bndry; cnt-1, cnt];
            
            faces = [faces; xp(cnt), cnt];
            if i<high
                ti = [ti; xp(cnt), cnt, xp(cnt)+1];
                faces = [faces; cnt, xp(cnt)+1];
            else 
                ti = [ti; xp(cnt), cnt, low];
                ti = [ti; cnt, cnt - neigh(n+1)+1, xp(cnt - neigh(n+1)+1) ];
                faces = [faces; cnt, low];
                faces = [faces; cnt, cnt - neigh(n+1)+1];
                face_bndry = [face_bndry; cnt, cnt - neigh(n+1)+1];
            end
%             plot_cells(ti,x);
%             if i==high
%                 ti = [ti; xp(i+1), i+1, low];
%             else
%                 ti = [ti; xp(i+1), i+1, i+2];
%             end
            
            if skip == n-2
                skip = 0;
            else
                skip = skip + 1;
            end
            
        else
              x(cnt+1,:) = [ cos( th ), sin( th ) ] + x(i,:);
              xp(cnt+1) = i;
              cnt = cnt + 1;
              
              ti = [ti; xp(cnt), cnt-1, cnt];
%               faces = [faces; xp(cnt), cnt-1];
              faces = [faces; cnt-1, cnt];
              faces = [faces; xp(cnt), cnt];
              face_bndry = [face_bndry; cnt-1, cnt];

              if i<high
                ti = [ti; xp(cnt), cnt, xp(cnt)+1];
                faces = [faces; cnt, xp(cnt)+1];
              else
                ti = [ti; cnt - neigh(n+1), cnt, low];
                ti = [ti; cnt, cnt - neigh(n+1)+1, xp(cnt - neigh(n+1)+1)];
                faces = [faces; cnt, cnt - neigh(n+1)+1];
                faces = [faces; cnt, xp(cnt - neigh(n+1)+1)];
                face_bndry = [face_bndry; cnt, cnt-neigh(n+1)+1];
              end
%             plot_cells(ti,x);
              if skip == n-2
                  skip = 0;
              else
                  skip = skip + 1;
              end


        end

    end
    
              
%                   [face_nodes, cell_faces,face_cells] = find_neighbors(ti);
%             for j = 1:size(faces,1)
%                 if (j==1); eval(['hold ', 'off']); end
%                 plot( x(faces(j,:),1), x(faces(j,:),2) ,'r-o','LineWidth',2);
%                 hold on
%             end
            
%             for j = 1:size(face_bndry,1)
%                 if (j==1); eval(['hold ', 'on']); end
%                 plot( x(face_bndry(j,:),1), x(face_bndry(j,:),2) ,'b-*','LineWidth',2);
%                 hold on
%             end

end

% Find neighboring cells
face_cell = zeros(size(faces));
for i = 1:length(faces);
    
    for n = 1:length(ti)
        nodes = [ti(n,:), ti(n,1)];
        if any( all(nodes(1:2)==faces(i,:)) || all(nodes(2:3)==faces(i,:)) || all(nodes(3:4)==faces(i,:))  )
            face_cell(i,1) = n;
        elseif any( all(nodes([2,1])==faces(i,:)) || all(nodes([3,2])==faces(i,:)) || all(nodes([4,3])==faces(i,:))  )
            face_cell(i,2) = n;
        end

    end
    
    if face_cell(i,1)==0 && face_cell(i,2)~=0 
        face_cell(i,:) = face_cell(i,[2,1]);
        face_cell(i,2) = -1;
    elseif face_cell(i,2)==0 && face_cell(i,1)~=0
        face_cell(i,2) = -1;
    elseif face_cell(i,1)>0 && face_cell(i,2)>0
        
    else
        fprintf('Error! Face %4.0f not found in triangulation\n',i);
    end
        
end

for i = 1:size(ti,1)
    II(i) = length( find(face_cell==i) );
end
I = find(II~=3);


% Find neighboring cells for bndry
face_bndrycell = zeros(size(face_bndry));
for i = 1:length(face_bndry);
    
    for n = 1:length(ti)
        nodes = [ti(n,:), ti(n,1)];
        if any( all(nodes(1:2)==face_bndry(i,:)) || all(nodes(2:3)==face_bndry(i,:)) || all(nodes(3:4)==face_bndry(i,:))  )
            face_bndrycell(i,1) = n;
        elseif any( all(nodes([2,1])==face_bndry(i,:)) || all(nodes([2,1])==face_bndry(i,:)) || all(nodes([4,3])==face_bndry(i,:))  )
            face_bndrycell(i,2) = n;
        end

    end
    
    if face_bndrycell(i,1)==0 && face_bndrycell(i,2)~=0 
        face_bndrycell(i,:) = face_bndrycell(i,[2,1]);
        face_bndrycell(i,2) = 1;
    elseif face_bndrycell(i,2)==0 && face_bndrycell(i,1)~=0
        face_bndrycell(i,2) = 1;
    elseif face_bndrycell(i,1)>0 && face_bndrycell(i,2)>0
        
    else
        fprintf('Error! Face %4.0f not found in triangulation\n',i);
    end
        
end

xi = x/(nrows+1); % scale to [-1,1] for x. y is also scale appropriately.
c_cntr = 1;


if write_mesh_file
%% Write *.mesh file

ncells = size(ti,1);
nfaces = size(faces,1);
nbndryfaces = size(face_bndry,1);
nvertex = size(x,1);

fid = fopen('equilateral.mesh','w');

% header: ncells, nfaces, nboundaryfaces, nvertex
fprintf(fid, '%7.0f %7.0f %7.0f %7.0f\n', ncells, nfaces, nbndryfaces, nvertex);

% vertices
for i = 1:nvertex
    fprintf(fid, '%23.15e %23.15e\n', xi(i,:));
end

% cell 1, cell 2, face node 1, face node 2
for i = 1:nfaces
    fprintf(fid, '%5.0f %5.0f %5.0f %5.0f\n', face_cell(i,:)-1, faces(i,:)-1);
end

% cell 1, bndry id, face node 1, face node 2
for i = 1:nbndryfaces
    fprintf(fid, '%5.0f %5.0f %5.0f %5.0f\n', face_bndrycell(i,:)-[1,0], face_bndry(i,:)-1);
end

% cell type ? =1
for i = 1:ncells
    fprintf(fid, '1\n');
end

fclose(fid);
end


