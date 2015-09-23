function [vertex, face, cell, vertex_cv, face_cv, cell_cv, cell_all, all_to_cv] = generate_mesh(imax, jmax, grid_type, neq, vertex_centered, r, grid_in)

perturb = 0.;
seed = 0;
rng(seed);

if (grid_type==1)
    [ti, xi, icell] = generate_equilateral_mesh(3);
%     load('circle_mesh-201.mat')
%     xi=x;


    for n=1:r-1
        [face_nodes] = find_neighbors(ti);
%         for j = 1:size(ti,1)
%             xi_c(j,1) = mean( xi(ti(j,:),1) );
%             xi_c(j,2) = mean( xi(ti(j,:),2) );
%         end
        for j = 1:size(face_nodes,1)
            xi_fc(j,1) = mean( xi(face_nodes(j,:),1) );
            xi_fc(j,2) = mean( xi(face_nodes(j,:),2) );
        end
        xi_new = [xi;xi_fc];
        [ti, xi] = delaunay_triangulation( xi_new(:,1),xi_new(:,2) );
    end
    minx = min(xi(:,1));
    range_x = max(xi(:,1))-minx;
    xi(:,1)=(xi(:,1)-minx)./range_x*0.5;
    
    miny = min(xi(:,2));
    range_y = max(xi(:,2)) - miny;
    xi(:,2) = (xi(:,2)-miny)./range_x*0.5;
    
    grid_center = mean(xi);
    xi(:,1) = xi(:,1) - grid_center(1);
    xi(:,2) = xi(:,2) - grid_center(2);
    
    ncells = size(xi,1);

    
    
    
    %% Build unstructured grid ----------------------------------------------------
%     [vertex, cell, face] = compute_grid_derived_data(xi(:,1),xi(:,2), grid_type, ti);
    
    [ti, vertex] = delaunay_triangulation( xi(:,1),xi(:,2) );
    vertex = [vertex,zeros(size(vertex(:,1)))];
    [face, cell] =  compute_face_data_from_triangulation(ti,vertex);
    
    % Find boundary nodes
    cnt = 1;
    for i = 1:numel(face)
        if ( face(i).cell_neg < 0 )
            fbndry(cnt,:) = face(i).nodes;
            cnt = cnt + 1;
        end
    end
    Ibndry = unique(fbndry);
    
    I = [1:size(vertex,1)]';
    Iinterior = setdiff( I, Ibndry);

    % Perturb interior nodes only
    [face_nodes] = find_neighbors(ti);
    fxi1(:,1:2) = xi( face_nodes(:,1),1:2);
    fxi2(:,1:2) = xi( face_nodes(:,2),1:2);
    ds = sqrt( (fxi2(:,1)-fxi1(:,1)).^2 + (fxi2(:,2)-fxi1(:,1)).^2  );
    dsmax = min(ds)*perturb;
    
    ninterior = length(Iinterior);
    del = dsmax*rand(ninterior,1);
    dth = 2*pi*rand(ninterior,1);
    vertex(Iinterior,:) = vertex(Iinterior,:) + [del.*cos(dth), del.*sin(dth),zeros(size(dth))];
    
elseif grid_type==-1
    
    [ti, xi] = read_mesh_file(grid_in);
%     minx = min(xi(:,1));
%     range_x = max(xi(:,1))-minx;
%     vertex(:,1)=(xi(:,1)-minx)./range_x*0.5;
%     
%     miny = min(xi(:,2));
%     range_y = max(xi(:,2)) - miny;
%     vertex(:,2) = (xi(:,2)-miny)./range_x*0.5;
%     
%     vertex(:,3) = 0;
    vertex = xi;
    [face, cell] =  compute_face_data_from_triangulation(ti,vertex);
    
elseif grid_type==-2
    mu = complex(-0.0475,0);
    te_angle = -2*(0.594689181*(1/2*0.298222773/sqrt(1) - 0.127125232 - 2*0.357907906*1 + 3*0.291984971*1*1 - 4*0.105174606*1*1*1));
        
    [ti, xi] = generate_kt_airfoil(mu, te_angle);

    vertex = xi;
    vertex(:,3) = 0;
    [face, cell] =  compute_face_data_from_triangulation(ti,vertex);
    
    
else
    


    x0 = 0;
    xL = 1;
    y0 = 0;
    yL = 1;



    % Generate grid points  -------------------------------------------------------
    x = zeros(imax, jmax);
    y = zeros(imax, jmax);
    % Cartesian
    xi = linspace(x0,xL,imax);
    eta = linspace(y0,yL,jmax);

    for n = 1:imax
      x(n,:) = xi(n);
      y(n,:) = eta;
    end

    rng(seed);
    dx = 1/(imax-1)*perturb;
    dy = 1/(jmax-1)*perturb;
    xdiff = 2*(rand( imax-2, jmax-2)-0.5)*dx;
    ydiff = 2*(rand( imax-2, jmax-2)-0.5)*dy;
    x(2:end-1,2:end-1) = x(2:end-1,2:end-1) + xdiff;
    y(2:end-1,2:end-1) = y(2:end-1,2:end-1) + ydiff;



    % Nozzle
    %eta = linspace(0,2*pi,jmax);
    %xi = linspace(x0,xL,imax);
    %
    %for n = 1:imax
    %  x(n,:) = xi(n);
    %  ymax = 1 + 0.5*cos(xi(n));
    %  y(n,:) = linspace(-ymax,ymax,jmax);
    %end



    %% Build unstructured grid ----------------------------------------------------
    [vertex, cell, face] = compute_grid_derived_data(x,y, grid_type);
    

    
end



%% Vertex Centered Setup
% If vertex centered compute cell centers and face centers to add to
% the grid afterwords the vertex is removed from the cell and face
% trianglation
if vertex_centered

  [cell_cv, face_cv, vertex_cv, cell_all, all_to_cv] = compute_vertex_centered_nonsense(cell, face, vertex);

else
    vertex_cv = vertex;
    cell_cv = cell;
    face_cv = face;
    cell_all = cell;
    all_to_cv = [1:length(cell.volume)]';
    
end


cell_cv = find_neighbor_cells( cell_cv, face_cv );

% Loop over cells
% for i = 1:cell_cv.ncells
%    % Loop over faces
%    nbrs = [];
%    nnbr = 0;
%    for j = 1:size(face_cv,2)
%        
%        % If there are two faces then find the one that's not the current
%        % cell
%        I = find(cell_cv.faces==j);
%        if length(I)==2
%            icell1 = rem(I(1),cell_cv.ncells);if icell1 == 0; icell1=cell_cv.ncells; end;
%            icell2 = rem(I(2),cell_cv.ncells);if icell2 == 0; icell2=cell_cv.ncells; end;
%            
%            if icell1 == i || icell2 == i
%                % Loop over the two faces
%                for ii = 1:length(I)
%                    icell = rem(I(ii),cell_cv.ncells); % Translate the linear index into a row index
%                    if icell == 0; icell=cell_cv.ncells; end;
% 
%                    %If the row index corresponds to not the current cell then
%                    %increment nnbr and add to the nbrs cell list
%                    if icell~=i
%                        nnbr = nnbr + 1;
%                        nbrs(1:nnbr) = [nbrs(1:nnbr-1), icell];
%                    end
%                end
%            end
%        end
%        
%    end
%    [C] = unique(nbrs);
%    cell_cv2.nnbr(i) = length(C);
%    cell_cv2.nbrs(i,1:cell_cv.nnbr(i)) = C;
%    cell_cv2.nbrs(i,1:cell_cv.nnbr(i))
%    cell_cv.nbrs(i,:)
% end


%% Write grid to file %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fid = fopen('grid_cv.vtk','w');
write_vtk_cv(vertex_cv, cell_cv, fid);
fclose(fid);
    
fid = fopen('grid.vtk','w');
write_vtk(vertex, cell, fid);
fclose(fid);

fid = fopen('grid.vtk','a');
if vertex_centered
  fprintf(fid,'CELL_DATA %8.0f\n',length(cell.volume));
  vol_ave = vertex_to_cell_average(cell_cv.volume',cell);
  write_vtk_data(vol_ave, 'volume', fid);
else
  fprintf(fid,'CELL_DATA %8.0f\n',length(cell.volume));
  write_vtk_data(cell.volume, 'volume', fid);
end
size(cell_cv.volume)
fclose(fid);

for n = 1:length(face)
  face(n).ul = zeros(1,neq);
  face(n).ur = zeros(1,neq);
end

