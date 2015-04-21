function [vertex, face, cell, vertex_cv, face_cv, cell_cv] = generate_mesh(imax, jmax, grid_type, neq, vertex_centered)

perturb = 0.0;
seed = 0;

if (grid_type==1)
    [ti, xi, icell] = generate_equilateral_mesh(3);
    [face_nodes] = find_neighbors(ti);
    for n=1:0
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
    xi(:,1)=(xi(:,1)-minx)./range_x;
    
    miny = min(xi(:,2));
    range_y = max(xi(:,2)) - miny;
    xi(:,2) = (xi(:,2)-miny)./range_x;
    
    ncells = size(xi,1);

    
    
    
    %% Build unstructured grid ----------------------------------------------------
    [vertex, cell, face] = compute_grid_derived_data(xi(:,1),xi(:,2), grid_type, ti);
        
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
    ncells = size(vertex,1);
    
    for i = 1:size(cell.nodes,1)
       n = cell.nodes(i,1);
       xcc(i,1) =  mean( vertex(cell.nodes(i,2:n+1),1) );
       xcc(i,2) =  mean( vertex(cell.nodes(i,2:n+1),2) );
       xcc(i,3) = 0;
    end

    for i = 1:length(face)
       xfc(i,1) =  mean( vertex(face(i).nodes,1) );
       xfc(i,2) =  mean( vertex(face(i).nodes,2) );
       xfc(i,3) = 0;
    end

    xvc = [vertex; xcc; xfc];
    [ti,xi] = delaunay_triangulation(xvc(:,1),xvc(:,2));

    % Remove the faces connected to the original nodes
    % Extract faces

    [face_nodes, cell_faces, face_cells] = find_neighbors(ti);
    nxfc = size(xfc,1);
    n = size(vertex,1)+size(xcc,1)+1;
    for i = n : size(xvc,1)
      face_connected_to_face = sort([i*ones(nxfc-1,1), [n:i-1, i+1:size(xvc,1)]'],2,'ascend');
      [f, IA, IB] = intersect(face_connected_to_face,face_nodes,'rows');
      for j = 1:length(IB)
         [ti, face_nodes, face_cells] = delaunay_face_swap(face_cells(IB(j),1), face_cells(IB(j),2), IB(j), ti, face_nodes, face_cells, xi(:,1:2));
      end

    end

    [vertex_cv, cell_all, face_all] = compute_grid_derived_data(xi(:,1),xi(:,2), grid_type, ti);

    % Now that all cells are connected to the original nodes, remove the
    % nodes from the cell connectivity

    for i = 1:ncells
        cell_nc.volume(i) = 0;


        I = find(ti==i);
        for j = 1:length(I)
            cell_nbr(j) = rem(I(j),size(ti,1));
            if (cell_nbr(j)==0); cell_nbr(j)=size(ti,1); end
            cell_nc.volume(i) = cell_nc.volume(i) + abs(cell_all.volume(cell_nbr(j)));
        end
        nodes = unique(ti(cell_nbr,:))';
        cell_nc.cell_type(i) = 2; % generic polygon
        cell_nc.nface(i) = length(nodes);

        % Update face mapping for new cell indexing
        cnt = 0; 
        fn_vec=[];
        for j = 1:length(I) % Loop over cells sharing node i
            nf = cell_all.nface(cell_nbr(j));
            f = cell_all.faces(cell_nbr(j),:);


            for jj = 1:length(f) % Loop over faces for the shared cells for node i
                fn = face_all(f(jj)).nodes;

                cp = face_all(f(jj)).cell_plus;
                cn = face_all(f(jj)).cell_neg;
                A  = face_all(f(jj)).area;
                normal = face_all(f(jj)).normal;

                face_add = 0; %logical to store addition of face
                if any(fn==i) % This means the face includes the current node
                   if (cn==-1) % means boundary node, add to cell face list
                       cnt = cnt + 1;
                       cell_nc.faces(i,cnt) = f(jj);
                       face_add = 1;
                   end
                else
                    % if does not include cell node, add face to list
                    cnt = cnt + 1;
                    cell_nc.faces(i,cnt) = f(jj);
                    face_add = 1;
                end

                % If face add
                if face_add
                   fn_vec = [fn_vec; fn];
                end



                % Transfer constant face parameters
                face_nc(f(jj)).nodes = fn;
                face_nc(f(jj)).area  = A;
                face_nc(f(jj)).normal = normal;

                % Update the cell plus and cell minus for the node centered
                % cell indexing
                if (cp==cell_nbr(j))
                  face_nc(f(jj)).cell_plus=i; %set cell plus to new node centered cell
                elseif (cn==cell_nbr(j))
                  face_nc(f(jj)).cell_neg=i;
                end

                if (cn==-1)
                  face_nc(f(jj)).cell_neg=-1;
                end

            end

        end

       % loop over and count the number of cell node i is connected to
       % if = 12 then interior node and remove the node from all cells and
       % form a single cell
       % if < 12 then boundary node

       % Set
       % cell.soln  - solution
       % cell.exact - exact solution
       % cell.nnbr - number of cell neighbors
       % cell.nbrs - cell neighbors

       for j = 1:size(fn_vec,1)-1
          [Ifn, lfn] = find(fn_vec(j,1)==fn_vec([1:j-1,j+1:end],1) );
          [Ifnn, lfnn] = find(fn_vec(j,2)==fn_vec([1:j-1,j+1:end],2) );
          if length(Ifnn)>0
              fn_vec(Ifnn+1,:) = fliplr( fn_vec(Ifnn+1,:) );
          end
              [Ifn2, lfn2] = find(fn_vec(j,2) == fn_vec([1:j-1,j+1:end],1) );
              rs1 = [j+1; Ifn2+1];
              rs2 = [Ifn2+1; j+1];
              fn_vec(rs1,:) = fn_vec(rs2,:);

       end

       if ~isempty(fn_vec)
         ti_nc(i,1:size(fn_vec,1)) = fn_vec(:,1)';
       end
    end

    cell_cv = cell_nc;
    face_cv = face_nc;
    cell_cv.nodes(:,2:size(ti_nc,2)+1) = ti_nc;
    for i = 1:size(ti_nc,1); 
        I = find(ti_nc(i,:)~=0);
        cell_cv.nodes(i,1) = length(I);
        cell_cv.xc(i,:) = vertex(i,1:2);
    end
    cell_cv.ncells = length(cell_cv.volume); %Find number of control volumes
    cell_cv.vtk_size = sum(cell_cv.nodes(:,1)) + cell_cv.ncells;

    
else
    vertex_cv = vertex;
    cell_cv = cell;
    face_cv = face;
    
end



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

