function [vertex, face, cell, vertex_cv, face_cv, cell_cv, cell_all, all_to_cv] = generate_mesh(imax, jmax, grid_type, neq, vertex_centered, r)

perturb = 0.0;
seed = 0;

if (grid_type==1)
%     [ti, xi, icell] = generate_equilateral_mesh(2);
    load('circle_mesh-201.mat')
    xi=x;
    
    
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
    xi(:,1)=(xi(:,1)-minx)./range_x;
    
    miny = min(xi(:,2));
    range_y = max(xi(:,2)) - miny;
    xi(:,2) = (xi(:,2)-miny)./range_x;
    
    ncells = size(xi,1);

    
    
    
    %% Build unstructured grid ----------------------------------------------------
%     [vertex, cell, face] = compute_grid_derived_data(xi(:,1),xi(:,2), grid_type, ti);
    
    [ti, vertex] = delaunay_triangulation( xi(:,1),xi(:,2) );
    vertex = [vertex,zeros(size(vertex(:,1)))];
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
    ncells = size(vertex,1);
    xi = vertex;
    lxi = ncells;
    
    
    for i = 1:size(cell.nodes,1)
       n = cell.nodes(i,1);
       xcc(i,1) =  mean( xi(cell.nodes(i,2:n+1),1) );
       xcc(i,2) =  mean( xi(cell.nodes(i,2:n+1),2) );
       xcc(i,3) = 0;
       xcc_node(i,1) = lxi+i;
    end
    xi = [xi; xcc];
    lxi = size(xi,1);
    
    for i = 1:length(face)
       xfc(i,1) =  mean( vertex(face(i).nodes,1) );
       xfc(i,2) =  mean( vertex(face(i).nodes,2) );
       xfc(i,3) = 0;
       xfc_node(i,1) = lxi+i;
    end
    xi = [xi; xfc];
    
    cnt = 1;
    lvertex = size(vertex,1);
    for i = 1:size(cell.nodes,1)
        nn = cell.nodes(i,1);
        nodes = cell.nodes(i,2:nn+1);
        ff = cell.faces(i,:);
        
        for j = 1:length(ff)
          fnodes(j,:) = face(ff(j)).nodes;
        end
        
        lxc = xcc_node(i);
        
        nodes = [nodes,nodes(1)];
       for j = 1:nn
           [~,~,Iface] = intersect( [nodes(j),nodes(j+1);nodes(j+1),nodes(j)],...
                                   fnodes, 'rows' );
           
           lfc = xfc_node(ff(Iface));
                               
           ti_tri(cnt,:) = [nodes(j), lfc, lxc];
           cnt = cnt + 1;
           
           ti_tri(cnt,:) = [lxc, lfc, nodes(j+1)];
           cnt = cnt + 1;
       end
     
       
    end
    ti = ti_tri;
    

   
    
%     xvc = [vertex(:,1:2); xcc; xfc];
%     [ti,xi] = delaunay_triangulation(xvc(:,1),xvc(:,2));

    % Remove the faces connected to the original nodes
    % Extract faces

%     [face_nodes, cell_faces, face_cells] = find_neighbors(ti);
%     nxfc = size(xfc,1);
%     n = size(vertex,1)+size(xcc,1)+1;
%     for i = n : size(xvc,1)
%       face_connected_to_face = sort([i*ones(nxfc-1,1), [n:i-1, i+1:size(xvc,1)]'],2,'ascend');
%       [f, IA, IB] = intersect(face_connected_to_face,face_nodes,'rows');
%       for j = 1:length(IB)
%          [ti, face_nodes, face_cells] = delaunay_face_swap(face_cells(IB(j),1), face_cells(IB(j),2), IB(j), ti, face_nodes, face_cells, xi(:,1:2));
%       end
% 
%     end
    
    ti = order_nodes(ti,xi);
%     [vertex_cv, cell_all, face_all] = compute_grid_derived_data(xi(:,1),xi(:,2), grid_type, ti);
    [face_all, cell_all] =  compute_face_data_from_triangulation(ti,xi);
    vertex_cv = xi;
    
    % Now that all cells are connected to the original nodes, remove the
    % nodes from the cell connectivity

    for i = 1:ncells
        cell_nc.volume(i) = 0;


        I = find(ti==i);
        for j = 1:length(I)
            cell_nbr(j) = rem(I(j),size(ti,1));
            if (cell_nbr(j)==0); cell_nbr(j)=size(ti,1); end
            cell_nc.volume(i) = cell_nc.volume(i) + abs(cell_all.volume(cell_nbr(j)));
            all_to_cv(i,j) = cell_nbr(j);
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
                face_nc(f(jj)).normal = normal;
                face_nc(f(jj)).area   = A;
                
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


    % Check order of ti_nc
    for i = 1:size(ti_nc,1)
        cn = ti_nc(i,:);
        Icn = find(cn~=0);
        lcn = length(Icn);
        xc = mean(vertex_cv(cn(Icn),1:2),1);
        A = compute_area(xc(1), xc(2), ...
                         vertex_cv(cn(1),1), vertex_cv(cn(1),2), ...
                         vertex_cv(cn(2),1), vertex_cv(cn(2),2) );
        if A<0
            ti_nc(i,1:lcn) = fliplr(ti_nc(i,1:lcn));
        end
        
    end
    
     [face_nodes, cell_faces, face_cells] = find_neighbors(ti_nc);
     cell_cv.volume = cell_nc.volume;
     cell_cv.cell_type = cell_nc.cell_type;
     cell_cv.faces = cell_faces;
     for i = 1:size(cell_faces,1)
         I = find(cell_faces(i,:)~=0);
         cell_cv.nface(i) = length(I);
     end
     
     

     for i = 1:length(face_cells)
         face_cv(i).nodes = face_nodes(i,:);
         fn = face_nodes(i,:);
         
         if (face_cells(i,2)==0);
             face_cells(i,2) = face_cells(i,1);
             face_cells(i,1) = -1;
            
             %Check to see if neg/pos nodes are correct
             cp = face_cells(i,2);
             cn = ti_nc(cp,:);
             Icn = find(cn~=0);
             cn = [cn(Icn),cn(1)];
             In = find(fn(1) == cn,1);
             if (cn(In+1) ~= fn(2))
                 face_nodes(i,:) = fliplr(face_nodes(i,:));
             end
             
         else
             %Check to see if neg/pos nodes are correct
             cp = face_cells(i,2);
             cn = ti_nc(cp,:);
             Icn = find(cn~=0);
             cn = [cn(Icn),cn(1)];
             In = find(fn(1) == cn,1);
             if (cn(In+1) ~= fn(2))
                 face_cells(i,:) = fliplr(face_cells(i,:));
             end
         end
         
         face_cv(i).cell_plus = face_cells(i,2);
         face_cv(i).cell_neg  = face_cells(i,1);
         
         ds_vec = vertex_cv(face_nodes(i,2),:) - vertex_cv(face_nodes(i,1),:);
         face_cv(i).area = sqrt( sum(ds_vec.^2) );
         face_cv(i).normal(1:2) = [ ds_vec(2), -ds_vec(1) ]/face_cv(i).area;
     end
     
%     cell_cv = cell_nc;
%     face_cv = face_nc;
    
    cell_cv.nodes(:,2:size(ti_nc,2)+1) = ti_nc;
    for i = 1:size(ti_nc,1);
        I = find(ti_nc(i,:)~=0);
        cell_cv.nodes(i,1) = length(I);
        cell_cv.xc(i,:) = mean(vertex_cv(ti_nc(i,I),1:2),1);
    end
    cell_cv.ncells = length(cell_cv.volume); %Find number of control volumes
    cell_cv.vtk_size = sum(cell_cv.nodes(:,1)) + cell_cv.ncells;



else
    vertex_cv = vertex;
    cell_cv = cell;
    face_cv = face;
    cell_all = cell;
    all_to_cv = [1:length(cell.volume)]';
    
end


for j = 1:length(face_cv)
   cell_nbrs(j,1) = face_cv(j).cell_plus;
   cell_nbrs(j,2) = face_cv(j).cell_neg;
end

for i = 1:cell_cv.ncells
    nnbr = 0;
    nbrs_pos = [];
    I = find(cell_nbrs==i);
    for j = 1:length(I)
        icell = rem(I(j),length(face_cv)); if icell == 0; icell=length(face_cv); end
        % Then left column of cell_nbrs
        nbrs_pos = [nbrs_pos; cell_nbrs(icell,:)];
    end
    
    Inbrs = find( nbrs_pos ~= i & nbrs_pos~=-1 );
    nbr = unique(nbrs_pos(Inbrs));
    cell_cv.nnbr(i) = length(nbr);
    cell_cv.nbrs(i,1:cell_cv.nnbr(i)) = nbr;
    
    
end

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

