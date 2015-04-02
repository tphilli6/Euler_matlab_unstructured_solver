function [vertex, face, cell] = generate_mesh(imax, jmax, grid_type, neq)

perturb = 0.0;
seed = 0;

if (grid_type==1)
    [ti, xi, icell] = generate_equilateral_mesh(3);
    minx = min(xi(:,1));
    range_x = max(xi(:,1))-minx;
    xi(:,1)=(xi(:,1)-minx)./range_x;
    
    miny = min(xi(:,2));
    range_y = max(xi(:,2)) - miny;
    xi(:,2) = (xi(:,2)-miny)./range_x;
    
    ncells = size(xi,1);

    
    
    
    %% Build unstructured grid ----------------------------------------------------
    [vertex, cell, face] = compute_grid_derived_data(xi(:,1),xi(:,2), grid_type, ti);
    
    % If vertex centered compute cell centers and face centers to add to
    % the grid afterwords the vertex is removed from the cell and face
    % trianglation
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
%     [vertex2, cell, face] = compute_grid_derived_data(xi(:,1),xi(:,2), grid_type, ti);
    
    % Remove the faces connected to the original nodes
    % Extract faces
%     face_vec = [];
%     for i = 1:length(face)
%         f = sort(face(i).nodes,2,'ascend');
%         face_vec = [face_vec; f];
%     end
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
    
    % Now that all cells are connected to the original nodes, remove the
    % nodes from the cell connectivity
    for i = 1:ncells
       % loop over and count the number of cell node i is connected to
       % if = 12 then interior node and remove the node from all cells and
       % form a single cell
       % if < 12 then boundary node 
    end
    
    
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



%% Write grid to file %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fid = fopen('grid.vtk','w');
write_vtk(vertex, cell, fid);
fclose(fid);

fid = fopen('grid.vtk','a');
fprintf(fid,'CELL_DATA %8.0f\n',length(cell.volume));
write_vtk_data(cell.volume, 'volume', fid);
fclose(fid);


for n = 1:length(face)
  face(n).ul = zeros(1,neq);
  face(n).ur = zeros(1,neq);
end

