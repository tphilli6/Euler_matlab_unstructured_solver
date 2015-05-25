function [cell_cv, face_cv, vertex_cv, cell_all, all_to_cv] = compute_vertex_centered_nonsense(cell, face, vertex)



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
