% This script builds the stencil for the kexact reconstruction and stores them in
% cell.stencil(i).cells(1:nstencil)
% cell.nunkowns

% Reconstruction types based on what's set in equation_setup.m
if strcmp(fit_type, 'lsq')
  nstencil = cell.nunknowns + 1;
 
elseif strcmp(fit_type, 'kexact')
  nstencil = cell.nunknowns;

end
nstencil

for n = 1:cell.ncells
 
  %cells(1:nstencil)=0
  stn_cnt= 1;
  cells(stn_cnt) = n; %assign current cell to stencil
  cell_dup = [1];
  nbrs_list = [n];
  while (stn_cnt<nstencil)

    for i = 1:stn_cnt;
      cur_cell = cells(i);
      nnbr = cell.nnbr(cur_cell);
      
      for j = 1:nnbr
%        [i,j,cell.nbrs(cur_cell,j), nbrs_list]
        if any(cell.nbrs(cur_cell,j)==nbrs_list)
          I= find(cell.nbrs(cur_cell,j)==nbrs_list);
          cell_dup(I) = cell_dup(I) + 1;
        else
          nbrs_list = [nbrs_list, cell.nbrs(cur_cell, j)];
          cell_dup  = [cell_dup, 1];
          stn_cnt = stn_cnt + 1;
        end
      end

    end

    [sorted, I] = sort(cell_dup,'descend');

    cells = nbrs_list( I(1:min(stn_cnt,nstencil)) );

  end
  cell.stencil(n).cells = cells;
  

end


% Write out the stencil in vtk format for inspection
% Need:
%       cell.ncells
%       cell.vtk_size
%       cell.nodes
%       cell.cell_type
for n = 1:cell.ncells

  cells = cell.stencil(n).cells;
  v_temp = [];
  cell_temp.ncells = length(cells);

  vrtx_list = [];
  vrtx_cnt  = 0;
  vertex_temp = [0,0,0];
  for i = 1:length(cells)
  
    cell_temp.cell_type(i) = cell.cell_type(cells(i)); %store cell type
    nvertex = cell.nodes(cells(i),1);%store nvertex
    vrtx = cell.nodes(cells(i),2:nvertex+1); %pull of vertex for cell

    % Build vertex list for new cells
    if isempty(vrtx_list) % if there is no current vertex list create one
      vrtx_list = vrtx; % list of absolute indicies
      vertex_temp(1:nvertex,:) = vertex(vrtx_list,:); % list of node locations
      vrtx_cnt = length(vrtx); % current number of nodes
      cell_temp.nodes = [nvertex,1:nvertex]; % store the nodes

    else %else check if the vertex has already been used, if not concatinate it to the list
      vrtx_out_cnt = 1;
      cell_temp.nodes(i,1) = nvertex;
      for j = 1:length(vrtx)

        I = find(vrtx(j) == vrtx_list);
        if length(I)>0
          vrtx_out(1,vrtx_out_cnt) = I(1);
          vrtx_out_cnt = vrtx_out_cnt + 1;

        else
          vrtx_list = [vrtx_list, vrtx(j)];
          vrtx_cnt = vrtx_cnt + 1;

          vrtx_out(1,vrtx_out_cnt) = vrtx_cnt;
          vrtx_out_cnt = vrtx_out_cnt + 1;
          
          vertex_temp = [vertex_temp; vertex(vrtx(j),:)];

        end
      end
      cell_temp.nodes(i,:) = [nvertex, vrtx_out];

    end

  end

  cell_temp.vtk_size = sum(cell_temp.nodes(:,1)) + cell_temp.ncells;

%write stencil to file
  fid = fopen(['grid-stencil-',num2str(n),'.vtk'],'w');
  write_vtk(vertex_temp, cell_temp, fid);
  fclose(fid);

end

