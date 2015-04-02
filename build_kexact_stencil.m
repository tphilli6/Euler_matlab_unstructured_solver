% This script builds the stencil for the kexact reconstruction and stores them in
% cell.stencil(i).cells(1:nstencil)
% cell.nunkowns

% Print status message
fprintf(['Building stencil for ',fit_type,' reconstruction\n'])

% Reconstruction types based on what's set in equation_setup.m
if strcmp(fit_type(1:3), 'lsq')
  extra = floor(cell.nunknowns*1.5);
  nstencil = cell.nunknowns + extra;
  nstencil = 13; % HARDCODE FOR TESTING
elseif strcmp(fit_type, 'kexact')
  nstencil = cell.nunknowns;
elseif strcmp(fit_type, 'extended')
  % requires manual setting of extra_terms prior to calling
  nstencil = cell.nunknowns + extra_terms;
end


for nn = 1:cell.ncells

  %cells(1:nstencil)=0
  stn_cnt= 1;
  cells(stn_cnt) = nn; %assign current cell to stencil
  cell_dup = [2];
  nbrs_list = [nn];

  while (stn_cnt<nstencil)

    for i = 1:stn_cnt;
      cur_cell = cells(i);
      nnbr = cell.nnbr(cur_cell);
      
      for j = 1:nnbr
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

    [~, I] = sort(cell_dup,'descend');
    cells = nbrs_list( I );

  end

    % Add on that chooses stencil based on condition number of reconstruction
    if (length(cell_dup)>nstencil)
      [sorted, I] = sort(cell_dup,'descend');

      matching_edges = find(sorted(nstencil)==sorted(2:end))+1;

      needed_matching = find(sorted(nstencil)==sorted(2:nstencil));
      n_needed = length(needed_matching);

      poss_sten_ind = nchoosek(matching_edges, n_needed);
      poss_sten = cells(poss_sten_ind);

      % Annoyingly, if n_needed is only one the poss_sten is a row vector otherwise it's an array
      if (n_needed == 1)
        poss_sten = poss_sten';
      end

      sten_base = cells(1:nstencil-n_needed);

      cond_check = 0;
      for nsten = 1:size(poss_sten,1)
        sten_check = [sten_base, poss_sten(nsten,:)];
        [~, cond_check(nsten)] = reconstruction_lhs(sten_check, ...
                                      cell.reconstruction,...
                                      fit_type,1);

      end
      [~,sten_choose] = max(cond_check);
      cells = [sten_base, poss_sten(sten_choose,:)];

    end


  cell.stencil(nn).cells = cells(1:nstencil);

end

% Write out the stencil in vtk format for inspection
% Need:
%       cell.ncells
%       cell.vtk_size
%       cell.nodes
%       cell.cell_type
if any(strcmp(who,'write_sten'))
    
    for nn = 1:cell.ncells

      cells = cell.stencil(nn).cells;
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
      fid = fopen(['grid-stencil-',num2str(nn),'.vtk'],'w');
      write_vtk(vertex_temp, cell_temp, fid);
      fclose(fid);

    end

end

