function [cell, face, vertex, icell] = construct_te_smooth_stencil(...
                                       cell_gbl, face_gbl, vertex_gbl, ...
                                       kexact_order, kexact_type,...
                                       fit_type, flux_integral_order, n )


nstencil = kexact_order + 1;
ilow  =   - nstencil;
ihigh = 1 + nstencil;

index = ilow:ihigh;
[eta, xi] = meshgrid(index, index);

for j = 1:length(index)
  for i = 1:length(index)
    x(i,j) = cell_gbl.map(n).x([xi(i,j),eta(i,j)]);
    y(i,j) = cell_gbl.map(n).y([xi(i,j),eta(i,j)]);
  end
end

if cell_gbl.cell_type(n) == 0
  grid_type = 0;
  imax = size(x,1);
  jmax = size(x,2);

  icell = nstencil + 1;
  nunknowns = (kexact_order+1)^2;

elseif cell_gbl.cell_type(n) == 1
  grid_type = 1;
  fprintf('Still need to do this one!\n')
end

[vertex, cell, face] = compute_grid_derived_data(x,y, grid_type);

setup_mapping;
setup_reconstruction;


ij_to_vec = @(i,j) (j-1)*(imax-1) + i; %i,j location to cell vector location


for j = 1:jmax-1
  jl = j-floor((kexact_order+1)/2);
  jh = jl + (kexact_order+1);

  if (jl < 1)
    jl = 1;
    jh = jl + (kexact_order+1);
  elseif (jh > imax - 1)
    jh = jmax-1;
    jl = jh - (kexact_order+1);
  end


  for i = 1:imax-1
    il = i-floor((kexact_order+1)/2);
    ih = il + (kexact_order+1);

    if (il < 1)
      il = 1;
      ih = il + (kexact_order+1);
    elseif (ih > imax - 1)
      ih = imax-1;
      il = ih - (kexact_order+1);
    end


    cnt = 1;
    for jj = jl:jh
      for ii = il:ih
        cells(cnt) = ij_to_vec(ii,jj);
        cnt = cnt + 1;
      end
    end

    cell.stencil(ij_to_vec(i,j)).cells = cells;

  end
end

cell.nunknowns = nunknowns;


