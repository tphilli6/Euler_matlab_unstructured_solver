function cell = find_neighbor_cells(cell, face)

for j = 1:length(face)
   cell_nbrs(j,1) = face(j).cell_plus;
   cell_nbrs(j,2) = face(j).cell_neg;
end

for i = 1:cell.ncells
    nnbr = 0;
    nbrs_pos = [];
    I = find(cell_nbrs==i);
    for j = 1:length(I)
        icell = rem(I(j),length(face)); if icell == 0; icell=length(face); end
        % Then left column of cell_nbrs
        nbrs_pos = [nbrs_pos; cell_nbrs(icell,:)];
    end
    
    Inbrs = find( nbrs_pos ~= i & nbrs_pos~=-1 );
    nbr = unique(nbrs_pos(Inbrs));
    cell.nnbr(i) = length(nbr);
    cell.nbrs(i,1:cell.nnbr(i)) = nbr;
    
    
end