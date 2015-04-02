function [face_out, cell_out] = remove_face(face, cell, fr)
%% This function removes the face given by fr and combines the corresponding cells in cell

cp = face(fr).cell_plus;
cn = face(fr).cell_neg;

%Keep cell_plus and remove cell_neg
ncell = size(cell.nodes,1);
Ikeep = [1:cn-1,cn+1:ncell];
ncp = cell.nodes(cp,1);
cp_nodes = cell.nodes(cp,2:ncp+1);
ncn = cell.nodes(cn,1);
cn_nodes = cell.nodes(cn,2:ncn+1);
