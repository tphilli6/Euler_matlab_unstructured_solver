function exact = compute_analytic_exact(cell, vertex, exact_order, analytic_soln, cell_tri_to_cv)


setup_mapping;

exact_tri = analytic_solution(vertex, cell, exact_order, analytic_soln);

exact = zeros(size(cell_tri_to_cv,1),length(exact_tri(1,:)));
for i = 1:size(cell_tri_to_cv,1)
    tri_to_cv = cell_tri_to_cv(i,:);
    I = find(tri_to_cv~=0);
    tri_to_cv = tri_to_cv(I);
 
    vol = 0;
   for j = 1:length(tri_to_cv)
       exact(i,:) = exact(i,:) + exact_tri(tri_to_cv(j),:)*cell.volume(tri_to_cv(j));      
       vol = vol + cell.volume(tri_to_cv(j));
   end
   exact(i,:) = exact(i,:)./vol;
end