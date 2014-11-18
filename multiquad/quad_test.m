%% quadrature test
dim=2;

for i = 1:20
    order=i;

    [x,w] = sparse_grid(dim, order);
    
    nint_exactness(x,w,order+1);
    
end