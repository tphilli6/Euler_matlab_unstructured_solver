function cnt = is_delaunay(ti,x)

% Check if Delaunay
cnt = 0;
for i = 1:size(ti,1)
    incir = incircle_pred(x(ti(i,1),1), x(ti(i,1),2),...
                          x(ti(i,2),1), x(ti(i,2),2),...
                          x(ti(i,3),1), x(ti(i,3),2),...
                          x(:,1), x(:,2));
                  
    I = sum(incir);
    if (I>0)
        cnt = cnt + 1;
    end
                  
end