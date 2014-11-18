function area = quad_area(x,y)

for j=1:size(x,2)-1
    for i =1:size(x,1)-1
        x0 = [x(i,j),y(i,j)];
        x1 = [x(i+1,j),y(i+1,j)];
        x2 = [x(i,j+1),y(i,j+1)];
        x3 = [x(i+1,j+1),y(i+1,j+1)];
        
        a = x1-x0;
        b = x2-x0;
        a1 = a(1)*b(2)-a(2)*b(1);
        
        a = x2-x3;
        b = x1-x3;
        a2 = a(1)*b(2)-a(2)*b(1);
        
        
        
        area(i,j) = (a1+a2)/2;
    end
end