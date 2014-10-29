function area = triangle_area(x1, x2, x3)
% Compute the angle of a triangle using the cross product of two vectors
% V1 = x1 - x2
% V2 = x3 - x2
% Area = 1/2* (V2 x V1) (could also be 1/2*(V1 x V2) but would be negative area.)


V1 = x1 - x2;
V2 = x3 - x2;
area_vec = 0.5* cross( V2, V1 );
area = area_vec(3);
