function [A] = KahanMatMul1(M1,M2)

for i = 1:size(M2,2)
    A(i) = KahanSum(M1.*M2(:,i)');
end