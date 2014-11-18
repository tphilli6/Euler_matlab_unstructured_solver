function [s] = KahanSum(data)

% s = sum(data);

[~,I] = sort(abs(data));
data = data(I);

s=0;
c=0;


for i = 1:length(data);
    y = data(i)-c;
    t = s + y;
    c = (t-s) - y;
    s = t;
end