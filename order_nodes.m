function tip = order_nodes(n, v)

[area] = compute_area(v(n(:,1),1), v(n(:,1),2), ...
                      v(n(:,2),1), v(n(:,2),2), ...
                      v(n(:,3),1), v(n(:,3),2) );
                  
I = find(area<0);
tip=n;
tip(I,[1,3]) = n(I,[3,1]);

% n=tip;
% [area] = compute_area(v(n(:,1),1), v(n(:,1),2), ...
%                       v(n(:,2),1), v(n(:,2),2), ...
%                       v(n(:,3),1), v(n(:,3),2) );
