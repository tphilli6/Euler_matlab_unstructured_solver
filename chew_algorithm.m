function [xadd_override, bfn_out, bfun_lbl_out, fu] = chew_algorithm(xadd, bfn, x, bfun_lbl, bfun)
% Uses chews algorithm to test if added node is too close to boundary
% xadd        = [x, y]
% bfn         = boundary face nodes (n,2)
% x           = x and y node locations
% bfun_lbl    = label used to call a boundary function if needed
% bfun(k).fun = function y=f(x) to place new boundary point.

% bfn = boundary face nodes
dx = (x(bfn(:,2),1)-x(bfn(:,1),1));
dy = (x(bfn(:,2),2)-x(bfn(:,1),2));
thad1 = dx<0;
thad2 = dx>0 & dy<0;

th = atan(dy./dx) + thad1*pi + thad2*2*pi;

r = sqrt( dx.^2 + dy.^2 )/2;
xc = [x(bfn(:,2),1)+x(bfn(:,1),1),  x(bfn(:,2),2)+x(bfn(:,1),2)]/2;

xcir = [r.*cos(th+pi/2),r.*sin(th+pi/2)]+xc;

incirc = incircle_pred(x(bfn(:,2),1), x(bfn(:,2),2),...
                      xcir(:,1), xcir(:,2),...
                      x(bfn(:,1),1), x(bfn(:,1),2),...
                      xadd(1), xadd(2)  );

                  
if any(incirc)                  
    [J] = find(incirc==1);
    I = J(1);
    xadd_override = bfun(bfun_lbl(I)).fun( x(bfn(I,1),:),x(bfn(I,2),:) );

    fu = bfn(I,:);% Face to force removal at point at because point could be outside of the domain
    
    Inot = find(incirc==0);
    bfn_out = [bfn([1:I-1],:)
       bfn(I,1), length(x(:,1))+1
       length(x(:,1))+1, bfn(I,2)
       bfn([I+1:end],:)];

    bfun_lbl_out = [bfun_lbl(1:I-1)
            bfun_lbl(I)
            bfun_lbl(I)
            bfun_lbl(I+1:end)];
         
else
    xadd_override = xadd;
    bfn_out = bfn;
    bfun_lbl_out = bfun_lbl;
    fu = [0,0];
    
end

