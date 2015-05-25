function [A, wij, Aeval] = compute_reconstruction_lhs(stencil,cell_num, moment, p, xc, ge_terms)

% p = [cell.reconstruction_param.px; cell.reconstruction_param.py];
% moment = cell.moment;

% for i = 1:size(cell.stencil,1)
pw = 1; %distance weighting parameter

if size(p,2)>1
  i = 1;
    xi = xc(cell_num,:);
%     xi = [0,0];
    I = find(stencil(i,:)==cell_num);
    for j = 1:size(stencil,2)
        xj = xc(stencil(i,j),:);


        rij = sqrt( (xj(1)-xi(1)).^2 + (xj(2)-xi(2)).^2 );
        wij(j,1)  = 1/rij^pw;
%         if (rij < eps) wij(j,1) = 1; end
        if (I==j) wij(j,1) = 1; end


        for nn = 1:size(p,2)
            m = p(2,nn);
            n = p(1,nn);

            s = 0;
            for l = 0:m
                for k = 0:n
                    cv = moment( stencil(i,j) ).cv(  n-k+1, m-l+1 );

                    s = s + factorial(m)/(factorial(l)*factorial(m-l))...
                    *factorial(n)/(factorial(k)*factorial(n-k))...
                    *(xj(1)-xi(1)).^k.*(xj(2)-xi(2)).^l.*cv;
                end
            end

        Aeval(j,nn) = s;
        A(j,nn)  = Aeval(j,nn)*wij(j,1);

        end
        
    end
    
%     wij = wij./max(wij);
%     wij = ones(size(wij));
%   B=A;
%   for j = 1:ge_terms
%       for i = j+1:size(A,1);
%           wij2(i,j) = B(i,j)./B(j,j);
%           B(i,j:end) = B(i,j:end) - wij2(i,j)*B(j,j:end);
%       end
%   end
%   A=B;
%   wij = wij2;





    I = find(stencil(i,:)==cell_num);
    A(I,:) = A(I,:)./A(I,1);
    row = A(I,:);
    cnt=1;
    for j = 1:size(A,1)
        if j~=I
            
            c = A(j,1)/row(1);
            
            Aout(cnt,:) = A(j,:)-c*row;

            wij_out(cnt,1) = wij(j,1);
            cnt = cnt + 1;
        else
            wij_out(cnt,1) = wij(j,1);
%             Aout(cnt,:) = zeros(size(row));
            Aout(cnt,:) = A(j,:);
            cnt = cnt + 1;
            
        end
    end
    A = Aout;
    wij = wij_out;
    
else
    A = 1;
    wij = 1;
    
end