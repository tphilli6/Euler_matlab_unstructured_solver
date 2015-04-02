function te_fun_test

coef = [1;1];
der = [1,0
       0,1];
% der = [0,1];
order = 2;
nx = 3;
nc = 2;
   
   [te_der] = te_1d(der, coef, order);

% te_der=[     1     1
%              2     1
%              0     2
%              0     3 ];
% %    
%    func.rho = @(x) 1+0.2*x.^4;
%    dfun.rho = @(x) 0.8*x.^3;
%    func.u = @(x) 10-3*x.^4;
%    dfun.u = @(x) -12*x.^3;
   
   exponent = [0,1,2,3,4];
%    Coef(1,:) = [1,-0.1,0.05,0.1];
   Coef(2,:) = [10, 2,-1,0.5,3];
   
   Coef(1,:) = [3,0,0,0,0];
%    Coef(2,:) = [100, 2,-1,0.5];
%    Coef(1,:) = [1,-0.1,0,0,0.2];
%    Coef(2,:) = [10, 1,0,0, 3];
         
   func.rho = Coef(1,:);
   func.u = Coef(2,:);
   func.exp = exponent;
%    func.rho = @(x) sum(Coef(1,:).*x.^exponent);
%    dfun.rho = @(x) sum(Coef(1,2:end).*exponent(2:end).*x.^(exponent(2:end)-1));
%    func.u = @(x) sum(Coef(2,:).*x.^exponent);
%    dfun.u = @(x) sum(Coef(2,2:end).*exponent(2:end).*x.^(exponent(2:end)-1));
   
   x00 = 0.5;
   dx0 = 0.01;
   x0vec = x00:dx0:x00+dx0*(size(te_der,1)-2);
%    x0vec = x00:dx0:x00+dx0*(40-1);
   for nn = 1:length(x0vec)
       
%        ns = [2,1];
%        scale = poly(x0vec(nn),Coef(ns(nn),:),exponent,0);
       for j = 1:size(te_der,2)
           for i = 1:size(te_der,1)-1
               der_mat(i,j) = poly(x0vec(nn),Coef(j,:),exponent,te_der(i,j));%./scale;
           end
       end
       row = prod(der_mat,2)';
       A(nn,:) = row;
   end
%    Ainv = (A'*A)^-1*A';

%    scale = max(A);
%    for i = 1:size(A,1)
%        A(i,:) = A(i,:)./scale;
%    end
   Ainv = A^-1;
   
   dx = 0.001;
   r=1.1;
   for i = 1:20

       for nn = 1:length(x0vec)
       
       x0 = x0vec(nn);
       
       h(i) = 1/r^(i-1);
       dxvec(i) = dx;
       x = dx*[0:nc+1]-x0;
       te(nn,i) = -(disc_scheme_2nd_order(x,func,nc) - cont_scheme(x,func,nc));

       end
       
       b = te(:,i)./dx.^2;
       te_coef(:,i) = Ainv*b;%.*scale';
%        R(i) = 1-sqrt( sum( (A*te_coef(:,i) - b).^2 )/(length(b)-length(te_coef(:,i))) );
   dx = dx/r;
   end
%    loglog(h,te);
   p = log(te(:,1:end-1)./te(:,2:end))./log(r);

   [N,D] = rat(te_coef(:,end),1e-3);
   pcoef = log( abs(te_coef(:,1:end-2) - te_coef(:,2:end-1))./abs(te_coef(:,2:end-1)-te_coef(:,3:end)))./log(r);
   te_coefs = te_coef(:,end) - (te_coef(:,end-1) - te_coef(:,end))/(r^1-1);
   [N,D] = rat(te_coefs,1e-3)
   semilogx(h./h(end),te_coef);

   
function f = poly(x,coef,p,der)
  if der == 0
        c = coef;
        f = sum(c.*x.^p);
  elseif der == 1
    c = coef(2:end).*p(2:end);
    f = sum( c.*x.^(p(2:end)-1) );
   
  elseif der == 2
    c = coef(3:end).*p(3:end).*(p(3:end)-1);
    f = sum( c.*x.^(p(3:end)-2) );  
  
  elseif der == 3
    c = coef(4:end).*p(4:end).*(p(4:end)-1).*(p(4:end)-2);
    f = sum( c.*x.^(p(4:end)-3) );  
    
  elseif der == 4
    c = coef(5:end).*p(5:end).*(p(5:end)-1).*(p(5:end)-2).*(p(5:end)-3);
    f = sum( c.*x.^(p(5:end)-4) );
    
  else
    f = 0;
    
  end
    
function [te] = te_1d(der, coef, order)
%% automatic truncation error calculation

% 1D Mass
% dx = [1, 0
%       0, 1];


 
% output
% coef = [1, 2, 1]';
% der  = [0, 1
%         1, 1
%         1, 0];
   
% dy = [1, 0
%       0, 1];
dmax = 1;
% order = 2;
% T = f(x+dx) = f(x) + df/dx*dx + d2f/dx2*dx^2 * dx^2/2 + d3f/dx3 * dx^3/6 + ...
T(1) = 0;
for n = 1:order+2
    T(n+1) = 1/factorial(n); %1/factorial(n)
end

for n = 1:order-dmax

    cnt = 1;
    for term = 1:size(der,1);

        for d = 1:size(der,2)
%             for dd = 1:size(der,2)

             if der(term,d)>0
              der_out(cnt,:) = der(term,:);
              der_out(cnt,d) = der(term,d)+2;
              cnt = cnt + 1;
%              end
             
             
            end
%              der_out(cnt,:) = der(term,:);
%              der_out(cnt,d) = der(term,d)+2;
%              cnt = cnt + 1;
        end
    end

    der = der_out;
    
end

term = 1;
der = sortrows(der);
for row = 1:size(der,1)-1
   cnt = 1;
   if (der(row,1)~=-1)
       te(term,:) = der(row,:);
       while der(row,:)== der(row+cnt,:)
         der(row+cnt,:) = -ones(1,size(der,2));
         cnt = cnt + 1;
       end
       te_coef(term,1) = cnt;
       term = term + 1;
   end

end

if (length(der(:,1))>1)
    if (der(end,1)~=-1)
       te = [te; der(end,:)];
       te_coef = [te_coef;1];
    end
else
    te = der;
%     te_coef = te_coef;
end

% [temp,I] = sort(te,2);
% [~,I2] = sortrows(temp);
% te = te(I2,:);

function [s] = cont_scheme(x,func, i)

rho = poly(x(i),func.rho,func.exp,0);
drdx = poly(x(i),func.rho,func.exp,1);

u = poly(x(i),func.u,func.exp,0);
dudx = poly(x(i),func.u,func.exp,1);

s = rho.*dudx + u.*drdx;
% s = drdx;

function [resid] = disc_scheme_2nd_order(x,func,i)

for n = 1:length(x)
    rho(n) = poly(x(n),func.rho,func.exp,0);
    u(n) = poly(x(n),func.u,func.exp,0);
end
resid = rho(2)*(u(3) - u(1))./(x(3)-x(1)) + u(2)*( rho(3) - rho(1) )./(x(3)-x(1));
% resid = ( rho(3) - rho(1) )./(x(3)-x(1));
% resid = ( -3/2*rho(i) + 2*rho(i+1) - 1/2*rho(i+2) )/(x(i+1)-x(i));
% resid = rho(i)*( u(i-2)/12 - 2*u(i-1)/3 + 2*u(i+1)/3 - u(i+2)/12) + ...
%         u(i)*( rho(i-2)/12 - 2*rho(i-1)/3 + 2*rho(i+1)/3 - rho(i+2)/12);