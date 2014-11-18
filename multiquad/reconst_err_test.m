function maxerr = reconst_err_test(order, scale)
% clc
% clear all
dimen=2;
% order=3;
% scale=10;
maxorder=order;

% skewed grid
    mapx = [0, 1, 0, 0.0]*scale;
    mapy = [0, 0, 1, -0.0]*scale;
    
    x_xi = @(xi,eta)mapx(1)+mapx(2).*xi+mapx(3).*eta+mapx(4)*xi.*eta;
    y_xi = @(xi,eta)mapy(1)+mapy(2)*xi+mapy(3)*eta+mapy(4)*xi.*eta;
    
    X = linspace(0,1,order+2);
    % X = [0,1,3,6];
    Y = X;
    [xi,eta]=meshgrid(X,Y);
    xi=xi';
    eta=eta';
    
    x = x_xi(xi,eta);
    y = y_xi(xi,eta);



% regular grid
%     mapx = [0, 2, 0, 0];
%     mapy = [0, 0, 2, 0];
% 
%     x_xi = @(xi,eta)mapx(1)+mapx(2).*xi+mapx(3).*eta+mapx(4)*xi.*eta;
%     y_xi = @(xi,eta)mapy(1)+mapy(2)*xi+mapy(3)*eta+mapy(4)*xi.*eta;
% 
% %     X = linspace(0,3,4);
%     X = [0,1,3,6];
%     Y = X;
%     [xi,eta]=meshgrid(X,Y);
%     
%     x = x_xi(xi,eta);
%     y = y_xi(xi,eta);
% 
%     x=x';
%     y=y';


mesh(x,y,zeros(size(x)))
view([0,0,90])
axis equal

Area = quad_area(x,y);
Area_xi = quad_area(xi',eta');

coefs_ones = ones(1,(maxorder+1).^2);
coefs_ex = zeros(size(coefs_ones));
coefs_ex((order+1)^2)=1;

pmat(:,:,1) = repmat([0:maxorder]',1,maxorder+1);
pmat(:,:,2) = pmat(:,:,1)';
p(1,:) = reshape(pmat(:,:,1),1,(maxorder+1)^2);
p(2,:) = reshape(pmat(:,:,2),1,(maxorder+1)^2);

genPfunc =@(x,y,c) c.*x.^p(1,:).*y.^p(2,:);
Pfunc = @(x,y)genPfunc(x,y,coefs_ex,p);

coefs_exint = coefs_ex./(p(1,:)+1)./(p(2,:)+1);
pint(1,:) = p(1,:)+1;
pint(2,:) = p(2,:)+1;

% Pfunc =@(x,y) (1 + x + x.^2+y+y.*x+y.*x.^2+y.^2+y.^2.*x+y.^2.*x.^2)/10;
ax = mapx(1); bx = mapx(2); cx = mapx(3); dx = mapx(4);
ay = mapy(1); by = mapy(2); cy = mapy(3); dy = mapy(4);
J = @(xi,eta)abs( bx.*cy - by.*cx - cx.*dy.*eta + cy.*dx.*eta + bx.*dy.*xi - by.*dx.*xi );
% Jac = @(xi,eta)abs((mapx(2)+mapx(4)*eta).*(mapy(3)+mapy(4)*xi)-(mapx(3)+mapx(4)*xi).*(mapy(2)+mapy(4)*eta));
Pfunc_xi = @(xi,eta) sum(coefs_exint.*x_xi(xi,eta).^pint(1,:).*y_xi(xi,eta).^pint(2,:));


soln = zeros(size(xi)-1);
for j = 1:size(x,2)-1
     for i=1:size(x,1)-1
         Pex_eta =@(eta) Pfunc_xi(xi(i+1,j),eta).*J(xi(i+1,j),eta)-Pfunc_xi(xi(i,j),eta).*J(xi(i,j),eta);
         soln(i,j) = (Pex_eta(eta(i,j+1))-Pex_eta(eta(i,j)))./Area(i,j);
        
    end
end
% soln_error = soln-soln2;
quad_order = sum((order)*2) + 1 - 2;
[xi_sg,w] = sparse_grid(dimen,quad_order);
xirange=[0,0
    1,0
    0,1
    1,1];

Pmat=@(x,y)genPfuncmat(x,y,coefs_ones,p);

%scale domain
minx = min(min(x));
miny = min(min(y));
scalex = (max(max(x))-min(min(x)))/2;
scaley = (max(max(y))-min(min(y)))/2;

solnmax = max(max(abs(soln)));
soln = soln/solnmax;

x = (x-minx)/scalex-1;
y = (y-miny)/scaley-1;

cnt=1;
for j=1:order+1
    for i=1:order+1
        xvec=[x(i,j),y(i,j)
              x(i+1,j),y(i+1,j)
              x(i,j+1),y(i,j+1)
              x(i+1,j+1),y(i+1,j+1)];
          
        A = [ones(4,1),xirange(:,1),xirange(:,2),xirange(:,1).*xirange(:,2)];
        mapxf = A\xvec(:,1);
        mapyf = A\xvec(:,2);
        
%         ax = mapxf(1); bx = mapxf(2); cx = mapxf(3); dx = mapxf(4);
%         ay = mapyf(1); by = mapyf(2); cy = mapyf(3); dy = mapyf(4);
%         Jf = @(xi,eta)abs( bx.*cy - by.*cx - cx.*dy.*eta + cy.*dx.*eta + bx.*dy.*xi - by.*dx.*xi );
        
        x_xif = @(xi)mapxf(1)+mapxf(2).*xi(1)+mapxf(3).*xi(2)+mapxf(4)*xi(1).*xi(2);
        y_xif = @(xi)mapyf(1)+mapyf(2)*xi(1)+mapyf(3)*xi(2)+mapyf(4)*xi(1).*xi(2);
%         Pmat=@(xi)[1,x_xif(xi),x_xif(xi).^2,y_xif(xi),x_xif(xi).*y_xif(xi),x_xif(xi).^2*y_xif(xi),y_xif(xi).^2,x_xif(xi).*y_xif(xi).^2,x_xif(xi).^2*y_xif(xi).^2];
        for ii=1:length(w)
          P(ii,:) = Pmat(x_xif(xi_sg(ii,:)),y_xif(xi_sg(ii,:)));
        end
        
        lhs_cell(cnt,:)=w'*P;
%         lhs_cell(cnt,:) = KahanMatMul1(w',P);
%         disp(max(abs(test-lhs_cell(cnt,:))))
        b(cnt,1) = soln(i,j);
        

        cnt=cnt+1;
    end
end
% disp(lhs_cell)
coefs_compspace=lhs_cell\b;



%Transform exact coefs to compare
transx=[(minx+scalex),scalex];
transy=[(miny+scaley),scaley];
coefs_ex_tran = transform_poly(coefs_ex,transx,transy,order)/solnmax*scale^2;

% transx=[minx,scalex];
% transy=[miny,scaley];
transx=[-(minx+scalex)/scalex,1/scalex];
transy=[-(miny+scaley)/scaley,1/scaley];
coefs = transform_poly(coefs_compspace,transx,transy,order)*solnmax/(scale^2);


% Compute solution using curve-fit results
% coefs=0.1*ones(9,1);
coefs_zero=zeros(maxorder+1,maxorder+1);
coefs=reshape(coefs,order+1,order+1);
coefs_zero(1:order+1,1:order+1)=coefs;
coefs=reshape(coefs_zero,1,(maxorder+1)^2);
Pfunc =@(x,y) genPfunc(x,y,coefs,p);
ax = mapx(1); bx = mapx(2); cx = mapx(3); dx = mapx(4);
ay = mapy(1); by = mapy(2); cy = mapy(3); dy = mapy(4);
J = @(xi,eta)abs( bx.*cy - by.*cx - cx.*dy.*eta + cy.*dx.*eta + bx.*dy.*xi - by.*dx.*xi );
% Jac = @(xi,eta)abs((mapx(2)+mapx(4)*eta).*(mapy(3)+mapy(4)*xi)-(mapx(3)+mapx(4)*xi).*(mapy(2)+mapy(4)*eta));
Pfunc_xi = @(xi,eta)Pfunc(x_xi(xi,eta),y_xi(xi,eta)).*J(xi,eta);


% Compute exact solution
% n=401;
% matlabpool(3)
% solntest = zeros(3,3);
% for j = 1:size(x,2)-1
% %     parfor i=1:size(x,1)-1
%      for i=1:size(x,1)-1
% 
% 
%         solntest(i,j) = quad2d( Pfunc_xi,xi(i,j),  xi(i,j+1),...
%                                    eta(i,j), eta(i+1,j),...
%                                    'AbsTol',1e-15,'RelTol',3e-14 )/Area(i,j);
% %         soln2(i,j) = quad2d( Pfunc,x_xi(xi(i,j),eta(i,j)),  x_xi(xi(i,j+1),eta(i,j+1)),...
% %                                    y_xi(xi(i,j),eta(i,j)),y_xi(xi(i+1,j),eta(i+1,j)),...
% %                                    'AbsTol',1e-10,'RelTol',1e-10)/Area(i,j);
%     end
% end

I=find(abs(coefs)>1e-12);
error = coefs(I)-coefs_ex(I);
error_trans = coefs_compspace - coefs_ex_tran;
% maxerr = max(abs(error)./coefs_ex(I));
maxerr = max(abs(error));
% maxerr = max(abs(error_trans));
% (solntest-soln)./soln
% maxerr=max(max(abs(solntest-soln)./soln));

end

function F=genPfunc(X,Y,c,p)

x = reshape(X,numel(X),1);
y = reshape(Y,numel(Y),1);

for i=1:numel(x)
f(i) = sum(c.*x(i).^p(1,:).*y(i).^p(2,:));
end

F = reshape(f,size(X));

end

function f=genPfuncmat(X,Y,c,p)

x = reshape(X,numel(X),1);
y = reshape(Y,numel(Y),1);

I=find(c~=0);

for i=1:numel(x)
f(i,:) = c(I).*x(i).^p(1,I).*y(i).^p(2,I);
end

% F = reshape(f,size(X));

end

function c=bincoef(n,k)
c=n.^k./factorial(k);
end

function coef_out = transform_poly(coef,xtrans,ytrans,ordermax)

x0=xtrans(1);
Sx=xtrans(2);
y0=ytrans(1);
Sy=ytrans(2);

%generate pascals triangle
Ptriangle(1,1)=1;
Ptriangle(2,1:2)=1;
for i=3:ordermax+1
 Ptriangle(i,1)=1;
 Ptriangle(i,2:i-1)=Ptriangle(i-1,1:end-1)+Ptriangle(i-1,2:end);
 Ptriangle(i,i)=1;    
end

coefmat = reshape(coef,ordermax+1,ordermax+1);

coef_out = zeros(ordermax+1,ordermax+1);
cnt=1;
for j=0:ordermax
    for i=0:ordermax
        %x expansion
        exp0=[i:-1:0];
        expS=[0:i];
        coefx = reshape(Ptriangle(i+1,1:i+1).*x0.^exp0.*Sx.^expS,length(exp0),1);
        
        %x expansion
        exp0=[j:-1:0];
        expS=[0:j];
        coefy = reshape(Ptriangle(j+1,1:j+1).*y0.^exp0.*Sy.^expS,1,length(exp0));
        
        cterm=coefmat(i+1,j+1)*coefx*coefy;
%         coef_out(1:i+1,1:j+1)=coef_out(1:i+1,1:j+1)+cterm;
        new_coef(1:i+1,1:j+1,cnt)=cterm;
        cnt = cnt + 1;
    end
end

for j=1:size(coef_out,2)
    for i = 1:size(coef_out,1)
        [coef_summed(i,j)] = KahanSum(new_coef(i,j,:));
    end
end

coef_out=reshape(coef_summed,(ordermax+1)^2,1);
% coef_out=reshape(coef_out,(ordermax+1)^2,1);


end