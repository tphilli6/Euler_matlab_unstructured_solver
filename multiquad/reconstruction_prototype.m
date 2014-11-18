clc
clear all

%% Define problem
ncells = 3;
xi(:,:) = meshgrid([0:1/ncells:1],[0:1/ncells:1])';
eta(:,:) = xi(:,:)';

% create physical grid
dsize=1;
x = 0 + dsize*xi + 0*eta + 0*xi.*eta;
y = 0 + 0*xi + dsize*eta + 0*xi.*eta;

F = @(x,y) 2+sin(x/7*pi) + cos(y/10*pi);
% F = @(x,y) 1 + x + x.^2 + y + x.*y;%  + x.^2.*y + y.^2 + x.*y.^2 + x.^2.*y.^2;
% F = @(x,y) 1 + x + y + x.*y;


%% compute local cell mapping
% polynomial curvfit
% x = cx1 + cx2*xi + cx3*eta + cx4*xi*eta
% y = cy1 + cy2*xi + cy3*eta + cy4*xi*eta
% Domain x = [0,1], y=[0,1]
A = [1, 0, 0, 0
     1, 1, 0, 0
     1, 0, 1, 0
     1, 1, 1, 1];
Ainv = A^-1;
for j = 1:ncells
    for i = 1:ncells
        cx(:,i,j) = Ainv*[x(i,j),x(i+1,j),x(i,j+1),x(i+1,j+1)]';
        cy(:,i,j) = Ainv*[y(i,j),y(i+1,j),y(i,j+1),y(i+1,j+1)]';
    end
end
gmap = @(x,y,c)c(1) + c(2)*x + c(3)*y + c(4).*x.*y;


%% Compute solution
[xi, w] = sparse_grid( 2, 4 );

for j = 1:ncells
    for i = 1:ncells
        
        soln(i,j) = F(gmap(xi(:,1),xi(:,2),cx(:,i,j)),gmap(xi(:,1),xi(:,2),cy(:,i,j)))'*w;%multiply by Jacobian but divided out when computing the average
        
    end
end
b = reshape(soln,[numel(soln),1]);

%% compute global cell mapping (to fit domain)
% specific to order and stencil
dx = 1/ncells;
dy = dx;
for j = 1:ncells
    for i = 1:ncells
        x0 = (i-1)*dx;
        y0 = (j-1)*dy;
        
        clocx(:,i,j) = [x0, dx];
        clocy(:,i,j) = [y0, dy];
        
    end
end
l_to_g = @(x,c) c(1) + x.*c(2);
g_to_l = @(x,c) (x-c(1))./c(2);
g_to_l_coef = @(c) [-c(1)/c(2), 1/c(2)];
l_to_g_coef = @(c) [c(1), c(2)];

%% reconstruction

[xi_cg, w_cg] = sparse_grid( 2, 2*ncells );

% build polynomial array
cnt = 1;
for j=1:ncells
    for i = 1:ncells
        p(1,cnt) = i-1;
        p(2,cnt) = j-1;
        cnt = cnt + 1;
    end
end
Peval = @(x,y,c) sum(c.*x.^p(1,:).*y.^p(2,:));
Arow = @(x,y) x.^p(1,:).*y.^p(2,:);

cnt = 1;
for j=1:ncells
    for i = 1:ncells
%         l_to_g(x(i,j),clocx(:,i,j))
%         g_to_l(x(i,j),clocx(:,i,j))
        %LHS for one cell
        for n = 1:length(w_cg)
            Alhs(n,:) = Arow(l_to_g(xi_cg(n,1),clocx(:,i,j)),l_to_g(xi_cg(n,2),clocy(:,i,j)));            
        end
        %Afit is computed for each cell and combined for which ever
        %reconstruction method is to be used
        Afit(cnt,:) = w_cg'*Alhs;
        cnt = cnt + 1;
    end
end


cfit = Afit\b;
% cfit = ones(size(cfit));
%% Shift coordinates for the global cell coordinates to local cell coordinates
[P,expnt] = create_poisson_triangle(ncells-1);
cnew = zeros(numel(cfit),1);
% for j = 1:ncells
%     for i = 1:ncells
i=2;
j=2;
        ax=l_to_g_coef(clocx(:,i,j));
        ay=l_to_g_coef(clocy(:,i,j));
        cnew = zeros(size(cfit));
        for n = 1:length(cfit)
            px = p(1,n);
            ii = 1:px+1;
            coefx = P(px+1,ii).*ax(1).^expnt(px+1,ii,1).*ax(2).^expnt(px+1,ii,2);
            
            py = p(2,n);
            ii = 1:py+1;
            coefy = P(py+1,ii).*ay(1).^expnt(py+1,ii,1).*ay(2).^expnt(py+1,ii,2);
            
            
%             coefins = coefx'*coefy;
%             coef = zeros(max(p(1,:)+1),max(p(2,:)+1));
%             coef(1:px+1,1:py+1) = coefins;
%             cnew(1:numel(coef)) = reshape(coef,[numel(coef),1]) + cnew(1:numel(coef));
            
            for jj = 1:py+1
                for ii = 1:px+1
                    cindx = ii + max(p(1,:)+1)*(jj-1);
                    cnew(cindx) = cfit(n)*coefx(ii)*coefy(jj) + cnew(cindx);                    
%                     cnew(cindx,i,j) = coefx(ii)*coefy(jj) + cnew(cindx,i,j);                    
                end
            end
        end
%         reshape(cnew,[3,3])



fprintf('Original polynomial evaluated in subdomain: %12.4e\n', Peval(l_to_g(1/2,clocx(:,i,j)),l_to_g(1/2,clocy(:,i,j)), cfit') );
fprintf('Original polynomial evaluated in cell domain: %12.4e\n', Peval(1/2, 1/2, cnew' ) );
fprintf('Transfer Error: %23.15e\n\n\n',  Peval(1/2, 1/2, cnew') -Peval(l_to_g(1/2,clocx(:,i,j)),l_to_g(1/2,clocy(:,i,j)), cfit') );
%     end
% end

%% Test fit
pfit = @(x,y) sum(cfit'.*x.^p(1,:).*y.^p(2,:));

for j = 1:ncells
    for i = 1:ncells
        for n = 1:length(w_cg);
            soln_temp(n) = pfit(l_to_g(xi_cg(n,1),clocx(:,i,j)),l_to_g(xi_cg(n,2),clocy(:,i,j))).*w_cg(n);%multiply by Jacobian but divided out when computing the average
        end
        soln_test(i,j) = sum(soln_temp);
%         soln_test(i,j) = soln_test(i,j) * dx*dy;
    end
end

error = sqrt( sum(reshape( (soln_test - soln).^2, [numel(soln),1]))/numel(soln) );
fprintf('Solution reconstruction error in average value: %12.4e\n',error);


