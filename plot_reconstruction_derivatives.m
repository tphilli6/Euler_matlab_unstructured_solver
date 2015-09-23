function l2_err = plot_reconstruction_derivatives(coef, px, py, xc, vol, exact_fun, cell_vertex)

cell_vertex = [cell_vertex; cell_vertex(1,:)];

plotstuff=0;
l = size(coef,1);

% c = coef(:,1)'; 
i = find( px>=0 & py>=0 );
f = @(x,c) sum( c(i).*(x(1)-xc(1)).^px.*(x(2)-xc(2)).^py );

i = find( px-1>=0 & py>=0 );
fx = @(x,c) sum(c(i).*px(i).*(x(1)-xc(1)).^max(px(i)-1, zeros(1,length(i))).*(x(2)-xc(2)).^(py(i)) );
i = find( px>=0 & py-1>=0 );
fy = @(x,c) sum(c(i).*py(i).*(x(1)-xc(1)).^(px(i)).*(x(2)-xc(2)).^max(py(i)-1, zeros(1,length(i))) );

i = find( px-2>=0 & py>=0 );
fxx = @(x,c) sum(c(i).*px(i).*(px(i)-1).*(x(1)-xc(1)).^(px(i)-2).*(x(2)-xc(2)).^(py(i)));
i = find( px-1>=0 & py-1>=0 );
fxy = @(x,c) sum(c(i).*px(i).*py(i).*(x(1)-xc(1)).^(px(i)-1).*(x(2)-xc(2)).^(py(i)-1));
i = find( px>=0 & py-2>=0 );
fyy = @(x,c) sum(c(i).*py(i).*(py(i)-1).*(x(1)-xc(1)).^(px(i)).*(x(2)-xc(2)).^(py(i)-2));

i = find( px-3>=0 & py>=0 );
fxxx = @(x,c) sum(c(i).*px(i).*(px(i)-1).*(px(i)-2).*(x(1)-xc(1)).^(px(i)-3).*(x(2)-xc(2)).^(py(i)));
i = find( px-2>=0 & py-1>=0 );
fxxy = @(x,c) sum(c(i).*px(i).*(px(i)-1).*(py(i)).*(x(1)-xc(1)).^(px(i)-2).*(x(2)-xc(2)).^(py(i)-1));
i = find( px-1>=0 & py-2>=0 );
fxyy = @(x,c) sum(c(i).*px(i).*(py(i)).*(py(i)-1).*(x(1)-xc(1)).^(px(i)-1).*(x(2)-xc(2)).^(py(i)-2));
i = find( px>=0 & py-3>=0 );
fyyy = @(x,c) sum(c(i).*py(i).*(py(i)-1).*(py(i)-2).*(x(1)-xc(1)).^(px(i)).*(x(2)-xc(2)).^(py(i)-3));

i = find( px-4>=0 & py>=0 );
fxxxx = @(x,c) sum(c(i).*px(i).*(px(i)-1).*(px(i)-2).*(px(i)-3).*(x(1)-xc(1)).^(px(i)-4).*(x(2)-xc(2)).^(py(i)));
i = find( px-3>=0 & py-1>=0 );
fxxxy = @(x,c) sum(c(i).*px(i).*(px(i)-1).*(px(i)-2).*py(i).*(x(1)-xc(1)).^(px(i)-3).*(x(2)-xc(2)).^(py(i)-1));
i = find( px-2>=0 & py-2>=0 );
fxxyy = @(x,c) sum(c(i).*px(i).*(px(i)-1).*(py(i)).*(py(i)-1).*(x(1)-xc(1)).^(px(i)-2).*(x(2)-xc(2)).^(py(i)-2));
i = find( px-1>=0 & py-3>=0 );
fxyyy = @(x,c) sum(c(i).*px(i).*(py(i)).*(py(i)-1).*(py(i)-2).*(x(1)-xc(1)).^(px(i)-1).*(x(2)-xc(2)).^(py(i)-3));
i = find( px>=0 & py-4>=0 );
fyyyy = @(x,c) sum(c(i).*py(i).*(py(i)-1).*(py(i)-2).*(py(i)-3).*(x(1)-xc(1)).^(px(i)).*(x(2)-xc(2)).^(py(i)-4));

dx = sqrt(vol);

%% Create evaluation array

% Create square array
n = 10;
xvec = linspace(-1*dx, 1*dx, n);
[y,x] = meshgrid(xvec,xvec);
x = reshape(x,[numel(x),1])+xc(1);
y = reshape(y,[numel(y),1])+xc(2);

xsquare = [x,y];

% Create face integral array
[ xq, wq ] = gauss_patterson( max( max(px,py)) ); adj = 0;

lxq = length(xq);
nvc = size(cell_vertex,1);
xline = zeros( (nvc-1)*(lxq-adj),2);
for i = 1:nvc-1;
    xn = cell_vertex(i:i+1,:);
      xline( (i-1)*(lxq-adj)+1:i*(lxq-adj), 1) = (xn(2,1) - xn(1,1)) .* xq(1:end-adj) + xn(1,1);
      xline( (i-1)*(lxq-adj)+1:i*(lxq-adj), 2) = (xn(2,2) - xn(1,2)) .* xq(1:end-adj) + xn(1,2);
end

% Create area integral array
[ xq, wq ] = sparse_grid( 2, max( max(px,py)), 'GP' );
xfit = [0, 0
        1, 0
        0, 1
        1, 1];
A = [ones(4,1), xfit(:,1), xfit(:,2), xfit(:,1).*xfit(:,2)];
Ainv = A^-1;
ffit = @(xi,ax,ay) [ax(1)+ax(2)*xi(:,1)+ax(3)*xi(:,2)+ax(4)*xi(:,1).*xi(:,2), ay(1)+ay(2)*xi(:,1)+ay(3)*xi(:,2)+ay(4)*xi(:,1).*xi(:,2)];

lxq = length(xq);
nvc = size(cell_vertex,1);
xarea = [];
for i = 1:nvc-1;
    xn = [cell_vertex(i:i+1,:); xc; xc];
    ax = Ainv*xn(:,1);
    ay = Ainv*xn(:,2);
    xtri = ffit(xq, ax, ay);
    xarea = [xarea; xtri];
end

% Plot eval points
% plot(cell_vertex(:,1), cell_vertex(:,2),'k-','LineWidth',5)
% axis equal
% hold on
% 
% nvc = size(cell_vertex,1);
% for i = 1:nvc-1;
%     xn = [cell_vertex(i:i+1,:); xc; cell_vertex(i,:)];
%     plot(xn(:,1), xn(:,2), 'k-')
% end
% plot(xsquare(:,1), xsquare(:,2), 'ro',...
%      xline(:,1), xline(:,2), 'g*',...
%      xarea(:,1), xarea(:,2), 'b*',...
%      xc(1), xc(2), 'ko');
% title(['Derivative Evaluation Points for k=',num2str(max(max(px,py)))])
% plot(xsquare(:,1), xsquare(:,2), 'ro');
% plot(xline(:,1), xline(:,2), 'g*');
% plot(xarea(:,1), xarea(:,2), 'b*');
% plot(xc(1), xc(2), 'ko')

x = xc;

for j = 1:size(coef,2);
      if j==1
          [fex_comp, f_ex_comp] = fun_rho(xc, exact_fun);
      elseif j==2
         [fex_comp, f_ex_comp] = fun_u(xc, exact_fun);
      elseif j==3
         [fex_comp, f_ex_comp] = fun_v(xc, exact_fun);         
      elseif j==4
         [fex_comp, f_ex_comp] = fun_p(xc, exact_fun);
      end
      fex_norm = [fex_comp, f_ex_comp]; 
    
    ctest = coef(:,j)';
    for i = 1:size(x,1);
      if j==1
          [fex(i,:), f_ex(i,:)] = fun_rho(x(i,:), exact_fun);
      elseif j==2
         [fex(i,:), f_ex(i,:)] = fun_u(x(i,:), exact_fun); 
      elseif j==3
         [fex(i,:), f_ex(i,:)] = fun_v(x(i,:), exact_fun);         
      elseif j==4
         [fex(i,:), f_ex(i,:)] = fun_p(x(i,:), exact_fun);
      end


        ftest(i,1) = f( x(i,:), ctest );

        f_test(i,1) = fx( x(i,:), ctest );
        f_test(i,2) = fy( x(i,:), ctest );

        f_test(i,3) = fxx( x(i,:), ctest );
        f_test(i,4) = fxy( x(i,:), ctest );
        f_test(i,5) = fyy( x(i,:), ctest );

        f_test(i,6) = fxxx( x(i,:), ctest );
        f_test(i,7) = fxxy( x(i,:), ctest );
        f_test(i,8) = fxyy( x(i,:), ctest );
        f_test(i,9) = fyyy( x(i,:), ctest );  

        f_test(i,10) = fxxxx( x(i,:), ctest );
        f_test(i,11) = fxxxy( x(i,:), ctest );
        f_test(i,12) = fxxyy( x(i,:), ctest );
        f_test(i,13) = fxyyy( x(i,:), ctest );
        f_test(i,14) = fyyyy( x(i,:), ctest );


    end

      l2_fun = sqrt( sum([fex.^2,f_ex.^2])/length(fex) );
    
      err_fun = ftest - fex;
      l2_err_fun(1) = sqrt( sum(err_fun.^2)/length(err_fun) );
      err_der = f_test - f_ex;
      l2_err_der(1,:) = sqrt( sum(err_der.^2,1)/size(err_der,1) );




  if plotstuff
    xmat = reshape(x(:,1),[n,n]);
    ymat = reshape(x(:,2),[n,n]);

    fexmat = reshape(fex, [n,n]);
    f_der_ex = reshape(f_ex, [n,n,size(f_ex,2)]);

    ftestmat = reshape(ftest, [n,n]);
    f_der_test = reshape(f_test, [n,n,size(f_ex,2)]);

    f_der_err_mat = f_der_test - f_der_ex;
    
    
    % for i = 1:size(f_der_ex,3)
        subplot(2,6,1)
        surf(xmat,ymat,ftestmat - fexmat);
        hold on
%         surf(xmat,ymat,ftestmat(:,:));
       plot3(cell_vertex(:,1),cell_vertex(:,2),ones(size(cell_vertex(:,1)))*max(max(ftestmat - fexmat)),'ko-');
       title(['L2 Error: ',num2str(l2_err_fun)]);
    %     title('f(x)');

        for jj = 1:11
            subplot(2,6,jj+1)
            surf(xmat,ymat,f_der_err_mat(:,:,jj));
            hold on
    %         surf(xmat,ymat,f_der_test(:,:,1));
            plot3(cell_vertex(:,1),cell_vertex(:,2),ones(size(cell_vertex(:,1)))*max(max(f_der_err_mat(:,:,jj))),'ko-');
            title(['L2 Error: ',num2str(l2_err_der(1,jj))])
        end
%         subplot(2,6,3)
%         surf(xmat,ymat,f_der_ex(:,:,2));
% 
%         subplot(2,6,4)
%         surf(xmat,ymat,f_der_ex(:,:,3));
% 
%         subplot(2,6,5)
%         surf(xmat,ymat,f_der_ex(:,:,4));
% 
%         subplot(2,6,6)
%         surf(xmat,ymat,f_der_ex(:,:,5));
% 
%         subplot(2,6,7)
%         surf(xmat,ymat,ftestmat(:,:));
%     %     title('f(x)');
% 
%         subplot(2,6,8)
%         surf(xmat,ymat,f_der_test(:,:,1));
% 
%         subplot(2,6,9)
%         surf(xmat,ymat,f_der_test(:,:,2));
% 
%         subplot(2,6,10)
%         surf(xmat,ymat,f_der_test(:,:,3));
% 
%         subplot(2,6,11)
%         surf(xmat,ymat,f_der_test(:,:,4));
% 
%         subplot(2,6,12)
%         surf(xmat,ymat,f_der_test(:,:,5));
    % end
    
  end
  
  l2_err(1,j,:) = [l2_err_fun', l2_err_der];

  
end



function [fex, f_ex] = fun_rho(x, exact_fun);

    fex(1,:) = exact_fun.rho(x);

    f_ex(1,1) = exact_fun.drdx(x);
    f_ex(1,2) = exact_fun.drdy(x);

    f_ex(1,3) = exact_fun.drdxx(x);
    f_ex(1,4) = exact_fun.drdxy(x);
    f_ex(1,5) = exact_fun.drdyy(x);

    f_ex(1,6) = exact_fun.drdxxx(x);
    f_ex(1,7) = exact_fun.drdxxy(x);
    f_ex(1,8) = exact_fun.drdxyy(x);
    f_ex(1,9) = exact_fun.drdyyy(x);

    f_ex(1,10) = exact_fun.drdxxxx(x);
    f_ex(1,11) = exact_fun.drdxxxy(x);
    f_ex(1,12) = exact_fun.drdxxyy(x);
    f_ex(1,13) = exact_fun.drdxyyy(x);
    f_ex(1,14) = exact_fun.drdyyyy(x);
    
    
function [fex, f_ex] = fun_u(x, exact_fun)

    fex(1,:) = exact_fun.u(x);

    f_ex(1,1) = exact_fun.dudx(x);
    f_ex(1,2) = exact_fun.dudy(x);

    f_ex(1,3) = exact_fun.dudxx(x);
    f_ex(1,4) = exact_fun.dudxy(x);
    f_ex(1,5) = exact_fun.dudyy(x);

    f_ex(1,6) = exact_fun.dudxxx(x);
    f_ex(1,7) = exact_fun.dudxxy(x);
    f_ex(1,8) = exact_fun.dudxyy(x);
    f_ex(1,9) = exact_fun.dudyyy(x);

    f_ex(1,10) = exact_fun.dudxxxx(x);
    f_ex(1,11) = exact_fun.dudxxxy(x);
    f_ex(1,12) = exact_fun.dudxxyy(x);
    f_ex(1,13) = exact_fun.dudxyyy(x);
    f_ex(1,14) = exact_fun.dudyyyy(x);
    
function [fex, f_ex] = fun_v(x, exact_fun)

    fex(1,:) = exact_fun.v(x);

    f_ex(1,1) = exact_fun.dvdx(x);
    f_ex(1,2) = exact_fun.dvdy(x);

    f_ex(1,3) = exact_fun.dvdxx(x);
    f_ex(1,4) = exact_fun.dvdxy(x);
    f_ex(1,5) = exact_fun.dvdyy(x);

    f_ex(1,6) = exact_fun.dvdxxx(x);
    f_ex(1,7) = exact_fun.dvdxxy(x);
    f_ex(1,8) = exact_fun.dvdxyy(x);
    f_ex(1,9) = exact_fun.dvdyyy(x);

    f_ex(1,10) = exact_fun.dvdxxxx(x);
    f_ex(1,11) = exact_fun.dvdxxxy(x);
    f_ex(1,12) = exact_fun.dvdxxyy(x);
    f_ex(1,13) = exact_fun.dvdxyyy(x);
    f_ex(1,14) = exact_fun.dvdyyyy(x);
    
function [fex, f_ex] = fun_p(x, exact_fun)

    fex(1,:) = exact_fun.p(x);

    f_ex(1,1) = exact_fun.dpdx(x);
    f_ex(1,2) = exact_fun.dpdy(x);

    f_ex(1,3) = exact_fun.dpdxx(x);
    f_ex(1,4) = exact_fun.dpdxy(x);
    f_ex(1,5) = exact_fun.dpdyy(x);

    f_ex(1,6) = exact_fun.dpdxxx(x);
    f_ex(1,7) = exact_fun.dpdxxy(x);
    f_ex(1,8) = exact_fun.dpdxyy(x);
    f_ex(1,9) = exact_fun.dpdyyy(x);

    f_ex(1,10) = exact_fun.dpdxxxx(x);
    f_ex(1,11) = exact_fun.dpdxxxy(x);
    f_ex(1,12) = exact_fun.dpdxxyy(x);
    f_ex(1,13) = exact_fun.dpdxyyy(x);
    f_ex(1,14) = exact_fun.dpdyyyy(x);
