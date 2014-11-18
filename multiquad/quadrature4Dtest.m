% This code tests multi-dimensional quadratures to determine the combination
% of monomials that can be exactly integrated.
% up to 4D
% Polynomial form: P = x^px*y^py*z^pz*t^pt
% Author: Tyrone Phillips

clear all
clc

%% Higher dimension quadrature test
dimen=2; % up to 4
maxorder = 10;
method = 'CC';

x1=1;
x0=0;

rule_fail = 0;

if dimen>2; usez=1; else usez=0; end
if dimen>3; uset=1; else uset=0; end

[xx,yy] = meshgrid([0.5:maxorder+2.5],[0.5:maxorder+2.5]);

for o = 2:maxorder

    %compute appropriate polynomial order for given quadrature order
    porder = o;
%     if (mod(o,2)==0); porder = o-1; end

    %compute quadrature points and weights
    if (dimen==1)
        [x,w] = curtis_clenshaw( o );
%         [x,w] = gauss_patterson( o );
    else
        [x,w] = sparse_grid(dimen,o,method);
    end
    


    %loop over exponent combinations
    for pt = (1:porder)*uset
        for pz = (1:porder)*usez
            pass = zeros(maxorder+2);

            for py=1:porder
               for px=1:porder


                    %compute exact integral
                    p=[px,py,pz,pt];
                    coefs=1;
                    int = p+1;
                    intcoefs=1;
                    for d=1:dimen
                        intcoefs = intcoefs(1:length(coefs))./int(:,d).*x1.^(p(:,d)+1)-intcoefs(1:length(coefs))./int(:,d).*x0.^(p(:,d)+1);
                    end
                    exactint = sum(intcoefs);



                    %integrate polynomial with quadrature points
                    for i = 1:length(w)
                        testcoefs = coefs;
                        for d = 1:dimen
                            testcoefs = testcoefs.*x(i,d).^p(:,d);
                        end
                        testint(i) = sum(testcoefs)*w(i);
                    end
                    testint = sum(testint);

                    %test accuracy and store results
                    err = testint - exactint;

                    if abs(err) > 1e-12
                       pass(px,py)=0;
                    else
                       pass(px,py)=1;
                    end

                    %rule test Curtis Clenshaw
                    if (px+py+pz+pt <= o+dimen-1)
                        if pass(px,py)==0
                        rule_fail = rule_fail + 1;
                        end
                    end
                    
                end
            end
            
%     hold off
%     mesh(xx,yy,zeros(size(xx)));
%     hold on
%     for j=1:size(pass,2);
%         for i=1:size(pass,1);
%             surf(xx(i:i+1,j:j+1),yy(i:i+1,j:j+1),pass(i,j)*ones(2,2))
%             text(xx(i,j)+0.25,yy(i,j)+0.25,pass(i,j),[num2str(i),',',num2str(j)])
%         end
%     end
%     view([0,0,90])
%     title(['Quadrature Order: ',num2str(o),' - pz=',num2str(pz)])
    
%     pause
            
        end
    end
    
    
    
end

if method=='CC'
    fprintf('Rule to specify quadrature order:\n     order = sum(exponents)+1-dimension\n\n');
elseif method=='GP'
    fprintf('Rule to specify quadrature order:\n     order = max(exponents)\n\n');
end
fprintf('Number of failures: %10.0f\n',rule_fail)