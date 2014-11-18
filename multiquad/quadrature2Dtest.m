% This code tests multi-dimensional quadratures to determine the combination
% of monomials that can be exactly integrated.
% 2D only.
% Polynomial form: P = x^px*y^py
% Author: Tyrone Phillips

clear all
clc

%% Higher dimension quadrature test
dimen=2; %do not change
maxorder = 10;
x1=1;
x0=0;
rule_fail = 0;

for o = 1:maxorder

    %compute appropriate polynomial order for given quadrature order
    porder = o;
%     if (mod(o,2)==0); porder = o-1; end

    %compute quadrature points and weights
    if (dimen==1)
        [x,w] = curtis_clenshaw( o );
    else
        [x,w] = sparse_grid(dimen,o,'CC');
    end
    
    pass = zeros(maxorder+2);
    [xx,yy] = meshgrid([0.5:maxorder+2.5],[0.5:maxorder+2.5]);
    
    %loop over exponent combinations
    for py=1:porder+1
        for px=1:porder+1
            
            
            %compute exact integral
            p=[px-1,py-1];
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
            
            %rule test
            if (sum(p) <= o)
                if pass(px,py)==0
                rule_fail = rule_fail + 1;
                end
            end
            
        end
    end

    
    
    hold off
    mesh(xx,yy,zeros(size(xx)));
    hold on
    for j=1:size(pass,2);
        for i=1:size(pass,1);
            surf(xx(i:i+1,j:j+1),yy(i:i+1,j:j+1),pass(i,j)*ones(2,2))
            text(xx(i,j)+0.25,yy(i,j)+0.25,pass(i,j),[num2str(i-1),',',num2str(j-1)])
        end
    end
    view([0,0,90])
    title(['Quadrature Order: ',num2str(o)])
    
    pause
    
end

fprintf('Rule to specify quadrature order:\n     order = sum(exponents)\n\n');
fprintf('Number of failures: %10.0f\n',rule_fail)