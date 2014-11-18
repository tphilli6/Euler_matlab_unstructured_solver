%% Sparse grid test
% Integrates a polynomial using a multidimensional quadrature for an 
% arbitrary number of dimensions.
% Polynomial form: P = 1 + x + x^2 + y + y*x + y*x^2 + y^2 + y^2*x + y^2*x^2
% Author: Tyrone Phillips

clc
clear all

%% Higher dimension quadrature test
dimen=1;
x1=1;
x0=0;
for o = 1:20

    %compute appropriate polynomial order for given quadrature order
    porder = o;
    if (mod(o,2)==0); porder = o-1; end %not necessary but included for consistency with Fortran version

    %set polynomial coefficients
    coefs = ones((porder+1).^dimen,1);
    
    %compute exponents
    n_unknowns = (porder+1)^dimen;
    for d=1:dimen
       spacing = (porder+1)^(d-1);
       for i = 1:porder+1:n_unknowns/spacing
          for j=0:porder
             p((i-1+j)*spacing+1:(i+j)*spacing,d) = j; 
          end
       end
    end

    %compute exact integral
    int = p+1;
    intcoefs=coefs;
    for d=1:dimen
        intcoefs = intcoefs(1:length(coefs))./int(:,d).*x1.^(p(:,d)+1)-intcoefs(1:length(coefs))./int(:,d).*x0.^(p(:,d)+1);
    end
    exactint = sum(intcoefs);

    %compute quadrature points and weights
    if (dimen==1)
        [x,w] = curtis_clenshaw( o );
    else
        [x,w] = sparse_grid(dimen,o);
    end

    %integrate polynomial with quadrature points
    for i = 1:length(w)
        testcoefs = coefs;
        for d = 1:dimen
            testcoefs = testcoefs.*x(i,d).^p(:,d);
        end
        testint(i) = sum(testcoefs)*w(i);
    end
    testint = sum(testint);

    %test accuracy and print results
    err = testint - exactint;

    if abs(err) > 1e-12
       fprintf(' Test failed at %2.0f  ',o)
       fprintf(' Exact: %23.15e Test:  %23.15e Error: %23.15e\n',exactint, testint,err);
    else
        fprintf(' Passed at %2.0f \n',o)
    end

end