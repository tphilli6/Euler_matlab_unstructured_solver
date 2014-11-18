%% Integration by parts test
% 
% P = @(y)100000+10000*y.^1;
% U = @(y)800 + 80*y.^1;
% 
% % PU = @(y) 80000000*y + 3200000*y.^5 + (800000*y.^9)/9;
% F = @(y) P(y).*U(y).^2;
% 
% order = 2;
% 
% [x8, w8] = curtis_clenshaw(2*order);
% 
% [x, w] = curtis_clenshaw(order);
% 
% int_cc = sum(F(x).*w);
% int_ccp1 = sum(F(x8).*w8);
% 
% % int_PU_ex2 = sum(P(x8).*U(x8).*w8);

%integration by parts test
coefs = [100000, 2000, 4000, -6000, 10000
         800,     -17,   40,    20,    80
         800,     -17,   40,    20,    80
         800,     -17,   40,    20,    80
         800,     -17,   40,    20,    80];
exponents = [0, 1, 2, 3, 4
             0, 1, 2, 3, 4
             0, 1, 2, 3, 4
             0, 1, 2, 3, 4
             0, 1, 2, 3, 4];
         
% coefs = [100000, 10000
%          800,  80
%          800,  80
%          800,  80
%          800,  80];
% exponents = [0, 4
%              0, 4
%              0, 4
%              0, 4
%              0, 4];
         
int_part1 = poly_int(0,1, coefs(1:3,:), exponents(1:3,:));
int_part2 = poly_int2(0,1, coefs(1:3,:), exponents(1:3,:));
int_part3 = poly_int3(0,1, coefs(1:3,:), exponents(1:3,:));
% int_part_1eq = poly_int3(0,1, coefs(1,:), exponents(1,:));

disp([int_part,int_part2,int_part3])
disp([int_part2-int_part,int_part3-int_part])

% fprintf('Curtis Clenshaw err     %23.15g, (%23.15g %%)\n', int_cc-int_part, (int_cc-int_part)/int_part*100 );
% fprintf('Curtis Clenshaw err2    %23.15g, (%23.15g %%)\n', int_ccp1-int_part, (int_ccp1-int_part)/int_part*100 );
% fprintf('Curtis Clenshaw err est %23.15g, (%23.15g %%)\n', int_cc-int_ccp1,(int_cc-int_ccp1)/int_ccp1*100);
% err_est = int_cc-int_ccp1;
% 
% disp(err)