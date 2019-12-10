% MECH 568 HW #4 Problem 10.2
% Author: Thomas Andreano
% 12/07/2019
%--------------------------------------------------------------------------
%
% Solves the following equation using 2nd order differences and a
% four grid multigrid method
%
%               d2u/dx2 - 6x = 0
%
% Domain: 0 <= x <= 1
% Boundary conditions: u(0) = 0, u(1) = 1
% Initial condition: u(x) = 0
% Grid size: M = 39
%
% Gauss-Seidel and the 3-step Richardson methiod are used.
%--------------------------------------------------------------------------

clear all;
close all;
clc;

nodes = 39;
dom = linspace(0, 1, nodes)';
deltaX = (dom(end) - dom(1))/(nodes);

bc = zeros(nodes,1);
bc(end) = 1;
g = 6*dom;

%create A matrix for approximation: u' = A*u
A_col_vec = zeros(nodes,1);
A_col_vec(1) = -2;
A_col_vec(2) = 1;
A_row_vec = zeros(1,nodes);
A_row_vec(1) = -2;
A_row_vec(2) = 1;
A = (toeplitz(A_col_vec, A_row_vec));

%create A2 matrix
A2_col_vec = zeros((nodes-1)/2,1);
A2_col_vec(1) = -2;
A2_col_vec(2) = 1;
A2_row_vec = zeros(1,(nodes-1)/2);
A2_row_vec(1) = -2;
A2_row_vec(2) = 1;
A2 = (toeplitz(A2_col_vec, A2_row_vec));

I = eye(nodes);
f = deltaX^2*g - bc;

tol = 1e-6;

%% Gauss Seidel
%--------------------------------------------------------------------------

H_GS = -(tril(A));
u_GS = zeros(nodes, 1);
convGS_th = max(eig(I + H_GS\A));
G1 = I + H_GS\A;

r = inf;
count = 0;
save_count = 1;

R_12 = zeros((nodes-1)/2, nodes);
for j = 1:(nodes-1)/2
   R_12(j, j*2) = 1;
end  
   
I_21 = zeros(nodes, (nodes-1)/2);
for j = 1:(nodes-1)/2
   I_21(j+j-1,j) =  1/2;
   I_21(j*2, j) = 1;
   I_21(j*2 + 1, j) = 1/2;
end

iter_mat = (I - I_21*(A2\R_12)*A)*G1;

while r > tol
   
    u_GS = iter_mat*u_GS - iter_mat*(A\f) + A\f;

    if count > 0
        r_old = r;
        r_new = max(abs(A*u_GS - f));

        if r_new < r_set/10^(save_count)
            u_save(:,:,save_count) = u_GS;
            save_count = save_count + 1;
        end

        r = r_new;
    else
        r_set = max(abs(A*u_GS - f));
        r = r_set;
    end
   
    L2_GS(count+1) = vecnorm(r,2);
    count = count + 1;

end

dom_final = linspace(0, 1, nodes+1);
u_GS(end + 1) = 1;
u_GS_1 = u_save(:,:,2);
u_GS_1(end + 1) = 1;
u_GS_2 = u_save(:,:,3);
u_GS_2(end + 1) = 1;
u_GS_3 = u_save(:,:,4);
u_GS_3(end + 1) = 1;

figure(1)
plot(dom_final,u_GS_1, dom_final, u_GS_2, dom_final, u_GS_3,...
    dom_final, u_GS, '--')
title('Gauss Seidel Method')
legend('Solution after 100x residual decrease',...
    'Solution after 1000x residual decrease',...
    'Solution after 10000x residual decrease', 'Final solution')

figure(2)
plot(linspace(1, count, count), log(L2_GS))
title('L2 Norm of Residual for Gauss Seidel')

%% 3-step Richardson
%--------------------------------------------------------------------------
H_R = -diag(diag(A));
u_R = zeros(nodes, 1);
N = 3;
convR_th = max(eig(I + H_R\A));
G1 = I + H_R\A;

for n = 1:3
    h(n) = 1/(.5*(-min(eig(A)) - max(eig(A)) + (min(eig(A)) - ...
        max(eig(A)))*cos((2*n - 1)*pi/(2*N))));
end

r = inf;
count = 0;
save_count = 1;
while r > tol
   
    for i = 1:3
         u_R = h(i)*(iter_mat*u_R - iter_mat*(A\f) + A\f);
    end
   
   
    if count > 0
        r_old = r;
        r_new = max(abs(A*u_R - f));

        if r_new < r_set/10^(save_count)
            u_save(:,:,save_count) = u_R;
            save_count = save_count + 1;
        end
       
         r = r_new;
    else
        r_set = max(abs(A*u_R - f));
        r = r_set;
    end
   
    L2_R(count+1) = vecnorm(r,2);
    count = count + 1;

end

u_R(end + 1) = 1;
u_R_1 = u_save(:,:,2);
u_R_1(end + 1) = 1;
u_R_2 = u_save(:,:,3);
u_R_2(end + 1) = 1;
u_R_3 = u_save(:,:,4);
u_R_3(end + 1) = 1;

figure(7)
plot(dom_final,u_R_1, dom_final, u_R_2, dom_final, u_R_3,...
    dom_final, u_R, '--')
title('Richardson Method')
legend('Solution after 100x residual decrease',...
    'Solution after 1000x residual decrease',...
    'Solution after 10000x residual decrease', 'Final solution')

figure(8)
plot(linspace(1, count, count), log(L2_R))
title('L2 Norm of Residual for Richardson Method')