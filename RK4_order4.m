% RK4 method for a 4th order non compact difference
% Author: Thomas Andreano
% Version 1: 10/23/2019
%-----------------------------------------------------------------------

function [u_RK4,error] = RK4_order4(courant, nodes, time)

a = 1;  %given by problem
sigma = 0.08; %variable in the exact solution, given

%solution with grid space of 100
%--------------------------------------------------------------------------
grid = linspace(0,1,nodes); %create a grid of n points from 0 to 1
deltaX = 1/nodes; %change in x determined by grid size of 50 and span of 0 to 1
h = courant*deltaX/a; %calculate time step h

%create matrix for AM2 solution and intermediate steps
u_RK4 = zeros(nodes,(1/h + 1)); 
u_hat = zeros(nodes,(1/h + 1));
u_tilde = zeros(nodes,(1/h + 1));
u_bar = zeros(nodes,(1/h + 1));

%create matrix for approximation: u' = A*u
A_col_vec = [0; -8; 1];
A_col_vec(end+1:nodes) =  0;

A_row_vec = [0 8 -1];
A_row_vec(end+1:nodes) =  0;

A = (toeplitz(A_col_vec, A_row_vec));
A(nodes,1) = 8;
A(1,nodes) = -8;
A(nodes-1,1) = 1;
A(1, nodes-1) = -1;
A = -1/(12*deltaX)*A;

for i = 1:nodes
    u_RK4(i,1) = exp((-0.5)*((grid(i) - 0.5)/sigma)^2); %calculate initial condition
end

%RK4 time marching approximation
for t = 2:(time/h + 1)
    
    u_hat(:,t) = u_RK4(:, t-1) + 0.5*h*A*u_RK4(:, t-1);

    u_tilde(:,t) = u_RK4(:, t-1) +  0.5*h*A*u_hat(:,t);

    u_bar(:,t) = u_RK4(:,t-1) + h*A*u_tilde(:,t);

    u_RK4(:,t) = u_RK4(:,t-1) + (1/6)*h*(A*u_RK4(:,t-1) + 2*(A*u_hat(:,t) + A*u_tilde(:,t)) + A*u_bar(:,t));

end

error = sqrt(sum(((u_RK4(:,1/h + 1) - u_RK4(:,1)).^2)/nodes));

end
