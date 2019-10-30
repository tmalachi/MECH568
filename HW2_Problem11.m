%MECH 568 HW 2 Problem 6.11
%Author: Thomas Andreano
%Version 1: 10/24/2019
%-------------------------------------------------------------------------
% Solution to the one-dimensional linear convection equation using 2nd
% order centered difference and 4th order Runge-Kutta with inflow-outflow
% boundary conditions

clear all;
close all;
clc;

a = 1;  %given by problem
nodes = 100;
dom_100 = linspace(0,1,nodes); %create a grid from 0 to 1
deltaX = 1/100; %change in x determined by grid size of 50 and span of 0 to 1
omega = 10*pi; %given by problem
h = deltaX; %based on a courant number of unity

%create matrices for solution and intermediate steps
u_RK4 = zeros(nodes,1); 
u_hat = zeros(nodes,1);
u_tilde = zeros(nodes,1);
u_bar = zeros(nodes,1);
bc = zeros(nodes,1);
bc(1) = 1/(2*deltaX);
u_exact = zeros(nodes,1);
error = zeros(1,1);

%Solve inflow boundary (node 0)
for i = 1:nodes
    u_RK4(:,1) = 1;
end

%create matrix for approximation: u' = A*u + (bc)
A_col_vec = [0; -1];
A_col_vec(end+1:nodes) =  0;

A_row_vec = [0 1];
A_row_vec(end+1:nodes) =  0;

A = (toeplitz(A_col_vec, A_row_vec));
A(nodes,nodes) = 2;
A(nodes,nodes -1) = -2;
A = -1/(2*deltaX)*A;

t = 2;
ss_cond = false;

%RK4 time marching approximation 
while ss_cond == false
    
    u_hat(:,t) = u_RK4(:,t-1) + 0.5*h*(A*u_RK4(:,t-1) + sin(omega*t*h)*bc);

    u_tilde(:,t) = u_RK4(:,t-1) +  0.5*h*(A*u_hat(:,t) + sin(omega*t*h)*bc);

    u_bar(:,t) = u_RK4(:,t-1) + h*(A*u_tilde(:,t) + sin(omega*t*h)*bc);

    u_RK4(:,t) = u_RK4(:,t-1) + (1/6)*h*((A*u_RK4(:, t-1) + sin(omega*t*h)*bc) + (2*((A*u_hat(:,t)+ sin(omega*t*h)*bc) + (A*u_tilde(:,t) + sin(omega*t*h)*bc)) + (A*u_bar(:,t) + sin(omega*t*h))));
    
    %calculate the exact solution
    for i = 1:nodes
        u_exact(i,t) = sin(omega*(t*h - i*deltaX));
    end
    
    %calculate the error
    error_100 = sqrt(sum(((u_RK4(:,t) - u_exact(:,t)).^2)/nodes));
    
    %check if the solution is significantly different than it was 2 seconds
    %before
    if t > (2/h)
       u_diff = u_RK4(:,t) - u_RK4(:,t - 2/h);
       
       B = abs(u_diff) < 1e-12;
       ss_cond = all(B);
    end
    
    t = t + 1;      
        
end

error(t-1) = error_100;

t = 2;
ss_cond = false;

%RK4 time marching approximation 
while ss_cond == false
    
    u_hat(:,t) = u_RK4(:,t-1) + 0.5*h*(A*u_RK4(:,t-1) + sin(omega*t*h)*bc);

    u_tilde(:,t) = u_RK4(:,t-1) +  0.5*h*(A*u_hat(:,t) + sin(omega*t*h)*bc);

    u_bar(:,t) = u_RK4(:,t-1) + h*(A*u_tilde(:,t) + sin(omega*t*h)*bc);

    u_RK4(:,t) = u_RK4(:,t-1) + (1/6)*h*((A*u_RK4(:, t-1) + sin(omega*t*h)*bc) + (2*((A*u_hat(:,t)+ sin(omega*t*h)*bc) + (A*u_tilde(:,t) + sin(omega*t*h)*bc)) + (A*u_bar(:,t) + sin(omega*t*h))));
    
    %calculate the exact solution
    for i = 1:nodes
        u_exact(i,t) = sin(omega*(t*h - i*deltaX));
    end
    
    %calculate the error
    error_200 = sqrt(sum(((u_RK4(:,t) - u_exact(:,t)).^2)/nodes));
    
    %check if the solution is significantly different than it was 2 seconds
    %before
    if t > (2/h)
       u_diff = u_RK4(:,t) - u_RK4(:,t - 2/h);
       
       B = abs(u_diff) < 1e-12;
       ss_cond = all(B);
    end
    
    t = t + 1;      
        
end


%RK4 time marching approximation 
while ss_cond == false
    
    u_hat(:,t) = u_RK4(:,t-1) + 0.5*h*(A*u_RK4(:,t-1) + sin(omega*t*h)*bc);

    u_tilde(:,t) = u_RK4(:,t-1) +  0.5*h*(A*u_hat(:,t) + sin(omega*t*h)*bc);

    u_bar(:,t) = u_RK4(:,t-1) + h*(A*u_tilde(:,t) + sin(omega*t*h)*bc);

    u_RK4(:,t) = u_RK4(:,t-1) + (1/6)*h*((A*u_RK4(:, t-1) + sin(omega*t*h)*bc) + (2*((A*u_hat(:,t)+ sin(omega*t*h)*bc) + (A*u_tilde(:,t) + sin(omega*t*h)*bc)) + (A*u_bar(:,t) + sin(omega*t*h))));
    
    %calculate the exact solution
    for i = 1:nodes
        u_exact(i,t) = sin(omega*(t*h - i*deltaX));
    end
    
    %calculate the error
    error_400 = sqrt(sum(((u_RK4(:,t) - u_exact(:,t)).^2)/nodes));
    
    %check if the solution is significantly different than it was 2 seconds
    %before
    if t > (2/h)
       u_diff = u_RK4(:,t) - u_RK4(:,t - 2/h);
       
       B = abs(u_diff) < 1e-12;
       ss_cond = all(B);
    end
    
    t = t + 1;      
        
end


x = [100, 200, 400];
error = [error_100, error_200, error_400];

%plot results
figure(1)
plot(dom, u_exact(:, t-1), dom, u_RK4(:,t-1), '--');
title('4th Order Runge-Kutta Approximation')
ylabel('u_n')
xlabel('Displacement')
legend('Exact', 'RK4')
xticks([0 .1 .2 .3 .4 .5 .6 .7 .8 .9 1.0])
xticklabels({'0','.1','.2','.3','.4', '.5', '.6', '.7', '.8', '.9', '1'})

figure(2)
plot(linspace(0,t, numel(error(1,:))), error(1,:));
xticks([0 time])
xticklabels({'0',time})

figure(3)
scatter(x, error);
set(gca,'yscale','log')
title('Error of RK4 with 2nd order spacial differencing')
ylabel('Error')
xlabel('Number of Nodes')
xticks([100 200 400])