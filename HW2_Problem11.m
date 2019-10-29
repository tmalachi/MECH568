%MECH 568 HW 2 Problem 6.11
%Author: Thomas Andreano
%Version 1: 10/24/2019
%-------------------------------------------------------------------------
% Solution to the one-dimensional linear convection equation using 2nd
% order centered difference and 4th order Runge-Kutta

clear all;
close all;
clc;

a = 1;  %given by problem
nodes = 100;
dom = linspace(0,1,nodes); %create a grid from 0 to 1
deltaX = 1/100; %change in x determined by grid size of 50 and span of 0 to 1
omega = 10*pi; %given by problem
h = deltaX; %based on a courant number of unity
time = 1000;

%create matrix for AM2 solution and intermediate steps
u_RK4 = zeros(nodes,time/h+1); 
u_hat = zeros(nodes,time/h+1);
u_tilde = zeros(nodes,time/h+1);
u_bar = zeros(nodes,time/h+1);
bc = zeros(nodes,1);
bc(1) = 1/(2*deltaX);
u_exact = zeros(nodes, time/h+1);

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


approx = animatedline();
exact = animatedline();

for t = 2:time/h+1
    
    %RK4 time marching approximation
    u_hat(:,t) = u_RK4(:,t-1) + 0.5*h*(A*u_RK4(:,t-1) + sin(omega*t*h)*bc);

    u_tilde(:,t) = u_RK4(:,t-1) +  0.5*h*(A*u_hat(:,t) + sin(omega*t*h)*bc);

    u_bar(:,t) = u_RK4(:,t-1) + h*(A*u_tilde(:,t) + sin(omega*t*h)*bc);

    u_RK4(:,t) = u_RK4(:,t-1) + (1/6)*h*((A*u_RK4(:, t-1) + sin(omega*t*h)*bc) + (2*((A*u_hat(:,t)+ sin(omega*t*h)*bc) + (A*u_tilde(:,t) + sin(omega*t*h)*bc)) + (A*u_bar(:,t) + sin(omega*t*h))));
    
    %Exact solution
    for i = 1:nodes
    u_exact(i,t) = sin(omega*(t*h - i*deltaX));
    addpoints(approx, u_RK4(i,t));
    addpoints(exact, u_exact(i,t));
    drawnow
    end
    
    
    
    title('4th Order Runge-Kutta Approximation')
    ylabel('u_n')   
    xlabel('Displacement')
    legend('Exact', 'RK4')
    xticks([0 .1 .2 .3 .4 .5 .6 .7 .8 .9 1.0])
    xticklabels({'0','.1','.2','.3','.4', '.5', '.6', '.7', '.8', '.9', '1'})

end


% %plot results
% plot(dom, u_exact, dom, u_RK4(:,time/h+1), '--');
% title('4th Order Runge-Kutta Approximation')
% ylabel('u_n')
% xlabel('Displacement')
% legend('Exact', 'RK4')
% xticks([0 .1 .2 .3 .4 .5 .6 .7 .8 .9 1.0])
% xticklabels({'0','.1','.2','.3','.4', '.5', '.6', '.7', '.8', '.9', '1'})
