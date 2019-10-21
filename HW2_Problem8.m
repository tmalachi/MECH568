% MECH 568 HW 2 Problem 6.8
% Author: Thomas Andreano
% Version 1: 10/17/2019
%-----------------------------------------------------------------------
% Solution to the one-dimensional linear convection equation using
% 4th-order (noncompact) difference in space and RK4 time marching.

clc;
clear all;
close all;

a = 1;  %given by problem
dom = linspace(0,1,50); %create a grid of 50 points from 0 to 1
sigma = 0.08; %variable in the exact solution, given
deltaX = 1/50; %change in x determined by grid size of 50 and span of 0 to 1
courantNum = 1; %given
h = courantNum*deltaX/a; %calculate time step h

%create matrix for AM2 solution and intermediate steps
u_RK4 = zeros(50,501); 
u_hat = zeros(50,501);
u_tilde = zeros(50,501);
u_bar = zeros(50,501);

A_col_vec = [0; -8; 1];
A_col_vec(end+1:50) =  0;

A_row_vec = [0 8 -1];
A_row_vec(end+1:50) =  0;

A = (toeplitz(A_col_vec, A_row_vec));
A(50,1) = 8;
A(1,50) = -8;
A(49,1) = 1;
A(1, 49) = -1;
A = 1/(12*deltaX)*A;

for i = 1:50
    u_RK4(i,1) = exp((-0.5)*((dom(i) - 0.5)/sigma)^2); %calculate initial condition
end

%RK4 time marching approximation
for t = 2:501
    
    u_hat(:,t) = u_RK4(:, t-1) + 0.5*h*A*u_RK4(:, t-1);

    u_tilde(:,t) = u_RK4(:, t-1) +  0.5*h*A*u_hat(:,t);

    u_bar(:,t) = u_RK4(:,t-1) + h*A*u_tilde(:,t);

    u_RK4(:,t) = u_RK4(:,t-1) + (1/6)*h*(A*u_RK4(:,t-1) + 2*(A*u_hat(:,t) + A*u_tilde(:,t)) + A*u_bar(:,t));

end

plot(dom, u_RK4(:,1), dom, u_RK4(:,51), dom, u_RK4(:,501));
title('4th Order Runge-Kutta Approximation')
ylabel('u_n')
xlabel('Displacement')
legend('Exact', 'RK4 at t = 1', 'RK4 at t = 10')
xticks([0 .1 .2 .3 .4 .5 .6 .7 .8 .9 1.0])
xticklabels({'0','.1','.2','.3','.4', '.5', '.6', '.7', '.8', '.9', '1'})
