% MECH 568 HW 2 Problem 6.7
% Author: Thomas Andreano
% Version 1: 10/10/2019
%-----------------------------------------------------------------------
% Solution to the one-dimensional linear convection equation using 2nd
% order centered difference. Uses explicit Euler, 2nd order
% Adams-Bashforth (AB2), implicit Euler, trapezoidal, and 4th order
% Runge-Kutta methods.

clc;
clear all;
close all;

a = 1;  %given by problem
dom = linspace(0,1,50); %create a grid of 50 points from 0 to 1
sigma = 0.08; %variable in the exact solution, given
deltaX = 1/50;

%Explicit Euler Method       u(n+1) = u(n) + h*u'(n)
%-------------------------------------------------------------------------
courantNum = 0.1;
h = courantNum*deltaX/a;

u_EE = zeros(50,501); %create matrix for EE solution

for i = 1:50
    u_EE(i,1) = exp((-0.5)*((dom(i) - 0.5)/sigma)^2); %calculate initial condition
end

%create matrix for approximation: u' = A*u
A_col_vec = [0; -1];
A_col_vec(end+1:50) =  0;

A_row_vec = [0 1];
A_row_vec(end+1:50) =  0;

A = (toeplitz(A_col_vec, A_row_vec)); 

%time-march solution
for t = 2:501
    u_EE(:,t) = u_EE(:,t-1) + h*A*u_EE(:,t-1); 
end

%2nd order Adams-Bashforth (AB2)      u_n+1 = u_n + .5*h*[3*u'_n - u'_n-1]
%-------------------------------------------------------------------------
u_AB2 = zeros(50,501); %create matrix for AB2 solution

for i = 1:50
    u_AB2(i,1) = exp((-0.5)*((dom(i) - 0.5)/sigma)^2); %calculate initial condition
end

%use explicit Euler method to approximate 2nd time-step solution
u_AB2(:,2) = u_AB2(:,1) + h*A*u_AB2(:,1);

%proceed with AB2 approximation
for t = 3:501
    u_AB2(:,t) = u_AB2(:,t-1) + 0.5*h*(3*A*u_AB2(:,t-1) - A*u_AB2(:,t-2));
end


%Implicit Euler Method          u_n+1 = u_n + h*u'_n+1
%-------------------------------------------------------------------------
courantNum = 1;
h = courantNum*deltaX/a;

u_IE = zeros(50,501); %create matrix for IE solution

for i = 1:50
    u_IE(i,1) = exp((-0.5)*((dom(i) - 0.5)/sigma)^2); %calculate initial condition
end

for t = 2:501
    u_IE(:,t) = inv(eye(50) - h*A)*u_IE(:,t-1);
end



%Trapezoidal (AM2)
%-------------------------------------------------------------------------
u_AM2 = zeros(50,501); %create matrix for AM2 solution

for i = 1:50
    u_AM2(i,1) = exp((-0.5)*((dom(i) - 0.5)/sigma)^2); %calculate initial condition
end

for t = 2:501
    u_AM2(:,t) = inv(eye(50) - 0.5*h*A)*((u_AM2(:,t-1)*(eye(501) + 0.5*h*A)));
end


%4th Order Runge-Kutta (RK4)
%-------------------------------------------------------------------------

plot(dom, u_EE(:,1), dom, u_EE(:,501), dom, u_AB2(:,501), '--',...
    dom, u_IE(:,501), '.', dom, u_AM2(:, 501));
title('AB2 Approximation')
ylabel('u_n')
xlabel('Displacement')
legend('Exact', 'EE', 'AB2')
xticks([0 .1 .2 .3 .4 .5 .6 .7 .8 .9 1.0])
xticklabels({'0','.1','.2','.3','.4', '.5', '.6', '.7', '.8', '.9', '1'})
