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

a = 1;
dom = linspace(0,1,50);
sigma = 0.08;
deltaX = 1/50;

u = zeros(50,500);

for i = 1:50
    u(i,1) = exp((-0.5)*((dom(i) - 0.5)/sigma)^2);
end

%Explicit Euler Method       u(n+1) = u(n) + h*u'(n)
%-------------------------------------------------------------------------
courantNum = 0.1;
h = courantNum*deltaX/a;

A_col_vec = [0; -1];
A_col_vec(end+1:50) =  0;

A_row_vec = [0 1];
A_row_vec(end+1:50) =  0;

A = (1/(2*deltaX))*(toeplitz(A_col_vec, A_row_vec));

for i = 1:500
    for n = 2:50
        u(n,i) = u(n-1,i) + A(n,n)*h*u(n-1,i);
    end
end

plot(dom, u(:,1), dom, u)
title('Explicit Euler Approximation')
ylabel('u_n')
xlabel('Node')
legend('Exact', 'Approximation')
xticks([0 .1 .2 .3 .4 .5 .6 .7 .8 .9 1.0])
xticklabels({'0','.1','.2','.3','.4', '.5', '.6', '.7', '.8', '.9', '1'})

