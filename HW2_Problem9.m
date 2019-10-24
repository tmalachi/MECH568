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
sigma = 0.08; %variable in the exact solution, given
courantNum = 1; %given

%solution with grid space of 100
%--------------------------------------------------------------------------
[u_100, error_100] = RK4_order2(1,100);
[u_200, error_200] = RK4_order2(1,200);
[u_400, error_400] = RK4_order2(1,400);

x = [100, 200, 400];
error = [error_100, error_200, error_400];

scatter(x, error);
set(gca,'yscale','log')
title('4th Order Runge-Kutta Approximation Error')
ylabel('Error')
xlabel('Number of Nodes')
xticks([100 200 400])