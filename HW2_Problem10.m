% MECH 568 HW 2 Problem 6.10
% Author: Thomas Andreano
% Version 1: 10/23/2019
%-----------------------------------------------------------------------
% Error of the 4th order differencing scheme with RK4 method
% at 100, 200, and 400 nodes.

clc;
clear all;
close all;

a = 1;  %given by problem
sigma = 0.08; %variable in the exact solution, given
courantNum = 1; %given

%solution with grid space of 100
%--------------------------------------------------------------------------
[u_100, error_100] = RK4_order4(1,100,1);
[u_200, error_200] = RK4_order4(1,200,1);
[u_400, error_400] = RK4_order4(1,400,1);

x = [100, 200, 400];
error = [error_100, error_200, error_400];

scatter(x, error);
set(gca,'yscale','log')
title('4th Order Runge-Kutta Approximation Error')
ylabel('Error')
xlabel('Number of Nodes')
xticks([100 200 400]) 
