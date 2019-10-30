% MECH 568 HW 2 Problem 6.8
% Author: Thomas Andreano
% Version 1: 10/17/2019
%-----------------------------------------------------------------------
% Solution to the one-dimensional linear convection equation using
% 4th-order (noncompact) difference in space and RK4 time marching.

clc;
clear all;
close all;

nodes = 50;
time = 10;
courant = 1;
grid = linspace(0,1,50);
[u_RK4, error] = RK4_order4(courant,nodes,time);

plot(grid, u_RK4(:,1), grid, u_RK4(:,51), grid, u_RK4(:,501));
title('RK4 Appoximation with 4th order spacial differencing')
ylabel('u_n')
xlabel('Displacement')
legend('Exact', 'RK4 at t = 1', 'RK4 at t = 10')
xticks([0 .1 .2 .3 .4 .5 .6 .7 .8 .9 1.0])
xticklabels({'0','.1','.2','.3','.4', '.5', '.6', '.7', '.8', '.9', '1'})