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

[u_100, error_100, exact_100] = RK4_6_11(100);
[u_200, error_200, exact_200] = RK4_6_11(200);
[u_400, error_400, exact_400] = RK4_6_11(400);

x = [100, 200, 400];
error = [error_100, error_200, error_400];

p = polyfit(x, log10(error_100), 1);
z = polyval(p, log10(x));

%plot results
figure(1)
plot(linspace(0,1,100), exact_100(:,end), linspace(0,1,100), u_100(:,end), '--');
title('4th Order Runge-Kutta Approximation 100 nodes')
ylabel('u_n')
xlabel('Displacement')
legend('Exact', 'RK4')
xticks([0 .1 .2 .3 .4 .5 .6 .7 .8 .9 1.0])
xticklabels({'0','.1','.2','.3','.4', '.5', '.6', '.7', '.8', '.9', '1'})

figure(2)
plot(linspace(0,1,200), exact_200(:,end), linspace(0,1,200), u_200(:,end), '--');
title('4th Order Runge-Kutta Approximation 200 nodes')
ylabel('u_n')
xlabel('Displacement')
legend('Exact', 'RK4')
xticks([0 .1 .2 .3 .4 .5 .6 .7 .8 .9 1.0])
xticklabels({'0','.1','.2','.3','.4', '.5', '.6', '.7', '.8', '.9', '1'})

figure(3)
plot(linspace(0,1,400), exact_400(:,end), linspace(0,1,400), u_400(:,end), '--');
title('4th Order Runge-Kutta Approximation 400 nodes')
ylabel('u_n')
xlabel('Displacement')
legend('Exact', 'RK4')
xticks([0 .1 .2 .3 .4 .5 .6 .7 .8 .9 1.0])
xticklabels({'0','.1','.2','.3','.4', '.5', '.6', '.7', '.8', '.9', '1'})

figure(4)
scatter(x, error);
set(gca,'yscale','log')
title('Error of RK4 with 2nd order spacial differencing')
ylabel('Error')
xlabel('Number of Nodes')
xticks([100 200 400])