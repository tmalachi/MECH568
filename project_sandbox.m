clc;
clear all;
close all;

% This script performs an approximation of the differential equation
%
%                      dy/dx - y = 0
%
% over the domain 0 < x < 1 with a grid size of N = 3
%% 

N = 3;
x = linspace(0, 1, N);

userPrompt = sprintf('%s\n%s \n%s \n%s \n%s \n%s\n', 'Please enter the number corresponding to the method you would like to use',...
    '1: Least Squares',...
    '2: Galerkin',...
    '3: Collocation', ...
    '4: All methods w/ comparison');

answer = inputdlg(userPrompt);      %user inputs desired method
method = str2double(answer{1});

%switch case for chosen method
if method < 4
    switch method
        case 1
            [M,d,a] = leastSquares();
        case 2
            [M,d,a] = galerkin();
        case 3
            [M,d,a] = collocation();
    end
    
    y_approx = ones(1,N);
    for i = 1:N
        y_approx(i) = 1 + a(1)*x(i) + a(2)*x(i)^2 + a(3)* x(i)^3;
    end
    
    x_exact = linspace(0,1);
    y = exp(x_exact);
    
    plot(x_exact,y,x,y_approx)
    legend('Exact', 'Approximation')
    xlabel('x')
    ylabel('y')
     
elseif method == 4      %choice 4 is a comparison of all 3 methods
        
    [M1,d1,a1] = leastSquares();
    [M2,d2,a2] = galerkin();
    [M3,d3,a3] = collocation();

    y_LS = ones(1,N);
    y_G = ones(1,N);
    y_C = ones(1,N);
    y = ones(1,N);
    
    for i = 1:N
        y_LS(i) = double(1 + a1(1)*x(i) + a1(2)*x(i)^2 + a1(3)* x(i)^3);
        y_G(i) = double(1 + a2(1)*x(i) + a2(2)*x(i)^2 + a2(3)* x(i)^3);
        y_C(i) = double(1 + a3(1)*x(i) + a3(2)*x(i)^2 + a3(3)* x(i)^3);
        y(i) = exp(x(i));
    end
    
    er_LS = norm(y - y_LS);
    er_G = norm(y - y_G);
    er_C = norm(y - y_C);    
    
    fprintf('              Comparison of Polynomial Coefficients')
    fprintf('\n---------------------------------------------------------------')
    fprintf('\n|     Method     |      a1      |      a2      |      a3      |')
    fprintf('\n|  Least Squares | %d | %d | %d |', a1(1), a1(2), a1(3))
    fprintf('\n|    Galerkin    | %d | %d | %d |', a2(1), a2(2), a2(3))
    fprintf('\n|   Collocation  |      %d       | %d | %d |', a3(1), a3(2), a3(3))
    fprintf('\n---------------------------------------------------------------\n')
    
    fprintf('\n\n                   Comparison of Approximate Solutions')
    fprintf('\n-------------------------------------------------------------------------------')
    fprintf('\n|    x    |  Least Squares  |    Galerkin    |  Collocation   |     Exact     |')
    fprintf('\n|    0    |        %d        |        %d       |       %d        |       %d       |', y_LS(1), y_G(1), y_C(1), y(1))
    fprintf('\n|   0.5   |  %d   |  %d  |  %d  | %d  |', y_LS(2), y_G(2), y_C(2), y(2))
    fprintf('\n|   1.0   |  %d   |  %d  |  %d  | %d  |', y_LS(3), y_G(3), y_C(3), y(3))
    fprintf('\n|  L2 err |  %d   |  %d  |  %d  |               |', er_LS, er_G, er_C)
    fprintf('\n-------------------------------------------------------------------------------\n')

else
    disp('Invalid input')
end
