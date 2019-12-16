function [M,d,a] = collocation()
%solves for the coefficients of a power series approximation using the
%collocation method which sets R(x_i) = 0

N = 3;
x = linspace(0,1,N);

%create symbolic function of residual at each node
syms R(a1, a2, a3);
Res1 = -1 + a1*(1 - x(1)) + a2*(2*x(1) - x(1)^2) + a3*(3*x(1)^2 - x(1)^3);
Res2 = -1 + a1*(1 - x(2)) + a2*(2*x(2) - x(2)^2) + a3*(3*x(2)^2 - x(2)^3);
Res3 = -1 + a1*(1 - x(3)) + a2*(2*x(3) - x(3)^2) + a3*(3*x(3)^2 - x(3)^3);

%get the coefficients
coeff_1 = double(coeffs(Res1));
coeff_1(3) = 0;
coeff_1(4) = 0;
coeff_2 = double(coeffs(Res2));
coeff_3 = double(coeffs(Res3));
coeff_3(4) = 0;

d = [-1*coeff_1(1) -1*coeff_2(1) -1*coeff_3(1)]';
M = [coeff_1(2:end); flip(coeff_2(2:end)); flip(coeff_3(2:end))];
a = linsolve(M,d);
end

