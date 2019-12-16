function [M,d,a] = galerkin()
%solves for the coefficients of a power series approximation using the
%Galerkin weighted residuals

%create the symbolic functions for the residual and weight fn.s
syms R(x,a1, a2, a3) W1(x) W2(x) W3(x);
R = -1 + a1*(1 - x) + a2*(2*x - x.^2) + a3*(3*x.^2 - x.^3);
W1 = 1;
W2 = x;
W3 = x^2;

%solve the integral of R*W over the domain and get the coefficients
S1 = int(R*W1,x,0,1);
coeff_1 = double(coeffs(S1));
S2 = int(R*W2,x,0,1);
coeff_2 = double(coeffs(S2));
S3 = int(R*W3,x,0,1);
coeff_3 = double(coeffs(S3));

d = [-1*coeff_1(1) -1*coeff_2(1) -1*coeff_3(1)]';
M = [flip(coeff_1(2:end)); flip(coeff_2(2:end)); flip(coeff_3(2:end))];
a = linsolve(M,d);
end

