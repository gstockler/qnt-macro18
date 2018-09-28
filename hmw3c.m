%% Chebyshev interpolation nodes and Chebyshev polynomial approximation
clear all;
clc

%% f=exp(1/x)
syms x
func=exp(1./x);

%order-3
[poly3]=cheb_poly(func,-1,1,3,4);
%error
error_3=abs(func-poly3);

%order-5
[poly5]=cheb_poly(func,-1,1,5,6);
%error
error_5=abs(func-poly5);

%order-10
[poly10]=cheb_poly(func,-1,1,10,11);
%error
error_10=abs(func-poly10);

%PLOTS
figure(1)
fplot(func,[-1 1],'b');
title('Approximation of $f=\exp{\frac{1}{x}}$ for $x\in [-1,1]$','interpreter','latex')
hold on
figure(1)
fplot(poly3,[-1 1],'r');
hold on
figure(1)
fplot(poly5,[-1 1],'y');
hold on
figure(1)
fplot(poly10,[-1 1],'m');
ylim([0 2]);
legend('$f=\exp{\frac{1}{x}}$', '3-cheb', '5-cheb', '10-cheb')
set(legend,'Interpreter','latex')
legend('show')
print -dpdf ex3c11.eps
hold off

figure(2)
fplot(error_3,[-1 1]);
title('Approx. Errors of $f=\exp{\frac{1}{x}}$ for $x\in [-1,1]$','interpreter','latex')
hold on
figure(2)
fplot(error_5,[-1 1]);
hold on
figure(2)
fplot(error_10,[-1 1]);
%ylim([-2 2]);
legend('3-cheb', '5-cheb', '10-cheb')
set(legend,'Interpreter','latex')
legend('show')
print -dpdf ex3c12.eps
hold off

%% runge function
clear all;
clc
syms x
func=1/(1+25.*x^2);

%order-3
[poly3]=cheb_poly(func,-1,1,3,4);
%error
error_3=abs(func-poly3);

%order-5
[poly5]=cheb_poly(func,-1,1,5,6);
%error
error_5=abs(func-poly5);

%order-10
[poly10]=cheb_poly(func,-1,1,10,11);
%error
error_10=abs(func-poly10);

%PLOTS
figure(1)
fplot(func,[-1 1],'b');
title('Approximation of $f=\frac{1}{1+25x^2}$ for $x\in [-1,1]$','interpreter','latex')
hold on
figure(1)
fplot(poly3,[-1 1],'r');
hold on
figure(1)
fplot(poly5,[-1 1],'y');
hold on
figure(1)
fplot(poly10,[-1 1],'m');
legend('$f=\frac{1}{1+25x^2}$', '3-cheb', '5-cheb', '10-cheb')
set(legend,'Interpreter','latex')
legend('show')
print -dpdf ex3c21.eps
hold off

figure(2)
fplot(error_3,[-1 1]);
title('Approx. Errors of $f=\frac{1}{1+25x^2}$ for $x\in [-1,1]$','interpreter','latex')
hold on
figure(2)
fplot(error_5,[-1 1]);
hold on
figure(2)
fplot(error_10,[-1 1]);
legend('3-cheb', '5-cheb', '10-cheb')
set(legend,'Interpreter','latex')
legend('show')
print -dpdf ex3c22.eps
hold off

%% ramp function
clear all;
clc
syms x
func=(x+abs(x))./2;

%order-3
[poly3]=cheb_poly(func,-1,1,3,4);
%error
error_3=abs(func-poly3);

%order-5
[poly5]=cheb_poly(func,-1,1,5,6);
%error
error_5=abs(func-poly5);

%order-10
[poly10]=cheb_poly(func,-1,1,10,11);
%error
error_10=abs(func-poly10);

%PLOTS
figure(1)
fplot(func,[-1 1],'b');
title('Approximation of $f=\frac{1}{1+25x^2}$ for $x\in [-1,1]$','interpreter','latex')
hold on
figure(1)
fplot(poly3,[-1 1],'r');
hold on
figure(1)
fplot(poly5,[-1 1],'y');
hold on
figure(1)
fplot(poly10,[-1 1],'m');
legend('$f=\frac{1}{1+25x^2}$', '3-cheb', '5-cheb', '10-cheb')
set(legend,'Interpreter','latex')
legend('show')
print -dpdf ex3c31.eps
hold off

figure(2)
fplot(error_3,[-1 1]);
title('Approx. Errors of $f=\frac{1}{1+25x^2}$ for $x\in [-1,1]$','interpreter','latex')
hold on
figure(2)
fplot(error_5,[-1 1]);
hold on
figure(2)
fplot(error_10,[-1 1]);
legend('3-cheb', '5-cheb', '10-cheb')
set(legend,'Interpreter','latex')
legend('show')
print -dpdf ex3c32.eps
hold off