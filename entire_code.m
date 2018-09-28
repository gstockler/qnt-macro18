%% Question 1
clear all;
%clc;

syms x
f=x.^(321/1000);
fplot(f, [0 4]);
hold on
title('Taylor series approximation of x^{.321} around x=1')

[TS_Approx1, TS_plot1] = taylorSeries(f,1,1,0,4);
hold on
[TS_Approx2, TS_plot2] = taylorSeries(f,1,2,0,4);
hold on
[TS_Approx5, TS_plot5] = taylorSeries(f,1,5,0,4);
hold on
[TS_Approx20, TS_plot20] = taylorSeries(f,1,20,0,4);
ylim([-5 5]);
legend('f=x^{.321}', 'f-approx order 1', 'f-approx order 2', 'f-approx order 5', 'f-approx order 20')

hold off
%print -dpdf ex1.eps

%% Question 2
clear all;
%clc;

syms x
f=(x+abs(x))./2;
fplot(f, [-2 6]);
hold on
title('Taylor series approximation of $f=\frac{x+|x|}{2}$ around x=2','interpreter','latex')

[TS_Approx1, TS_plot1] = taylorSeries(f,2,1,-2,6);
hold on
[TS_Approx2, TS_plot2] = taylorSeries(f,2,2,-2,6);
hold on
[TS_Approx5, TS_plot5] = taylorSeries(f,2,5,-2,6);
hold on
[TS_Approx20, TS_plot20] = taylorSeries(f,2,20,-2,6);
ylim([-2 6]);
legend('$f=\frac{x+|x|}{2}$', 'f-approx order 1', 'f-approx order 2', 'f-approx order 5', 'f-approx order 20')
set(legend,'Interpreter','latex');
hold off
%print -dpdf ex2.eps

%% Question 3 - part a

%% Evenly spaced interpolation nodes and polynomial approximation
clc;
clear all;

%% f=exp(1/x)
syms x a0 a1 a2 a3 
f=exp(1./x);

%cubic polynomial
p=poly2sym([a0 a1 a2 a3]);
x=interpolation_n(-1,1,4);
p_ev=eval(p);
f_ev=eval(f);
H=f_ev==p_ev;
[a0,a1,a2,a3]=solve(H,[a0,a1,a2,a3]);
p=poly2sym([a0 a1 a2 a3]);
%error
error_3=abs(f-p);

%5-monomial
syms x a0 a1 a2 a3 a4 a5
f=exp(1./x);
p2=poly2sym([a0 a1 a2 a3 a4 a5]);
x=interpolation_n(-1,1,6);
p_ev=eval(p2);
f_ev=eval(f);
H=f_ev==p_ev;
[a0,a1,a2,a3,a4,a5]=solve(H,[a0,a1,a2,a3,a4,a5]);
p2=poly2sym([a0,a1,a2,a3,a4,a5]);
%error
error_5=abs(f-p2);

%10-monomial
syms x a0 a1 a2 a3 a4 a5 a6 a7 a8 a9 a10
f=exp(1./x);
p3=poly2sym([a0 a1 a2 a3 a4 a5 a6 a7 a8 a9 a10]);
x=interpolation_n(-1,1,11);
p_ev=eval(p3);
f_ev=eval(f);
H=f_ev==p_ev;
[a0,a1,a2,a3,a4,a5,a6, a7, a8, a9, a10]=solve(H,[a0,a1,a2,a3,a4,a5,a6 a7 a8 a9 a10]);
p3=poly2sym([a0,a1,a2,a3,a4,a5,a6 a7 a8 a9 a10]);
%error
error_20=abs(f-p3);

%PLOTS
figure(1)
fplot(f,[-1 1],'b');
title('Approximation of $f=\exp{\frac{1}{x}}$ for $x\in [-1,1]$','interpreter','latex')
hold on
figure(1)
fplot(p,[-1 1],'r');
hold on
figure(1)
fplot(p2,[-1 1],'y');
hold on
figure(1)
fplot(p3,[-1 1],'m');
ylim([-5 5]);;
legend('$f=\exp{\frac{1}{x}}$', '3-monomial', '5-monomial', '10-monomial')
set(legend,'Interpreter','latex')
legend('show')
print -dpdf ex3a11.eps
hold off

figure(2)
fplot(error_3,[-1 1]);
title('Approx. Errors of $f=\exp{\frac{1}{x}}$ for $x\in [-1,1]$','interpreter','latex')
hold on
figure(2)
fplot(error_5,[-1 1]);
hold on
figure(2)
fplot(error_20,[-1 1]);
ylim([-5 5]);
legend('3-monomial', '5-monomial', '10-monomial')
set(legend,'Interpreter','latex')
legend('show')
print -dpdf ex3a12.eps
hold off


%% runge function
clear all;
clc
syms x a0 a1 a2 a3 
f=1/(1+25.*x^2);
%3-monomial
p=poly2sym([a0 a1 a2 a3]);
x=interpolation_n(-1,1,4);
p_ev=eval(p);
f_ev=eval(f);
H=f_ev==p_ev;
[a0,a1,a2,a3]=solve(H,[a0,a1,a2,a3]);
p=poly2sym([a0 a1 a2 a3]);
%error
error_3=abs(f-p);

%5-monomial
syms x a0 a1 a2 a3 a4 a5
f=1/(1+25*x^2);
p2=poly2sym([a0 a1 a2 a3 a4 a5]);
x=interpolation_n(-1,1,6);
p_ev=eval(p2);
f_ev=eval(f);
H=f_ev==p_ev;
[a0,a1,a2,a3,a4,a5]=solve(H,[a0,a1,a2,a3,a4,a5]);
p2=poly2sym([a0,a1,a2,a3,a4,a5]);

%error
error_5=abs(f-p2);


%10-monomial
syms x a0 a1 a2 a3 a4 a5 a6 a7 a8 a9 a10
f=1/(1+25*x^2);
p3=poly2sym([a0 a1 a2 a3 a4 a5 a6 a7 a8 a9 a10]);
x=interpolation_n(-1,1,11);
p_ev=eval(p3);
f_ev=eval(f);
H=f_ev==p_ev;
[a0,a1,a2,a3,a4,a5,a6, a7, a8, a9, a10]=solve(H,[a0,a1,a2,a3,a4,a5,a6 a7 a8 a9 a10]);
p3=poly2sym([a0,a1,a2,a3,a4,a5,a6 a7 a8 a9 a10]);
%error
error_20=abs(f-p3);

%PLOTS
figure(1)
fplot(f,[-1 1],'b');
title('Approximation of $f=\frac{1}{1+25x^2}$ for $x\in [-1,1]$','interpreter','latex')
hold on
figure(1)
fplot(p,[-1 1],'r');
hold on
figure(1)
fplot(p2,[-1 1],'y');
hold on
figure(1)
fplot(p3,[-1 1],'m');
legend('$f=\frac{1}{1+25x^2}$', '3-monomial', '5-monomial', '10-monomial')
set(legend,'Interpreter','latex')
legend('show')
print -dpdf ex3a21.eps
hold off

figure(2)
fplot(error_3,[-1 1]);
title('Approx. Errors of $f=\frac{1}{1+25x^2}$ for $x\in [-1,1]$','interpreter','latex')
hold on
figure(2)
fplot(error_5,[-1 1]);
hold on
figure(2)
fplot(error_20,[-1 1]);
legend('3-monomial', '5-monomial', '10-monomial')
set(legend,'Interpreter','latex')
legend('show')
print -dpdf ex3a22.eps
hold off

%% ramp function
clear all;
clc
syms x a0 a1 a2 a3 
f=(x+abs(x))./2;
%3-monomial
p=poly2sym([a0 a1 a2 a3]);
x=interpolation_n(-1,1,4);
p_ev=eval(p);
f_ev=eval(f);
H=f_ev==p_ev;
[a0,a1,a2,a3]=solve(H,[a0,a1,a2,a3]);
p=poly2sym([a0 a1 a2 a3]);
%error
error_3=abs(f-p);

%5-monomial
syms x a0 a1 a2 a3 a4 a5
f=(x+abs(x))./2;
p2=poly2sym([a0 a1 a2 a3 a4 a5]);
x=interpolation_n(-1,1,6);
p_ev=eval(p2);
f_ev=eval(f);
H=f_ev==p_ev;
[a0,a1,a2,a3,a4,a5]=solve(H,[a0,a1,a2,a3,a4,a5]);
p2=poly2sym([a0,a1,a2,a3,a4,a5]);
%error
error_5=abs(f-p2);

%10-monomial
syms x a0 a1 a2 a3 a4 a5 a6 a7 a8 a9 a10
f=(x+abs(x))./2;
p3=poly2sym([a0 a1 a2 a3 a4 a5 a6 a7 a8 a9 a10]);
x=interpolation_n(-1,1,11);
p_ev=eval(p3);
f_ev=eval(f);
H=f_ev==p_ev;
[a0,a1,a2,a3,a4,a5,a6, a7, a8, a9, a10]=solve(H,[a0,a1,a2,a3,a4,a5,a6 a7 a8 a9 a10]);
p3=poly2sym([a0,a1,a2,a3,a4,a5,a6 a7 a8 a9 a10]);
%error
error_20=abs(f-p3);

%PLOTS
figure(1)
fplot(f,[-1 1],'b');
title('Approximation of $f=\frac{x+|x|}{2}$ for $x\in [-1,1]$','interpreter','latex')
hold on
figure(1)
fplot(p,[-1 1],'r');
hold on
figure(1)
fplot(p2,[-1 1],'y');
hold on
figure(1)
fplot(p3,[-1 1],'m');
legend('$f=\frac{x+|x|}{2}$', '3-monomial', '5-monomial', '10-monomial')
set(legend,'Interpreter','latex')
legend('show')
print -dpdf ex3a31.eps
hold off

figure(2)
fplot(error_3,[-1 1]);
title('Approx. Errors of $f=\frac{x+|x|}{2}$ for $x\in [-1,1]$','interpreter','latex')
hold on
figure(2)
fplot(error_5,[-1 1]);
hold on
figure(2)
fplot(error_20,[-1 1]);
legend('3-monomial', '5-monomial', '10-monomial')
set(legend,'Interpreter','latex')
legend('show')
print -dpdf ex3a32.eps
hold off

%% Question 3 - part b
%% Chebyshev interpolation nodes and polynomial approximation
clear all;
clc

%% f=exp(1/x)
syms x a0 a1 a2 a3 
f=exp(1./x);

%3-monomial
p=poly2sym([a0 a1 a2 a3]);
x=interpolation_cheb(-1,1,4);
p_ev=eval(p);
f_ev=eval(f);
H=f_ev==p_ev;
[a0,a1,a2,a3]=solve(H,[a0,a1,a2,a3]);
p=poly2sym([a0 a1 a2 a3]);
%error
error_3=abs(f-p);

%5-monomial
syms x a0 a1 a2 a3 a4 a5
f=exp(1./x);
p2=poly2sym([a0 a1 a2 a3 a4 a5]);
x=interpolation_cheb(-1,1,6);
p_ev=eval(p2);
f_ev=eval(f);
H=f_ev==p_ev;
[a0,a1,a2,a3,a4,a5]=solve(H,[a0,a1,a2,a3,a4,a5]);
p2=poly2sym([a0,a1,a2,a3,a4,a5]);
%error
error_5=abs(f-p2);

%10-mononial
syms x a0 a1 a2 a3 a4 a5 a6 a7 a8 a9 a10
f=exp(1./x);
p3=poly2sym([a0 a1 a2 a3 a4 a5 a6 a7 a8 a9 a10]);
x=interpolation_cheb(-1,1,11);
p_ev=eval(p3);
f_ev=eval(f);
H=f_ev==p_ev;
[a0,a1,a2,a3,a4,a5,a6 a7 a8 a9 a10]=solve(H,[a0,a1,a2,a3,a4,a5,a6 a7 a8 a9 a10]);
p3=poly2sym([a0,a1,a2,a3,a4,a5,a6 a7 a8 a9 a10]);

%errors
error_20=abs(f-p3);

%PLOTS
figure(1)
fplot(f,[-1 1],'b');
title('Approximation of $f=\exp{\frac{1}{x}}$ for $x\in [-1,1]$','interpreter','latex')
hold on
figure(1)
fplot(p,[-1 1],'r');
hold on
figure(1)
fplot(p2,[-1 1],'y');
hold on
figure(1)
fplot(p3,[-1 1],'m');
ylim([-2 2]);
legend('$f=\exp{\frac{1}{x}}$', '3-monomial', '5-monomial', '10-monomial')
set(legend,'Interpreter','latex')
legend('show')
print -dpdf ex3b11.eps
hold off

figure(2)
fplot(error_3,[-1 1]);
title('Approx. Errors of $f=\exp{\frac{1}{x}}$ for $x\in [-1,1]$','interpreter','latex')
hold on
figure(2)
fplot(error_5,[-1 1]);
hold on
figure(2)
fplot(error_20,[-1 1]);
ylim([-2 2]);
legend('3-monomial', '5-monomial', '10-monomial')
set(legend,'Interpreter','latex')
legend('show')
print -dpdf ex3b12.eps
hold off

%% runge function
clear all;
clc
syms x a0 a1 a2 a3 
f=1/(1+25.*x^2);
%3-monomial
p=poly2sym([a0 a1 a2 a3]);
x=interpolation_cheb(-1,1,4);
p_ev=eval(p);
f_ev=eval(f);
H=f_ev==p_ev;
[a0,a1,a2,a3]=solve(H,[a0,a1,a2,a3]);
p=poly2sym([a0 a1 a2 a3]);
%error
error_3=abs(f-p);

%5-monomial
syms x a0 a1 a2 a3 a4 a5
f=1/(1+25*x^2);
p2=poly2sym([a0 a1 a2 a3 a4 a5]);
x=interpolation_cheb(-1,1,6);
p_ev=eval(p2);
f_ev=eval(f);
H=f_ev==p_ev;
[a0,a1,a2,a3,a4,a5]=solve(H,[a0,a1,a2,a3,a4,a5]);
p2=poly2sym([a0,a1,a2,a3,a4,a5]);

%error
error_5=abs(f-p2);


%10-monomial
syms x a0 a1 a2 a3 a4 a5 a6 a7 a8 a9 a10
f=1/(1+25*x^2);
p3=poly2sym([a0 a1 a2 a3 a4 a5 a6 a7 a8 a9 a10]);
x=interpolation_cheb(-1,1,11);
p_ev=eval(p3);
f_ev=eval(f);
H=f_ev==p_ev;
[a0,a1,a2,a3,a4,a5,a6, a7, a8, a9, a10]=solve(H,[a0,a1,a2,a3,a4,a5,a6 a7 a8 a9 a10]);
p3=poly2sym([a0,a1,a2,a3,a4,a5,a6 a7 a8 a9 a10]);
%error
error_20=abs(f-p3);

%PLOTS
figure(1)
fplot(f,[-1 1],'b');
title('Approximation of $f=\frac{1}{1+25x^2}$ for $x\in [-1,1]$','interpreter','latex')
hold on
figure(1)
fplot(p,[-1 1],'r');
hold on
figure(1)
fplot(p2,[-1 1],'y');
hold on
figure(1)
fplot(p3,[-1 1],'m');
legend('$f=\frac{1}{1+25x^2}$', '3-monomial', '5-monomial', '10-monomial')
set(legend,'Interpreter','latex')
legend('show')
print -dpdf ex3b21.eps
hold off

figure(2)
fplot(error_3,[-1 1]);
title('Approx. Errors of $f=\frac{1}{1+25x^2}$ for $x\in [-1,1]$','interpreter','latex')
hold on
figure(2)
fplot(error_5,[-1 1]);
hold on
figure(2)
fplot(error_20,[-1 1]);
legend('3-monomial', '5-monomial', '10-monomial')
set(legend,'Interpreter','latex')
legend('show')
print -dpdf ex3b22.eps
hold off

%% ramp function
clear all;
clc
syms x a0 a1 a2 a3 
f=(x+abs(x))./2;
%3-monomial
p=poly2sym([a0 a1 a2 a3]);
x=interpolation_cheb(-1,1,4);
p_ev=eval(p);
f_ev=eval(f);
H=f_ev==p_ev;
[a0,a1,a2,a3]=solve(H,[a0,a1,a2,a3]);
p=poly2sym([a0 a1 a2 a3]);
%error
error_3=abs(f-p);

%5-monomial
syms x a0 a1 a2 a3 a4 a5
f=(x+abs(x))./2;
p2=poly2sym([a0 a1 a2 a3 a4 a5]);
x=interpolation_cheb(-1,1,6);
p_ev=eval(p2);
f_ev=eval(f);
H=f_ev==p_ev;
[a0,a1,a2,a3,a4,a5]=solve(H,[a0,a1,a2,a3,a4,a5]);
p2=poly2sym([a0,a1,a2,a3,a4,a5]);
%error
error_5=abs(f-p2);

%10-monomial
syms x a0 a1 a2 a3 a4 a5 a6 a7 a8 a9 a10
f=(x+abs(x))./2;
p3=poly2sym([a0 a1 a2 a3 a4 a5 a6 a7 a8 a9 a10]);
x=interpolation_cheb(-1,1,11);
p_ev=eval(p3);
f_ev=eval(f);
H=f_ev==p_ev;
[a0,a1,a2,a3,a4,a5,a6, a7, a8, a9, a10]=solve(H,[a0,a1,a2,a3,a4,a5,a6 a7 a8 a9 a10]);
p3=poly2sym([a0,a1,a2,a3,a4,a5,a6 a7 a8 a9 a10]);
%error
error_20=abs(f-p3);

%PLOTS
figure(1)
fplot(f,[-1 1],'b');
title('Approximation of $f=\frac{x+|x|}{2}$ for $x\in [-1,1]$','interpreter','latex')
hold on
figure(1)
fplot(p,[-1 1],'r');
hold on
figure(1)
fplot(p2,[-1 1],'y');
hold on
figure(1)
fplot(p3,[-1 1],'m');
legend('$f=\frac{x+|x|}{2}$', '3-monomial', '5-monomial', '10-monomial')
set(legend,'Interpreter','latex')
legend('show')
print -dpdf ex3b31.eps
hold off

figure(2)
fplot(error_3,[-1 1]);
title('Approx. Errors of $f=\frac{x+|x|}{2}$ for $x\in [-1,1]$','interpreter','latex')
hold on
figure(2)
fplot(error_5,[-1 1]);
hold on
figure(2)
fplot(error_20,[-1 1]);
legend('3-monomial', '5-monomial', '10-monomial')
set(legend,'Interpreter','latex')
legend('show')
print -dpdf ex3b32.eps
hold off

%% Question 3 - part c
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

%% Question 4
%% Chebyshev interpolation and approximation of probability function
clear all;
clc

%% first combination rho1=1/0.2
syms x alpha rho1 rho2
func=exp(-alpha*x)/(rho1+rho2*exp(-alpha*x));

%first combination
alpha=1;
rho2=1/100;
rho1=1/0.2;
func1=eval(func);

[poly1_3]=cheb_poly(func1,0,10,3,4);
error1_3=abs(func1-poly1_3);

[poly1_5]=cheb_poly(func1,0,10,5,6);
error1_5=abs(func1-poly1_5);

[poly1_10]=cheb_poly(func1,0,10,10,11);
error1_10=abs(func1-poly1_10);

%% second combination: rho1=1/0.25
alpha=1;
rho2=1/100;
rho1=1/0.25;
func2=eval(func);

[poly2_3]=cheb_poly(func2,0,10,3,4);
error2_3=abs(func2-poly2_3);

[poly2_5]=cheb_poly(func2,0,10,5,6);
error2_5=abs(func2-poly2_5);

[poly2_10]=cheb_poly(func2,0,10,10,11);
error2_10=abs(func2-poly2_10);

%% plots
%order-3
figure(1)
fplot(func1,[0 10],'r:');
title('Chebyshev Approximation - order 3')
hold on
figure(1)
fplot(func2,[0 10],'b:');
hold on
figure(1)
fplot(poly1_3,[0 10],'r');
hold on
figure(1)
fplot(poly2_3,[0 10],'b');
hold on
figure(1)
fplot(error1_3,[0 10],'r--');
hold on 
figure(1)
fplot(error2_3,[0 10],'b--');
legend('f1:$\rho_1=\frac{1}{0.2}$ ','f2:$\rho_1=\frac{1}{0.25}$',  'approx. with $\rho_1=\frac{1}{0.2}$', 'approx. with $\rho_1=\frac{1}{0.25}$','error with $\rho_1=\frac{1}{0.2}$','error with $\rho_1=\frac{1}{0.25}$')
set(legend,'Interpreter','latex')
legend('show')
print -dpdf ex4a1.eps
hold off

%order-5
figure(2)
fplot(func1,[0 10],'r:');
title('Chebyshev Approximation - order 5')
hold on
figure(2)
fplot(func2,[0 10],'b:');
hold on
figure(2)
fplot(poly1_5,[0 10],'r');
hold on
figure(2)
fplot(poly2_5,[0 10],'b');
hold on
figure(2)
fplot(error1_5,[0 10],'r--');
hold on 
figure(2)
fplot(error2_5,[0 10],'b--');
legend('f1:$\rho_1=\frac{1}{0.2}$ ','f2:$\rho_1=\frac{1}{0.25}$',  'approx. with $\rho_1=\frac{1}{0.2}$', 'approx. with $\rho_1=\frac{1}{0.25}$','error with $\rho_1=\frac{1}{0.2}$','error with $\rho_1=\frac{1}{0.25}$')
set(legend,'Interpreter','latex')
legend('show')
print -dpdf ex4a2.eps
hold off

%order-10
figure(3)
fplot(func1,[0 10],'r:');
title('Chebyshev Approximation - order 10')
hold on
figure(3)
fplot(func2,[0 10],'b:');
hold on
figure(3)
fplot(poly1_10,[0 10],'r');
hold on
figure(3)
fplot(poly2_10,[0 10],'b');
hold on
figure(3)
fplot(error1_10,[0 10],'r--');
hold on 
figure(3)
fplot(error2_10,[0 10],'b--');
legend('f1:$\rho_1=\frac{1}{0.2}$ ','f2:$\rho_1=\frac{1}{0.25}$',  'approx. with $\rho_1=\frac{1}{0.2}$', 'approx. with $\rho_1=\frac{1}{0.25}$','error with $\rho_1=\frac{1}{0.2}$','error with $\rho_1=\frac{1}{0.25}$')
set(legend,'Interpreter','latex')
legend('show')
print -dpdf ex4a3.eps
hold off

%% Question 5 - first combination
%% Multivariate CES Function Approximation - sigma=0.25
clc
clear all;

syms k h alpha sigma
f=((1-alpha)*k^((sigma-1)/sigma)+alpha*h^((sigma-1)/sigma))^(sigma/(sigma-1));
alpha=0.5;
sigma=0.25;
f_e=eval(f);

%% A) Show that sigma is the elasticity of substitution (ES)
n=diff(f_e,k)/diff(f_e,h);
d=k/h;
es=diff(log(n))/diff(log(d));
ES=-(1/es);
if sigma == ES;
    disp(ES,'Sigma is the ES')
end

%% B) Labor share
%It's defined as the marginal productivity of labor
syms h k
f_h=diff(f,h);
h_share=(h*f_h)/f;
he_share=eval(h_share);

%% C) Chebyshec Approximation
syms k h
syms x y
k=x;
h=y;
f=eval(f_e);
[x]=interpolation_cheb(0,10,20);
[y]=interpolation_cheb(0,10,20);
wk=eval(f);
xk=x;
yk=y;
m=20;
%% order 3
n=3;
coef=zeros(1,(n+1));
for i=1:n+1
    c=wk(1:m).*cos((i-1).*acos(xk(1,(1:m)))).*cos((i-1).*acos(yk(1,(1:m))));
    d1=cos((i-1)*acos(xk(1,(1:m)))).* cos((i-1)*acos(xk(1,(1:m))));
    d2=cos((i-1)*acos(yk(1,(1:m)))).* cos((i-1)*acos(yk(1,(1:m))));
    coef(1,i)=sum(c)./(sum(d1)*sum(d2)); 
end
syms x y
psii1=cos((0:n).*acos(2*((x)/(10))-1));
psii2=cos((0:n).*acos(2*((y)/(10))-1));
p=coef(1,(1:n+1)).*psii1(1,(1:n+1)).*psii2(1,(1:n+1));
poly_3=sum(p);

%plotting
[x y]=meshgrid(0:1:10);

figure(1);
ax1=subplot(1,3,2);
z_p=eval(poly_3);
mesh(ax1,x,y,z_p);
title(ax1,'approximation');

z_f=eval(f);
ax2=subplot(1,3,1);
mesh(ax2,x,y,z_f);
title(ax2,'function');

z_error=abs(z_f-z_p);
ax3=subplot(1,3,3);
mesh(ax3,x,y,z_error);
title(ax3,'error');

suptitle('3rd-order Chebyshev polynomial')

print -dpdf ex5c3.eps

%% order 4
n=4;
coef=zeros(1,(n+1));
for i=1:n+1
    c=wk(1:m).*cos((i-1).*acos(xk(1,(1:m)))).*cos((i-1).*acos(yk(1,(1:m))));
    d1=cos((i-1)*acos(xk(1,(1:m)))).* cos((i-1)*acos(xk(1,(1:m))));
    d2=cos((i-1)*acos(yk(1,(1:m)))).* cos((i-1)*acos(yk(1,(1:m))));
    coef(1,i)=sum(c)./(sum(d1)*sum(d2)); 
end
syms x y
psii1=cos((0:n).*acos(2*((x)/(10))-1));
psii2=cos((0:n).*acos(2*((y)/(10))-1));
p=coef(1,(1:n+1)).*psii1(1,(1:n+1)).*psii2(1,(1:n+1));
poly_4=sum(p);
%plotting
[x y]=meshgrid(0:1:10);
figure(2)
z_p=eval(poly_4);
mesh(x,y,z_p);
hold on
z_f=eval(f);
mesh(x,y,z_f);
legend('Function approx order 4', 'Function')
hold off

%% order 5
n=5;
coef=zeros(1,(n+1));
for i=1:n+1
    c=wk(1:m).*cos((i-1).*acos(xk(1,(1:m)))).*cos((i-1).*acos(yk(1,(1:m))));
    d1=cos((i-1)*acos(xk(1,(1:m)))).* cos((i-1)*acos(xk(1,(1:m))));
    d2=cos((i-1)*acos(yk(1,(1:m)))).* cos((i-1)*acos(yk(1,(1:m))));
    coef(1,i)=sum(c)./(sum(d1)*sum(d2)); 
end
syms x y
psii1=cos((0:n).*acos(2*((x)/(10))-1));
psii2=cos((0:n).*acos(2*((y)/(10))-1));
p=coef(1,(1:n+1)).*psii1(1,(1:n+1)).*psii2(1,(1:n+1));
poly_5=sum(p);
%plotting
[x y]=meshgrid(0:1:10);
figure(3)
z_p=eval(poly_5);
mesh(x,y,z_p);
hold on
z_f=eval(f);
mesh(x,y,z_f);
legend('Function approx order 5', 'Function')
hold off

%% order 6
n=6;
coef=zeros(1,(n+1));
for i=1:n+1
    c=wk(1:m).*cos((i-1).*acos(xk(1,(1:m)))).*cos((i-1).*acos(yk(1,(1:m))));
    d1=cos((i-1)*acos(xk(1,(1:m)))).* cos((i-1)*acos(xk(1,(1:m))));
    d2=cos((i-1)*acos(yk(1,(1:m)))).* cos((i-1)*acos(yk(1,(1:m))));
    coef(1,i)=sum(c)./(sum(d1)*sum(d2)); 
end
syms x y
psii1=cos((0:n).*acos(2*((x)/(10))-1));
psii2=cos((0:n).*acos(2*((y)/(10))-1));
p=coef(1,(1:n+1)).*psii1(1,(1:n+1)).*psii2(1,(1:n+1));
poly_6=sum(p);
%plotting

[x y]=meshgrid(0:1:10);

figure(4);
ax1=subplot(1,3,2);
z_p=eval(poly_6);
mesh(ax1,x,y,z_p);
title(ax1,'approximation');

z_f=eval(f);
ax2=subplot(1,3,1);
mesh(ax2,x,y,z_f);
title(ax2,'function');

z_error=abs(z_f-z_p);
ax3=subplot(1,3,3);
mesh(ax3,x,y,z_error);
title(ax3,'error');

suptitle('6th-order Chebyshev polynomial')

print -dpdf ex5c6.eps
%% order 7
n=7;
coef=zeros(1,(n+1));
for i=1:n+1
    c=wk(1:m).*cos((i-1).*acos(xk(1,(1:m)))).*cos((i-1).*acos(yk(1,(1:m))));
    d1=cos((i-1)*acos(xk(1,(1:m)))).* cos((i-1)*acos(xk(1,(1:m))));
    d2=cos((i-1)*acos(yk(1,(1:m)))).* cos((i-1)*acos(yk(1,(1:m))));
    coef(1,i)=sum(c)./(sum(d1)*sum(d2)); 
end
syms x y
psii1=cos((0:n).*acos(2*((x)/(10))-1));
psii2=cos((0:n).*acos(2*((y)/(10))-1));
p=coef(1,(1:n+1)).*psii1(1,(1:n+1)).*psii2(1,(1:n+1));
poly_7=sum(p);
%plotting
[x y]=meshgrid(0:1:10);
figure(5)
z_p=eval(poly_7);
mesh(x,y,z_p);
hold on
z_f=eval(f);
mesh(x,y,z_f);
legend('Function approx order 7', 'Function')
hold off

%% order 8
n=8;
coef=zeros(1,(n+1));
for i=1:n+1
    c=wk(1:m).*cos((i-1).*acos(xk(1,(1:m)))).*cos((i-1).*acos(yk(1,(1:m))));
    d1=cos((i-1)*acos(xk(1,(1:m)))).* cos((i-1)*acos(xk(1,(1:m))));
    d2=cos((i-1)*acos(yk(1,(1:m)))).* cos((i-1)*acos(yk(1,(1:m))));
    coef(1,i)=sum(c)./(sum(d1)*sum(d2)); 
end
syms x y
psii1=cos((0:n).*acos(2*((x)/(10))-1));
psii2=cos((0:n).*acos(2*((y)/(10))-1));
p=coef(1,(1:n+1)).*psii1(1,(1:n+1)).*psii2(1,(1:n+1));
poly_8=sum(p);
%plotting
[x y]=meshgrid(0:1:10);
figure(6)
z_p=eval(poly_8);
mesh(x,y,z_p);
hold on
z_f=eval(f);
mesh(x,y,z_f);
legend('Function approx order 8', 'Function')
hold off

%% order 9
n=9;
coef=zeros(1,(n+1));
for i=1:n+1
    c=wk(1:m).*cos((i-1).*acos(xk(1,(1:m)))).*cos((i-1).*acos(yk(1,(1:m))));
    d1=cos((i-1)*acos(xk(1,(1:m)))).* cos((i-1)*acos(xk(1,(1:m))));
    d2=cos((i-1)*acos(yk(1,(1:m)))).* cos((i-1)*acos(yk(1,(1:m))));
    coef(1,i)=sum(c)./(sum(d1)*sum(d2)); 
end
syms x y
psii1=cos((0:n).*acos(2*((x)/(10))-1));
psii2=cos((0:n).*acos(2*((y)/(10))-1));
p=coef(1,(1:n+1)).*psii1(1,(1:n+1)).*psii2(1,(1:n+1));
poly_9=sum(p);
%plotting
[x y]=meshgrid(0:1:10);
figure(7)

ax1=subplot(1,3,2);
z_p=eval(poly_9);
mesh(ax1,x,y,z_p);
title(ax1,'approximation');

z_f=eval(f);
ax2=subplot(1,3,1);
mesh(ax2,x,y,z_f);
title(ax2,'function');

z_error=abs(z_f-z_p);
ax3=subplot(1,3,3);
mesh(ax3,x,y,z_error);
title(ax3,'error');

suptitle('9th-order Chebyshev polynomial')

print -dpdf ex5c9.eps

%% order 10
n=10;
coef=zeros(1,(n+1));
for i=1:n+1
    c=wk(1:m).*cos((i-1).*acos(xk(1,(1:m)))).*cos((i-1).*acos(yk(1,(1:m))));
    d1=cos((i-1)*acos(xk(1,(1:m)))).* cos((i-1)*acos(xk(1,(1:m))));
    d2=cos((i-1)*acos(yk(1,(1:m)))).* cos((i-1)*acos(yk(1,(1:m))));
    coef(1,i)=sum(c)./(sum(d1)*sum(d2)); 
end
syms x y
psii1=cos((0:n).*acos(2*((x)/(10))-1));
psii2=cos((0:n).*acos(2*((y)/(10))-1));
p=coef(1,(1:n+1)).*psii1(1,(1:n+1)).*psii2(1,(1:n+1));
poly_15=sum(p);
%plotting
[x y]=meshgrid(0:1:10);
figure(8)
z_p=eval(poly_10);
mesh(x,y,z_p);
hold on
z_f=eval(f);
mesh(x,y,z_f);
legend('Function approx order 10', 'Function')
hold off

%% order 11
n=11;
coef=zeros(1,(n+1));
for i=1:n+1
    c=wk(1:m).*cos((i-1).*acos(xk(1,(1:m)))).*cos((i-1).*acos(yk(1,(1:m))));
    d1=cos((i-1)*acos(xk(1,(1:m)))).* cos((i-1)*acos(xk(1,(1:m))));
    d2=cos((i-1)*acos(yk(1,(1:m)))).* cos((i-1)*acos(yk(1,(1:m))));
    coef(1,i)=sum(c)./(sum(d1)*sum(d2)); 
end
syms x y
psii1=cos((0:n).*acos(2*((x)/(10))-1));
psii2=cos((0:n).*acos(2*((y)/(10))-1));
p=coef(1,(1:n+1)).*psii1(1,(1:n+1)).*psii2(1,(1:n+1));
poly_11=sum(p);
%plotting
[x y]=meshgrid(0:1:10);
figure(9)
z_p=eval(poly_11);
mesh(x,y,z_p);
hold on
z_f=eval(f);
mesh(x,y,z_f);
legend('Function approx order 11', 'Function')
hold off

%% order 12
n=12;
coef=zeros(1,(n+1));
for i=1:n+1
    c=wk(1:m).*cos((i-1).*acos(xk(1,(1:m)))).*cos((i-1).*acos(yk(1,(1:m))));
    d1=cos((i-1)*acos(xk(1,(1:m)))).* cos((i-1)*acos(xk(1,(1:m))));
    d2=cos((i-1)*acos(yk(1,(1:m)))).* cos((i-1)*acos(yk(1,(1:m))));
    coef(1,i)=sum(c)./(sum(d1)*sum(d2)); 
end
syms x y
psii1=cos((0:n).*acos(2*((x)/(10))-1));
psii2=cos((0:n).*acos(2*((y)/(10))-1));
p=coef(1,(1:n+1)).*psii1(1,(1:n+1)).*psii2(1,(1:n+1));
poly_12=sum(p);
%plotting
[x y]=meshgrid(0:1:10);
figure(10);

ax1=subplot(1,3,2);
z_p=eval(poly_12);
mesh(ax1,x,y,z_p);
title(ax1,'approximation');

z_f=eval(f);
ax2=subplot(1,3,1);
mesh(ax2,x,y,z_f);
title(ax2,'function');

z_error=abs(z_f-z_p);
ax3=subplot(1,3,3);
mesh(ax3,x,y,z_error);
title(ax3,'error');

suptitle('12th-order Chebyshev polynomial');

print -dpdf ex5c12.eps
%% order 13
n=13;
coef=zeros(1,(n+1));
for i=1:n+1
    c=wk(1:m).*cos((i-1).*acos(xk(1,(1:m)))).*cos((i-1).*acos(yk(1,(1:m))));
    d1=cos((i-1)*acos(xk(1,(1:m)))).* cos((i-1)*acos(xk(1,(1:m))));
    d2=cos((i-1)*acos(yk(1,(1:m)))).* cos((i-1)*acos(yk(1,(1:m))));
    coef(1,i)=sum(c)./(sum(d1)*sum(d2)); 
end
syms x y
psii1=cos((0:n).*acos(2*((x)/(10))-1));
psii2=cos((0:n).*acos(2*((y)/(10))-1));
p=coef(1,(1:n+1)).*psii1(1,(1:n+1)).*psii2(1,(1:n+1));
poly_13=sum(p);
%plotting
[x y]=meshgrid(0:1:10);
figure(11)
z_p=eval(poly_13);
mesh(x,y,z_p);
hold on
z_f=eval(f);
mesh(x,y,z_f);
legend('Function approx order 13', 'Function')
hold off

%% order 14
n=14;
coef=zeros(1,(n+1));
for i=1:n+1
    c=wk(1:m).*cos((i-1).*acos(xk(1,(1:m)))).*cos((i-1).*acos(yk(1,(1:m))));
    d1=cos((i-1)*acos(xk(1,(1:m)))).* cos((i-1)*acos(xk(1,(1:m))));
    d2=cos((i-1)*acos(yk(1,(1:m)))).* cos((i-1)*acos(yk(1,(1:m))));
    coef(1,i)=sum(c)./(sum(d1)*sum(d2)); 
end
syms x y
psii1=cos((0:n).*acos(2*((x)/(10))-1));
psii2=cos((0:n).*acos(2*((y)/(10))-1));
p=coef(1,(1:n+1)).*psii1(1,(1:n+1)).*psii2(1,(1:n+1));
poly_14=sum(p);
%plotting
[x y]=meshgrid(0:1:10);
figure(12)
z_p=eval(poly_14);
mesh(x,y,z_p);
hold on
z_f=eval(f);
mesh(x,y,z_f);
legend('Function approx order 14', 'Function')
hold off

%% order 15
n=15;
coef=zeros(1,(n+1));
for i=1:n+1
    c=wk(1:m).*cos((i-1).*acos(xk(1,(1:m)))).*cos((i-1).*acos(yk(1,(1:m))));
    d1=cos((i-1)*acos(xk(1,(1:m)))).* cos((i-1)*acos(xk(1,(1:m))));
    d2=cos((i-1)*acos(yk(1,(1:m)))).* cos((i-1)*acos(yk(1,(1:m))));
    coef(1,i)=sum(c)./(sum(d1)*sum(d2)); 
end
syms x y
psii1=cos((0:n).*acos(2*((x)/(10))-1));
psii2=cos((0:n).*acos(2*((y)/(10))-1));
p=coef(1,(1:n+1)).*psii1(1,(1:n+1)).*psii2(1,(1:n+1));
poly_15=sum(p);
%plotting
[x y]=meshgrid(0:1:10);
figure(13)

ax1=subplot(1,3,2);
z_p=eval(poly_15);
mesh(ax1,x,y,z_p);
title(ax1,'approximation');

z_f=eval(f);
ax2=subplot(1,3,1);
mesh(ax2,x,y,z_f);
title(ax2,'function');

z_error=abs(z_f-z_p);
ax3=subplot(1,3,3);
mesh(ax3,x,y,z_error);
title(ax3,'error');

suptitle('15th-order Chebyshev polynomial')

print -dpdf ex5c15.eps

%% D) Isoquants
%isoquants with exact function and approximation of order 3
%function contour
Z=subs(f_e,'h','H');
Z=subs(Z,'k','K');
[K,H] = meshgrid(0:10, 0:10);
Z=(0.5.*(K).^((-0.75)/0.25)+(0.5).*(H).^(-0.75/0.25)).^(0.25/-0.75);
[C,h] = contour(K,H,Z);
print -dpdf ex5d1.eps
%approx
poly_3s=simplify(poly_3);
poly_3s=subs(poly_3s,'x','K');
poly_3s=subs(poly_3s,'y','H');
Z3=simplify(poly_3s);
Z3=(H/5 - 1)*(K/450 - 1/90) + (1540388770512533*((2*H^2)/25 - (4*H)/5 + 1)*((2*K^2)/25 - (4*K)/5 + 1))/36893488147419103232 + (19343440530901*((4*H^3)/125 - (12*H^2)/25 + (9*H)/5 - 1)*((4*K^3)/125 - (12*K^2)/25 + (9*K)/5 - 1))/147573952589676412928 + 1/4
[K,H] = meshgrid(0:10, 0:10);
[C,h] = contour(K,H,Z3);
print -dpdf ex5d2.eps

%approx vs function contour
[K,H] = meshgrid(0:10, 0:10);
[C,h] = contour(K,H,Z3);
hold on
[C,h] = contour(K,H,Z);
print -dpdf ex5d3.eps

%error contour
e=simplify(abs(f_e-poly_3));
e=subs(e,'x','K');
e=subs(e,'y','H');
e=subs(e,'k','K');
e=subs(e,'h','H');
e=abs((H./5 - 1).*(K./450 - 1/90) + (1540388770512533.*((2.*H.^2)./25 - (4.*H)./5 + 1).*((2.*K.^2)./25 - (4.*K)./5 + 1))./36893488147419103232 - 2^(1/3)./((H.^3 + K.^3)./(H.^3.*K.^3)).^(1/3) + (19343440530901.*((4.*H.^3)./125 - (12.*H.^2)./25 + (9.*H)./5 - 1).*((4.*K.^3)./125 - (12.*K.^2)./25 + (9.*K)./5 - 1))./147573952589676412928 + 1/4)
[K,H] = meshgrid(0:10, 0:10);
[C,h] = contour(K,H,e);
hold on
[C,h] = contour(K,H,Z3);
print -dpdf ex5d4.eps
hold off

%% Question 5 - second combination
%% Multivariate CES Function Approximation - sigma=5.0
clc
clear all;

syms k h alpha sigma
f=((1-alpha)*k^((sigma-1)/sigma)+alpha*h^((sigma-1)/sigma))^(sigma/(sigma-1));
alpha=0.5;
sigma=5;
f_e=eval(f);

%% A) Show that sigma is the elasticity of substitution (ES)
n=diff(f_e,k)/diff(f_e,h);
d=k/h;
es=diff(log(n))/diff(log(d));
ES=-(1/es);
if sigma == ES;
    disp(ES,'Sigma is the ES')
end

%% B) Labor share
%It's defined as the marginal productivity of labor
syms h k
f_h=diff(f,h);
h_share=(h*f_h)/f;
he_share=eval(h_share);

%% C) Chebyshec Approximation
syms k h
syms x y
k=x;
h=y;
f=eval(f_e);
[x]=interpolation_cheb(0,10,20);
[y]=interpolation_cheb(0,10,20);
wk=eval(f);
xk=x;
yk=y;
m=20;
%% order 3
n=3;
coef=zeros(1,(n+1));
for i=1:n+1
    c=wk(1:m).*cos((i-1).*acos(xk(1,(1:m)))).*cos((i-1).*acos(yk(1,(1:m))));
    d1=cos((i-1)*acos(xk(1,(1:m)))).* cos((i-1)*acos(xk(1,(1:m))));
    d2=cos((i-1)*acos(yk(1,(1:m)))).* cos((i-1)*acos(yk(1,(1:m))));
    coef(1,i)=sum(c)./(sum(d1)*sum(d2)); 
end
syms x y
psii1=cos((0:n).*acos(2*((x)/(10))-1));
psii2=cos((0:n).*acos(2*((y)/(10))-1));
p=coef(1,(1:n+1)).*psii1(1,(1:n+1)).*psii2(1,(1:n+1));
poly_3=sum(p);

%plotting
[x y]=meshgrid(0:1:10);

figure(1);
ax1=subplot(1,3,2);
z_p=eval(poly_3);
mesh(ax1,x,y,z_p);
title(ax1,'approximation');

z_f=eval(f);
ax2=subplot(1,3,1);
mesh(ax2,x,y,z_f);
title(ax2,'function');

z_error=abs(z_f-z_p);
ax3=subplot(1,3,3);
mesh(ax3,x,y,z_error);
title(ax3,'error');

suptitle('3rd-order Chebyshev polynomial')

print -dpdf ex5c3_2.eps

%% order 4
n=4;
coef=zeros(1,(n+1));
for i=1:n+1
    c=wk(1:m).*cos((i-1).*acos(xk(1,(1:m)))).*cos((i-1).*acos(yk(1,(1:m))));
    d1=cos((i-1)*acos(xk(1,(1:m)))).* cos((i-1)*acos(xk(1,(1:m))));
    d2=cos((i-1)*acos(yk(1,(1:m)))).* cos((i-1)*acos(yk(1,(1:m))));
    coef(1,i)=sum(c)./(sum(d1)*sum(d2)); 
end
syms x y
psii1=cos((0:n).*acos(2*((x)/(10))-1));
psii2=cos((0:n).*acos(2*((y)/(10))-1));
p=coef(1,(1:n+1)).*psii1(1,(1:n+1)).*psii2(1,(1:n+1));
poly_4=sum(p);

%% order 5
n=5;
coef=zeros(1,(n+1));
for i=1:n+1
    c=wk(1:m).*cos((i-1).*acos(xk(1,(1:m)))).*cos((i-1).*acos(yk(1,(1:m))));
    d1=cos((i-1)*acos(xk(1,(1:m)))).* cos((i-1)*acos(xk(1,(1:m))));
    d2=cos((i-1)*acos(yk(1,(1:m)))).* cos((i-1)*acos(yk(1,(1:m))));
    coef(1,i)=sum(c)./(sum(d1)*sum(d2)); 
end
syms x y
psii1=cos((0:n).*acos(2*((x)/(10))-1));
psii2=cos((0:n).*acos(2*((y)/(10))-1));
p=coef(1,(1:n+1)).*psii1(1,(1:n+1)).*psii2(1,(1:n+1));
poly_5=sum(p);

%% order 6
n=6;
coef=zeros(1,(n+1));
for i=1:n+1
    c=wk(1:m).*cos((i-1).*acos(xk(1,(1:m)))).*cos((i-1).*acos(yk(1,(1:m))));
    d1=cos((i-1)*acos(xk(1,(1:m)))).* cos((i-1)*acos(xk(1,(1:m))));
    d2=cos((i-1)*acos(yk(1,(1:m)))).* cos((i-1)*acos(yk(1,(1:m))));
    coef(1,i)=sum(c)./(sum(d1)*sum(d2)); 
end
syms x y
psii1=cos((0:n).*acos(2*((x)/(10))-1));
psii2=cos((0:n).*acos(2*((y)/(10))-1));
p=coef(1,(1:n+1)).*psii1(1,(1:n+1)).*psii2(1,(1:n+1));
poly_6=sum(p);
%plotting

[x y]=meshgrid(0:1:10);

figure(4);
ax1=subplot(1,3,2);
z_p=eval(poly_6);
mesh(ax1,x,y,z_p);
title(ax1,'approximation');

z_f=eval(f);
ax2=subplot(1,3,1);
mesh(ax2,x,y,z_f);
title(ax2,'function');

z_error=abs(z_f-z_p);
ax3=subplot(1,3,3);
mesh(ax3,x,y,z_error);
title(ax3,'error');

suptitle('6th-order Chebyshev polynomial')

print -dpdf ex5c6_2.eps
%% order 7
n=7;
coef=zeros(1,(n+1));
for i=1:n+1
    c=wk(1:m).*cos((i-1).*acos(xk(1,(1:m)))).*cos((i-1).*acos(yk(1,(1:m))));
    d1=cos((i-1)*acos(xk(1,(1:m)))).* cos((i-1)*acos(xk(1,(1:m))));
    d2=cos((i-1)*acos(yk(1,(1:m)))).* cos((i-1)*acos(yk(1,(1:m))));
    coef(1,i)=sum(c)./(sum(d1)*sum(d2)); 
end
syms x y
psii1=cos((0:n).*acos(2*((x)/(10))-1));
psii2=cos((0:n).*acos(2*((y)/(10))-1));
p=coef(1,(1:n+1)).*psii1(1,(1:n+1)).*psii2(1,(1:n+1));
poly_7=sum(p);

%% order 8
n=8;
coef=zeros(1,(n+1));
for i=1:n+1
    c=wk(1:m).*cos((i-1).*acos(xk(1,(1:m)))).*cos((i-1).*acos(yk(1,(1:m))));
    d1=cos((i-1)*acos(xk(1,(1:m)))).* cos((i-1)*acos(xk(1,(1:m))));
    d2=cos((i-1)*acos(yk(1,(1:m)))).* cos((i-1)*acos(yk(1,(1:m))));
    coef(1,i)=sum(c)./(sum(d1)*sum(d2)); 
end
syms x y
psii1=cos((0:n).*acos(2*((x)/(10))-1));
psii2=cos((0:n).*acos(2*((y)/(10))-1));
p=coef(1,(1:n+1)).*psii1(1,(1:n+1)).*psii2(1,(1:n+1));
poly_8=sum(p);
%% order 9
n=9;
coef=zeros(1,(n+1));
for i=1:n+1
    c=wk(1:m).*cos((i-1).*acos(xk(1,(1:m)))).*cos((i-1).*acos(yk(1,(1:m))));
    d1=cos((i-1)*acos(xk(1,(1:m)))).* cos((i-1)*acos(xk(1,(1:m))));
    d2=cos((i-1)*acos(yk(1,(1:m)))).* cos((i-1)*acos(yk(1,(1:m))));
    coef(1,i)=sum(c)./(sum(d1)*sum(d2)); 
end
syms x y
psii1=cos((0:n).*acos(2*((x)/(10))-1));
psii2=cos((0:n).*acos(2*((y)/(10))-1));
p=coef(1,(1:n+1)).*psii1(1,(1:n+1)).*psii2(1,(1:n+1));
poly_9=sum(p);
%plotting
[x y]=meshgrid(0:1:10);
figure(7)

ax1=subplot(1,3,2);
z_p=eval(poly_9);
mesh(ax1,x,y,z_p);
title(ax1,'approximation');

z_f=eval(f);
ax2=subplot(1,3,1);
mesh(ax2,x,y,z_f);
title(ax2,'function');

z_error=abs(z_f-z_p);
ax3=subplot(1,3,3);
mesh(ax3,x,y,z_error);
title(ax3,'error');

suptitle('9th-order Chebyshev polynomial')

print -dpdf ex5c9_2.eps

%% order 10
n=10;
coef=zeros(1,(n+1));
for i=1:n+1
    c=wk(1:m).*cos((i-1).*acos(xk(1,(1:m)))).*cos((i-1).*acos(yk(1,(1:m))));
    d1=cos((i-1)*acos(xk(1,(1:m)))).* cos((i-1)*acos(xk(1,(1:m))));
    d2=cos((i-1)*acos(yk(1,(1:m)))).* cos((i-1)*acos(yk(1,(1:m))));
    coef(1,i)=sum(c)./(sum(d1)*sum(d2)); 
end
syms x y
psii1=cos((0:n).*acos(2*((x)/(10))-1));
psii2=cos((0:n).*acos(2*((y)/(10))-1));
p=coef(1,(1:n+1)).*psii1(1,(1:n+1)).*psii2(1,(1:n+1));
poly_15=sum(p);

%% order 11
n=11;
coef=zeros(1,(n+1));
for i=1:n+1
    c=wk(1:m).*cos((i-1).*acos(xk(1,(1:m)))).*cos((i-1).*acos(yk(1,(1:m))));
    d1=cos((i-1)*acos(xk(1,(1:m)))).* cos((i-1)*acos(xk(1,(1:m))));
    d2=cos((i-1)*acos(yk(1,(1:m)))).* cos((i-1)*acos(yk(1,(1:m))));
    coef(1,i)=sum(c)./(sum(d1)*sum(d2)); 
end
syms x y
psii1=cos((0:n).*acos(2*((x)/(10))-1));
psii2=cos((0:n).*acos(2*((y)/(10))-1));
p=coef(1,(1:n+1)).*psii1(1,(1:n+1)).*psii2(1,(1:n+1));
poly_11=sum(p);

%% order 12
n=12;
coef=zeros(1,(n+1));
for i=1:n+1
    c=wk(1:m).*cos((i-1).*acos(xk(1,(1:m)))).*cos((i-1).*acos(yk(1,(1:m))));
    d1=cos((i-1)*acos(xk(1,(1:m)))).* cos((i-1)*acos(xk(1,(1:m))));
    d2=cos((i-1)*acos(yk(1,(1:m)))).* cos((i-1)*acos(yk(1,(1:m))));
    coef(1,i)=sum(c)./(sum(d1)*sum(d2)); 
end
syms x y
psii1=cos((0:n).*acos(2*((x)/(10))-1));
psii2=cos((0:n).*acos(2*((y)/(10))-1));
p=coef(1,(1:n+1)).*psii1(1,(1:n+1)).*psii2(1,(1:n+1));
poly_12=sum(p);
%plotting
[x y]=meshgrid(0:1:10);
figure(10);

ax1=subplot(1,3,2);
z_p=eval(poly_12);
mesh(ax1,x,y,z_p);
title(ax1,'approximation');

z_f=eval(f);
ax2=subplot(1,3,1);
mesh(ax2,x,y,z_f);
title(ax2,'function');

z_error=abs(z_f-z_p);
ax3=subplot(1,3,3);
mesh(ax3,x,y,z_error);
title(ax3,'error');

suptitle('12th-order Chebyshev polynomial');

print -dpdf ex5c12_2.eps
%% order 13
n=13;
coef=zeros(1,(n+1));
for i=1:n+1
    c=wk(1:m).*cos((i-1).*acos(xk(1,(1:m)))).*cos((i-1).*acos(yk(1,(1:m))));
    d1=cos((i-1)*acos(xk(1,(1:m)))).* cos((i-1)*acos(xk(1,(1:m))));
    d2=cos((i-1)*acos(yk(1,(1:m)))).* cos((i-1)*acos(yk(1,(1:m))));
    coef(1,i)=sum(c)./(sum(d1)*sum(d2)); 
end
syms x y
psii1=cos((0:n).*acos(2*((x)/(10))-1));
psii2=cos((0:n).*acos(2*((y)/(10))-1));
p=coef(1,(1:n+1)).*psii1(1,(1:n+1)).*psii2(1,(1:n+1));
poly_13=sum(p);

%% order 14
n=14;
coef=zeros(1,(n+1));
for i=1:n+1
    c=wk(1:m).*cos((i-1).*acos(xk(1,(1:m)))).*cos((i-1).*acos(yk(1,(1:m))));
    d1=cos((i-1)*acos(xk(1,(1:m)))).* cos((i-1)*acos(xk(1,(1:m))));
    d2=cos((i-1)*acos(yk(1,(1:m)))).* cos((i-1)*acos(yk(1,(1:m))));
    coef(1,i)=sum(c)./(sum(d1)*sum(d2)); 
end
syms x y
psii1=cos((0:n).*acos(2*((x)/(10))-1));
psii2=cos((0:n).*acos(2*((y)/(10))-1));
p=coef(1,(1:n+1)).*psii1(1,(1:n+1)).*psii2(1,(1:n+1));
poly_14=sum(p);

%% order 15
n=15;
coef=zeros(1,(n+1));
for i=1:n+1
    c=wk(1:m).*cos((i-1).*acos(xk(1,(1:m)))).*cos((i-1).*acos(yk(1,(1:m))));
    d1=cos((i-1)*acos(xk(1,(1:m)))).* cos((i-1)*acos(xk(1,(1:m))));
    d2=cos((i-1)*acos(yk(1,(1:m)))).* cos((i-1)*acos(yk(1,(1:m))));
    coef(1,i)=sum(c)./(sum(d1)*sum(d2)); 
end
syms x y
psii1=cos((0:n).*acos(2*((x)/(10))-1));
psii2=cos((0:n).*acos(2*((y)/(10))-1));
p=coef(1,(1:n+1)).*psii1(1,(1:n+1)).*psii2(1,(1:n+1));
poly_15=sum(p);
%plotting
[x y]=meshgrid(0:1:10);

figure(13)
z_f=eval(f);
ax1=subplot(1,3,1);
mesh(ax1,x,y,z_f);
title(ax1,'function');

ax2=subplot(1,3,2);
z_p=eval(poly_15);
mesh(ax2,x,y,z_p);
title(ax2,'approximation');

z_error=abs(z_f-z_p);
ax3=subplot(1,3,3);
mesh(ax3,x,y,z_error);
title(ax3,'error');

suptitle('15th-order Chebyshev polynomial')

print -dpdf ex5c15_2.eps

%% D) Isoquants
%isoquants with exact function and approximation of order 3
%function contour
Z=subs(f_e,'h','H');
Z=subs(Z,'k','K');
Z =(H.^(4/5)./2 + K.^(4/5)./2).^(5/4);
[K,H] = meshgrid(0:10, 0:10);
[C,h] = contour(K,H,Z);
print -dpdf ex5d1_2.eps
%approx
poly_3s=simplify(poly_3);
poly_3s=subs(poly_3s,'x','K');
poly_3s=subs(poly_3s,'y','H');
Z3=simplify(poly_3s);
Z3 =(H/5 - 1)*(K/450 - 1/90) + (3080777541025067*((2*H^2)/25 - (4*H)/5 + 1)*((2*K^2)/25 - (4*K)/5 + 1))/73786976294838206464 + (4951920775910657*((4*H^3)/125 - (12*H^2)/25 + (9*H)/5 - 1)*((4*K^3)/125 - (12*K^2)/25 + (9*K)/5 - 1))/37778931862957161709568 + 1/4;
[K,H] = meshgrid(0:10, 0:10);
[C,h] = contour(K,H,Z3);
print -dpdf ex5d2_2.eps

%approx vs function contour
[K,H] = meshgrid(0:10, 0:10);
[C,h] = contour(K,H,Z3);
hold on
[C,h] = contour(K,H,Z);
print -dpdf ex5d3_2.eps

%error contour
e=simplify(abs(f_e-poly_3));
e=subs(e,'x','K');
e=subs(e,'y','H');
e=subs(e,'k','K');
e=subs(e,'h','H');
e=abs((3080777541025067.*((2.*H.^2)./25 - (4.*H)./5 + 1).*((2.*K.^2)./25 - (4.*K)./5 + 1))./73786976294838206464 - K./450 - H./450 + (4951920775910657.*((4.*H.^3)./125 - (12.*H.^2)./25 + (9.*H)./5 - 1).*((4.*K.^3)./125 - (12.*K.^2)./25 + (9.*K)./5 - 1))./37778931862957161709568 - (2.^(3/4).*(H.^(4/5) + K.^(4/5)).^(5/4))./4 + (H.*K)./2250 + 47/180);

[K,H] = meshgrid(0:10, 0:10);
[C,h] = contour(K,H,e);
hold on
[C,h] = contour(K,H,Z3);
print -dpdf ex5d4_2.eps
hold off

%% Question 5 - third combination

%% Multivariate CES Function Approximation - sigma=1
clc
clear all;

syms k h alpha sigma
f=((1-alpha)*k^((sigma-1)/sigma)+alpha*h^((sigma-1)/sigma))^(sigma/(sigma-1));
alpha=0.5;
sigma=0.99;
f_e=eval(f);

%% A) Show that sigma is the elasticity of substitution (ES)
n=diff(f_e,k)/diff(f_e,h);
d=k/h;
es=diff(log(n))/diff(log(d));
ES=-(1/es);
if sigma == ES;
    disp(ES,'Sigma is the ES')
end

%% B) Labor share
%It's defined as the marginal productivity of labor
syms h k
f_h=diff(f,h);
h_share=(h*f_h)/f;
he_share=eval(h_share);

%% C) Chebyshec Approximation
syms k h
syms x y
k=x;
h=y;
f=eval(f_e);
%f=eval(f_e) not used since it's constant
[x]=interpolation_cheb(0,10,20);
[y]=interpolation_cheb(0,10,20);
wk=eval(f);
%wk=ones(1,20);
xk=x;
yk=y;
m=20;
%% order 3
n=3;
coef=zeros(1,(n+1));
for i=1:n+1
    c=wk(1:m).*cos((i-1).*acos(xk(1,(1:m)))).*cos((i-1).*acos(yk(1,(1:m))));
    d1=cos((i-1)*acos(xk(1,(1:m)))).* cos((i-1)*acos(xk(1,(1:m))));
    d2=cos((i-1)*acos(yk(1,(1:m)))).* cos((i-1)*acos(yk(1,(1:m))));
    coef(1,i)=sum(c)./(sum(d1)*sum(d2)); 
end
syms x y
psii1=cos((0:n).*acos(2*((x)/(10))-1));
psii2=cos((0:n).*acos(2*((y)/(10))-1));
p=coef(1,(1:n+1)).*psii1(1,(1:n+1)).*psii2(1,(1:n+1));
poly_3=sum(p);

%plotting
[x y]=meshgrid(0:1:10);

figure(1);
ax1=subplot(1,3,2);
z_p=eval(poly_3);
mesh(ax1,x,y,z_p);
title(ax1,'approximation');

z_f=ones(11,11);
ax2=subplot(1,3,1);
mesh(ax2,x,y,z_f);
title(ax2,'function');

z_error=abs(z_f-z_p);
ax3=subplot(1,3,3);
mesh(ax3,x,y,z_error);
title(ax3,'error');

suptitle('3rd-order Chebyshev polynomial')

print -dpdf ex5c3_3.eps

%% order 4
n=4;
coef=zeros(1,(n+1));
for i=1:n+1
    c=wk(1:m).*cos((i-1).*acos(xk(1,(1:m)))).*cos((i-1).*acos(yk(1,(1:m))));
    d1=cos((i-1)*acos(xk(1,(1:m)))).* cos((i-1)*acos(xk(1,(1:m))));
    d2=cos((i-1)*acos(yk(1,(1:m)))).* cos((i-1)*acos(yk(1,(1:m))));
    coef(1,i)=sum(c)./(sum(d1)*sum(d2)); 
end
syms x y
psii1=cos((0:n).*acos(2*((x)/(10))-1));
psii2=cos((0:n).*acos(2*((y)/(10))-1));
p=coef(1,(1:n+1)).*psii1(1,(1:n+1)).*psii2(1,(1:n+1));
poly_4=sum(p);

%% order 5
n=5;
coef=zeros(1,(n+1));
for i=1:n+1
    c=wk(1:m).*cos((i-1).*acos(xk(1,(1:m)))).*cos((i-1).*acos(yk(1,(1:m))));
    d1=cos((i-1)*acos(xk(1,(1:m)))).* cos((i-1)*acos(xk(1,(1:m))));
    d2=cos((i-1)*acos(yk(1,(1:m)))).* cos((i-1)*acos(yk(1,(1:m))));
    coef(1,i)=sum(c)./(sum(d1)*sum(d2)); 
end
syms x y
psii1=cos((0:n).*acos(2*((x)/(10))-1));
psii2=cos((0:n).*acos(2*((y)/(10))-1));
p=coef(1,(1:n+1)).*psii1(1,(1:n+1)).*psii2(1,(1:n+1));
poly_5=sum(p);


%% order 6
n=6;
coef=zeros(1,(n+1));
for i=1:n+1
    c=wk(1:m).*cos((i-1).*acos(xk(1,(1:m)))).*cos((i-1).*acos(yk(1,(1:m))));
    d1=cos((i-1)*acos(xk(1,(1:m)))).* cos((i-1)*acos(xk(1,(1:m))));
    d2=cos((i-1)*acos(yk(1,(1:m)))).* cos((i-1)*acos(yk(1,(1:m))));
    coef(1,i)=sum(c)./(sum(d1)*sum(d2)); 
end
syms x y
psii1=cos((0:n).*acos(2*((x)/(10))-1));
psii2=cos((0:n).*acos(2*((y)/(10))-1));
p=coef(1,(1:n+1)).*psii1(1,(1:n+1)).*psii2(1,(1:n+1));
poly_6=sum(p);
%plotting

[x y]=meshgrid(0:1:10);

figure(4);
ax1=subplot(1,3,2);
z_p=eval(poly_6);
mesh(ax1,x,y,z_p);
title(ax1,'approximation');

z_f=ones(11,11);
ax2=subplot(1,3,1);
mesh(ax2,x,y,z_f);
title(ax2,'function');

z_error=abs(z_f-z_p);
ax3=subplot(1,3,3);
mesh(ax3,x,y,z_error);
title(ax3,'error');

suptitle('6th-order Chebyshev polynomial')

print -dpdf ex5c6_3.eps
%% order 7
n=7;
coef=zeros(1,(n+1));
for i=1:n+1
    c=wk(1:m).*cos((i-1).*acos(xk(1,(1:m)))).*cos((i-1).*acos(yk(1,(1:m))));
    d1=cos((i-1)*acos(xk(1,(1:m)))).* cos((i-1)*acos(xk(1,(1:m))));
    d2=cos((i-1)*acos(yk(1,(1:m)))).* cos((i-1)*acos(yk(1,(1:m))));
    coef(1,i)=sum(c)./(sum(d1)*sum(d2)); 
end
syms x y
psii1=cos((0:n).*acos(2*((x)/(10))-1));
psii2=cos((0:n).*acos(2*((y)/(10))-1));
p=coef(1,(1:n+1)).*psii1(1,(1:n+1)).*psii2(1,(1:n+1));
poly_7=sum(p);

%% order 8
n=8;
coef=zeros(1,(n+1));
for i=1:n+1
    c=wk(1:m).*cos((i-1).*acos(xk(1,(1:m)))).*cos((i-1).*acos(yk(1,(1:m))));
    d1=cos((i-1)*acos(xk(1,(1:m)))).* cos((i-1)*acos(xk(1,(1:m))));
    d2=cos((i-1)*acos(yk(1,(1:m)))).* cos((i-1)*acos(yk(1,(1:m))));
    coef(1,i)=sum(c)./(sum(d1)*sum(d2)); 
end
syms x y
psii1=cos((0:n).*acos(2*((x)/(10))-1));
psii2=cos((0:n).*acos(2*((y)/(10))-1));
p=coef(1,(1:n+1)).*psii1(1,(1:n+1)).*psii2(1,(1:n+1));
poly_8=sum(p);


%% order 9
n=9;
coef=zeros(1,(n+1));
for i=1:n+1
    c=wk(1:m).*cos((i-1).*acos(xk(1,(1:m)))).*cos((i-1).*acos(yk(1,(1:m))));
    d1=cos((i-1)*acos(xk(1,(1:m)))).* cos((i-1)*acos(xk(1,(1:m))));
    d2=cos((i-1)*acos(yk(1,(1:m)))).* cos((i-1)*acos(yk(1,(1:m))));
    coef(1,i)=sum(c)./(sum(d1)*sum(d2)); 
end
syms x y
psii1=cos((0:n).*acos(2*((x)/(10))-1));
psii2=cos((0:n).*acos(2*((y)/(10))-1));
p=coef(1,(1:n+1)).*psii1(1,(1:n+1)).*psii2(1,(1:n+1));
poly_9=sum(p);
%plotting
[x y]=meshgrid(0:1:10);
figure(7)

ax1=subplot(1,3,2);
z_p=eval(poly_9);
mesh(ax1,x,y,z_p);
title(ax1,'approximation');

z_f=ones(11,11);
ax2=subplot(1,3,1);
mesh(ax2,x,y,z_f);
title(ax2,'function');

z_error=abs(z_f-z_p);
ax3=subplot(1,3,3);
mesh(ax3,x,y,z_error);
title(ax3,'error');

suptitle('9th-order Chebyshev polynomial')

print -dpdf ex5c9_3.eps

%% order 10
n=10;
coef=zeros(1,(n+1));
for i=1:n+1
    c=wk(1:m).*cos((i-1).*acos(xk(1,(1:m)))).*cos((i-1).*acos(yk(1,(1:m))));
    d1=cos((i-1)*acos(xk(1,(1:m)))).* cos((i-1)*acos(xk(1,(1:m))));
    d2=cos((i-1)*acos(yk(1,(1:m)))).* cos((i-1)*acos(yk(1,(1:m))));
    coef(1,i)=sum(c)./(sum(d1)*sum(d2)); 
end
syms x y
psii1=cos((0:n).*acos(2*((x)/(10))-1));
psii2=cos((0:n).*acos(2*((y)/(10))-1));
p=coef(1,(1:n+1)).*psii1(1,(1:n+1)).*psii2(1,(1:n+1));
poly_15=sum(p);

%% order 11
n=11;
coef=zeros(1,(n+1));
for i=1:n+1
    c=wk(1:m).*cos((i-1).*acos(xk(1,(1:m)))).*cos((i-1).*acos(yk(1,(1:m))));
    d1=cos((i-1)*acos(xk(1,(1:m)))).* cos((i-1)*acos(xk(1,(1:m))));
    d2=cos((i-1)*acos(yk(1,(1:m)))).* cos((i-1)*acos(yk(1,(1:m))));
    coef(1,i)=sum(c)./(sum(d1)*sum(d2)); 
end
syms x y
psii1=cos((0:n).*acos(2*((x)/(10))-1));
psii2=cos((0:n).*acos(2*((y)/(10))-1));
p=coef(1,(1:n+1)).*psii1(1,(1:n+1)).*psii2(1,(1:n+1));
poly_11=sum(p);


%% order 12
n=12;
coef=zeros(1,(n+1));
for i=1:n+1
    c=wk(1:m).*cos((i-1).*acos(xk(1,(1:m)))).*cos((i-1).*acos(yk(1,(1:m))));
    d1=cos((i-1)*acos(xk(1,(1:m)))).* cos((i-1)*acos(xk(1,(1:m))));
    d2=cos((i-1)*acos(yk(1,(1:m)))).* cos((i-1)*acos(yk(1,(1:m))));
    coef(1,i)=sum(c)./(sum(d1)*sum(d2)); 
end
syms x y
psii1=cos((0:n).*acos(2*((x)/(10))-1));
psii2=cos((0:n).*acos(2*((y)/(10))-1));
p=coef(1,(1:n+1)).*psii1(1,(1:n+1)).*psii2(1,(1:n+1));
poly_12=sum(p);
%plotting
[x y]=meshgrid(0:1:10);
figure(10);

ax1=subplot(1,3,2);
z_p=eval(poly_12);
mesh(ax1,x,y,z_p);
title(ax1,'approximation');

z_f=ones(11,11);
ax2=subplot(1,3,1);
mesh(ax2,x,y,z_f);
title(ax2,'function');

z_error=abs(z_f-z_p);
ax3=subplot(1,3,3);
mesh(ax3,x,y,z_error);
title(ax3,'error');

suptitle('12th-order Chebyshev polynomial');

print -dpdf ex5c12_3.eps
%% order 13
n=13;
coef=zeros(1,(n+1));
for i=1:n+1
    c=wk(1:m).*cos((i-1).*acos(xk(1,(1:m)))).*cos((i-1).*acos(yk(1,(1:m))));
    d1=cos((i-1)*acos(xk(1,(1:m)))).* cos((i-1)*acos(xk(1,(1:m))));
    d2=cos((i-1)*acos(yk(1,(1:m)))).* cos((i-1)*acos(yk(1,(1:m))));
    coef(1,i)=sum(c)./(sum(d1)*sum(d2)); 
end
syms x y
psii1=cos((0:n).*acos(2*((x)/(10))-1));
psii2=cos((0:n).*acos(2*((y)/(10))-1));
p=coef(1,(1:n+1)).*psii1(1,(1:n+1)).*psii2(1,(1:n+1));
poly_13=sum(p);
%plotting
[x y]=meshgrid(0:1:10);
figure(11)
z_p=eval(poly_13);
mesh(x,y,z_p);
hold on
z_f=eval(f);
mesh(x,y,z_f);
legend('Function approx order 13', 'Function')
hold off

%% order 14
n=14;
coef=zeros(1,(n+1));
for i=1:n+1
    c=wk(1:m).*cos((i-1).*acos(xk(1,(1:m)))).*cos((i-1).*acos(yk(1,(1:m))));
    d1=cos((i-1)*acos(xk(1,(1:m)))).* cos((i-1)*acos(xk(1,(1:m))));
    d2=cos((i-1)*acos(yk(1,(1:m)))).* cos((i-1)*acos(yk(1,(1:m))));
    coef(1,i)=sum(c)./(sum(d1)*sum(d2)); 
end
syms x y
psii1=cos((0:n).*acos(2*((x)/(10))-1));
psii2=cos((0:n).*acos(2*((y)/(10))-1));
p=coef(1,(1:n+1)).*psii1(1,(1:n+1)).*psii2(1,(1:n+1));
poly_14=sum(p);
%plotting
[x y]=meshgrid(0:1:10);
figure(12)
z_p=eval(poly_14);
mesh(x,y,z_p);
hold on
z_f=eval(f);
mesh(x,y,z_f);
legend('Function approx order 14', 'Function')
hold off

%% order 15
n=15;
coef=zeros(1,(n+1));
for i=1:n+1
    c=wk(1:m).*cos((i-1).*acos(xk(1,(1:m)))).*cos((i-1).*acos(yk(1,(1:m))));
    d1=cos((i-1)*acos(xk(1,(1:m)))).* cos((i-1)*acos(xk(1,(1:m))));
    d2=cos((i-1)*acos(yk(1,(1:m)))).* cos((i-1)*acos(yk(1,(1:m))));
    coef(1,i)=sum(c)./(sum(d1)*sum(d2)); 
end
syms x y
psii1=cos((0:n).*acos(2*((x)/(10))-1));
psii2=cos((0:n).*acos(2*((y)/(10))-1));
p=coef(1,(1:n+1)).*psii1(1,(1:n+1)).*psii2(1,(1:n+1));
poly_15=sum(p);
%plotting
[x y]=meshgrid(0:1:10);
figure(13)

ax1=subplot(1,3,2);
z_p=eval(poly_15);
mesh(ax1,x,y,z_p);
title(ax1,'approximation');

z_f=ones(11,11);
ax2=subplot(1,3,1);
mesh(ax2,x,y,z_f);
title(ax2,'function');

z_error=abs(z_f-z_p);
ax3=subplot(1,3,3);
mesh(ax3,x,y,z_error);
title(ax3,'error');

suptitle('15th-order Chebyshev polynomial')

print -dpdf ex5c15.eps


%% D) Isoquants
%isoquants with exact function and approximation of order 3
%function contour
[K,H] = meshgrid(0:10, 0:10);
Z=subs(f_e,'h','H');
Z=subs(Z,'k','K');
Z =1./(1./(2.*H.^(1./99)) + 1./(2.*K.^(1/99))).^99;
[K,H] = meshgrid(0:10, 0:10);
[C,h] = contour(K,H,Z);
print -dpdf ex5d1_3.eps

%approx
poly_3s=simplify(poly_3);
poly_3s=subs(poly_3s,'x','K');
poly_3s=subs(poly_3s,'y','H');
Z3=simplify(poly_3s);
Z3 =(H./5 - 1).*(K./450 - 1./90) + (3080777541025067.*((2.*H.^2)./25 - (4.*H)./5 + 1).*((2.*K.^2)./25 - (4.*K)./5 + 1))./73786976294838206464 + (4951920775910657.*((4.*H.^3)./125 - (12.*H.^2)./25 + (9.*H)./5 - 1).*((4.*K.^3)./125 - (12.*K.^2)./25 + (9.*K)./5 - 1))./37778931862957161709568 + 1/4;
[K,H] = meshgrid(0:10, 0:10);
[C,h] = contour(K,H,Z3);

print -dpdf ex5d2_3.eps

%approx vs function contour
[K,H] = meshgrid(0:10, 0:10);
[C,h] = contour(K,H,Z3);
hold on
[C,h] = contour(K,H,Z);
print -dpdf ex5d3_3.eps

%error contour
e=simplify(abs(f_e-poly_3));
e=subs(e,'x','K');
e=subs(e,'y','H');
e=subs(e,'k','K');
e=subs(e,'h','H');
e = abs((H./5 - 1).*(K./450 - 1./90) + (3080777541025067.*((2.*H.^2)./25 - (4.*H)./5 + 1).*((2.*K.^2)./25 - (4.*K)./5 + 1))./73786976294838206464 + (4951920775910657.*((4.*H.^3)./125 - (12.*H.^2)./25 + (9.*H)./5 - 1).*((4.*K.^3)./125 - (12.*K.^2)./25 + (9.*K)./5 - 1))./37778931862957161709568 - (633825300114114700748351602688.*H.*K)./(H.^(1./99) + K.^(1./99)).^99 + 1/4);
[K,H] = meshgrid(0:10, 0:10);
[C,h] = contour(K,H,e);
hold on
[C,h] = contour(K,H,Z3);
print -dpdf ex5d4_3.eps
hold off

