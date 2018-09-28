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