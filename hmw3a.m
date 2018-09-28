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