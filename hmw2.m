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
print -dpdf ex2.eps
 