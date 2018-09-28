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
print -dpdf ex1.eps
 