function [TS_Approx, TS_plot] = taylorSeries(Fun,a,N,xl,xu)
syms x
TS_Approx = subs(Fun,a); 
for n = 1:N
    derivative = diff(Fun,n); 
    TS_Approx = TS_Approx +subs(derivative,a)*(x-a)^n/factorial(n);
end
TS_plot=fplot(TS_Approx, [xl xu]);
end

