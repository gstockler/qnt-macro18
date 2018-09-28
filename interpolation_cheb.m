function [x]=interpolation_cheb(a,b,n)
x = (a+b)/2 + ((b-a)/2)*cos((pi/n)*((1:n)-(1/2))); 
end
