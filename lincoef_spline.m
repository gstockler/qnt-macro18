function [theta1,theta2] = lincoef_spline(x,y)
% x: knots
% y: image of x - f(x)
n=length(x); %number of knots;
for i=1:n-1
    theta1(i)=y(i)-x(i)*(y(i+1)-y(i))/(x(i+1)-x(i));
    theta2(i)=(y(i+1)-y(i))/(x(i+1)-x(i));
end


end
