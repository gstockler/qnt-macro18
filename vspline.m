function [Tv] = vspline(a)
% Calculates the spline of the value function evaluated at some asset in
% the state space
% Use to get policy function for asset, since ga=argmax
% For CRRA utility -> need to adapt!


global beta w r a0 sigma a_nk y0 theta1 theta2 
%ind=piecewise_locate(a_nk,a);

ind=max(sum(a>a_nk),1);% location of asset on knot-subintervals
b0=theta1(ind); %spline intercept
b1=theta2(ind); %spline slope
%spline= @(a) V0(ind)-a_nk(ind)*(V0(ind+1)-V0(ind))/(a_nk(ind+1)-a_nk(ind))+(V0(ind+1)-V0(ind))/(a_nk(ind+1)-a_nk(ind))*a;
spline=b0+b1*a;
c_nk=w*y0+(1+r)*a0-a; %consumption
Tv = (c_nk.^(1-sigma)-1)./(1-sigma)+beta*spline; %value funvtion T-mapping
Tv =- Tv; %since will use minimization command
