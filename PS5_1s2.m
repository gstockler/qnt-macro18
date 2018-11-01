%% VFI II: Continous Method - Piecewise linear splines
%% Setting up the model
clear
close all;

global beta sigma a_nk nk w r
% parameters of the model
r = .04;
rho =.06;
beta = 1/(1+rho);
cbar = 100;
sigma =2;
sigma_y = .1; 
gamma =0;

%global na ny w y Pi
% normalization
w=1;

% Markov Process (need to make it more general)
y=[1-sigma_y;1+sigma_y];
Pi=[(1+gamma)/2,(1-gamma)/2; (1-gamma)/2,(1+gamma)/2];
ny=length(y);


% State-space bounds
a_max=50;
% lower bound: depends on the borrowing constraint
% 1= natural limit, else is no borrowing at all
[a_min]=borrowing_const(1,y,r);

% For iteration:
maxits = 300;
tol = 0.001;
dif = .5;

%% Spline set-up
% Setting the number of knots
nk=20;

% knots: allocating them in the bounds
a_nk=linspace(a_min,a_max,nk);
%% Using spline function vspline.m
% Initial guess
global V_nk theta1 theta2 y0 a0
V_its=zeros(nk,1);

% using only first shocke
y0=y(1,1);

its=0;

while dif>tol & its<maxits
    V_nk=V_its;
    [theta1,theta2] = lincoef_spline(a_nk,V_nk);
   
for i=1:nk
    a0=a_nk(1,i);
    aprime = fminbnd(@vspline,a_min,a_max);
    V_its(i,1)=-vspline(aprime);
    ga(i,1)=aprime;
end
dif=max(abs(V_its-V_nk));
its=its+1;
     
end


