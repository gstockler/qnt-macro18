%% Question 2 - no tax
clear all;
clc

% productivity parameter
eta=[1, 1.5, 2.5, 3];

% Random initial endowment
y0=.001 + (.009-.001).*rand(100,1);
y0(y0>=0.0055 & y0<=0.0087)=.001;

%y0(y0<=0.0055)=0.001; 
%y0(y0>=0.0087)=0.001;

% Guess for the interest rate
r= 0.004;

% Setting parameters values
sigma=3;
k=4;
v=4;
beta=0.99;
tau=0;
T0=0;
T1=0;

%% for least-productive workers: eta_y =1 (z=1)
sol1=zeros(100,7);
for n=1:100
      i=y0(n,1); %random initial endowment
      j=eta(1,1); %productivity
      %randoming the shock
      s = randi([1, 2], 1); 
      e_y=[-.05;.05]; 
      eps=e_y(s,1); 
      %system of equations based on FOCs and borrowing constraint: unkowns
      %a, h0, h1, lambda
       f=@(a,h0,h1,lambda)[ ((1-tau)*j*h0 + i + T0 -a)^(-sigma)*(1-tau)*j-k*(h0)^(1/v);
            beta*(((1-tau)*j*h1)+(1+r)*a + T1)^(-sigma)*(1-tau)*j-k*(h1)^(1/v);
          beta*(((1-tau)*j*h1)+(1+r)*a + T1)^(-sigma)*(1+r)-lambda-((1-tau)*j*h0 + i + T0 -a)^(-sigma);
           ((1-tau)*j*h0 + i + T0 -a) + (1/(1+r))*((1-tau)*(j+eps)*h1 + (1+r)*a + T1) - i - (1+r)*(j+eps)*h1];
       fun=@(u) f(u(1),u(2),u(3),u(4));
      xGuess= ([0.05,0.1,0.1, 1]);
      u1 = fsolve(fun, xGuess) ;
      sol1(n,1:7)=[u1,i,j,eps];
 
end

%% for second-least-productive workers: eta_y =1.5 (z=2)
sol2=zeros(100,7);
for n=1:100
   
      i=y0(n,1);
      j=eta(1,2);
      s = randi([1, 2], 1);
      e_y=[-.05;.05];
      eps=e_y(s,1);
            %a = x(1);
            %h0 = x(2);
            %h1 = x(3);
            %lambda = x(4);   
           
       f=@(a,h0,h1,lambda)[ ((1-tau)*j*h0 + i + T0 -a)^(-sigma)*(1-tau)*j-k*(h0)^(1/v);
            beta*(((1-tau)*j*h1)+(1+r)*a + T1)^(-sigma)*(1-tau)*j-k*(h1)^(1/v);
          beta*(((1-tau)*j*h1)+(1+r)*a + T1)^(-sigma)*(1+r)-lambda-((1-tau)*j*h0 + i + T0 -a)^(-sigma);
           ((1-tau)*j*h0 + i + T0 -a) + (1/(1+r))*((1-tau)*(j+eps)*h1 + (1+r)*a + T1) - i - (1+r)*(j+eps)*h1];
       fun=@(u) f(u(1),u(2),u(3),u(4));
      xGuess= ([0.05,0.1,0.1, 1]);
      u1 = fsolve(fun, xGuess) ;
      sol2(n,1:7)=[u1,i,j,eps];
 
end

%% for second-most-productive workers: eta_y =2.5 (z=3)
sol3=zeros(100,7);
for n=1:100
      i=y0(n,1);
      j=eta(1,3);
      s = randi([1, 2], 1);
      e_y=[-.05;.05];
      eps=e_y(s,1);   
       f=@(a,h0,h1,lambda)[ ((1-tau)*j*h0 + i + T0 -a)^(-sigma)*(1-tau)*j-k*(h0)^(1/v);
            beta*(((1-tau)*j*h1)+(1+r)*a + T1)^(-sigma)*(1-tau)*j-k*(h1)^(1/v);
          beta*(((1-tau)*j*h1)+(1+r)*a + T1)^(-sigma)*(1+r)-lambda-((1-tau)*j*h0 + i + T0 -a)^(-sigma);
           ((1-tau)*j*h0 + i + T0 -a) + (1/(1+r))*((1-tau)*(j+eps)*h1 + (1+r)*a + T1) - i - (1+r)*(j+eps)*h1];
       fun=@(u) f(u(1),u(2),u(3),u(4));
      xGuess= ([0.05,0.1,0.1, 1]);
      u1 = fsolve(fun, xGuess) ;
      sol3(n,1:7)=[u1,i,j,eps];
 
end

%% for most-productive workers: eta_y =3 (z=4)
sol4=zeros(100,7);
for n=1:100
      i=y0(n,1);
      j=eta(1,4);
      s = randi([1, 2], 1);
      e_y=[-.05;.05];
      eps=e_y(s,1);  
       f=@(a,h0,h1,lambda)[ ((1-tau)*j*h0 + i + T0 -a)^(-sigma)*(1-tau)*j-k*(h0)^(1/v);
            beta*(((1-tau)*j*h1)+(1+r)*a + T1)^(-sigma)*(1-tau)*j-k*(h1)^(1/v);
          beta*(((1-tau)*j*h1)+(1+r)*a + T1)^(-sigma)*(1+r)-lambda-((1-tau)*j*h0 + i + T0 -a)^(-sigma);
           ((1-tau)*j*h0 + i + T0 -a) + (1/(1+r))*((1-tau)*(j+eps)*h1 + (1+r)*a + T1) - i - (1+r)*(j+eps)*h1];
       fun=@(u) f(u(1),u(2),u(3),u(4));
      xGuess= ([0.05,0.1,0.1, 1]);
      u1 = fsolve(fun, xGuess) ;
      sol4(n,1:7)=[u1,i,j,eps];
 
end

%% In the aggregate 
% concatenating all solutions in a matrix
sol=[sol1;sol2;sol3;sol4]; %sol=[a, h0, h1, lambda, y0, eta_y, eps]
a=sol(:,1);
h0=sol(:,2);
h1=sol(:,3);
eps_y=sol(:,7);
% Calculating consumption in both periods and savings rate:
% c0=(1-tau)*j*h0+i+T0-a,c1=(1-tau)*(j+eps)*h1+(1+r)*a+T1, sv=a/(i+j*h0*(1-tau))
%c0_f=(1-tau).*eta.h0+y0_f+T0-a;

c0=(1-tau).*sol(:,6).*sol(:,2)+sol(:,5)+T0-sol(:,1);
c1=(1-tau).*(sol(:,6)+sol(:,7)).*sol(:,3)+(1+r).*sol(:,1)+T1;
sv=sol(:,1)./(sol(:,5)+sol(:,6).*sol(:,2).*(1-tau));

%% Graphs
figure(1)
subplot(2,2,1)
plot(a)
title('optimal savings');
 
subplot(2,2,2)
plot(c0)
title('optimal consumption-1st period');

subplot(2,2,3)
plot(c1)
title('optimal consumption-2nd period');

subplot(2,2,4)
plot(y0)
title('initial income distribution (per productivity)');

print -dpdf q2_1.eps
