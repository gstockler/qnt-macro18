%% PART II.2 - The infnitely-lived households economy
% Quadratic utility

%% Setting up the model
clear
close all;

% parameters of the model
r = .04;
rho =.06;
beta = 1/(1+rho);
cbar = 100;
sigma =2;
sigma_y = .1; 
gamma =0;

% normalization
w=1;

% Markov Process (need to make it more general)
y=[1-sigma_y;1+sigma_y];
Pi=[(1+gamma)/2,(1-gamma)/2; (1-gamma)/2,(1+gamma)/2];
ny=length(y);

% functional forms
%u1  = @(c) c.^(1-sigma)/(1-sigma)-1/(1-sigma);   % preferences
u2  = @(c) .5.*(c-cbar).^2 ;  % preferences
%% A) LIQUIDITY CONSTRAINTS
% State-space bounds
a_max=30;
% lower bound: depends on the borrowing constraint
% 1= natural limit (precautionary savings - prudence);
% else = no borrowing at all (liquidity constraints);
[a_min]=borrowing_const(2,y,r);

% Discretizing state space
na=100; %number of grids
agrid=zeros(na,1);
for i=1:na
    agrid(i)=a_min+(i-1).*(a_max-a_min)/(na-1);
    agrid(i,1)=agrid(i)';
end

% converting grid into matrices: amat varies along columns, apmat varies along rows.
[amat,apmat] = ndgrid(agrid,agrid);     % ap denotes a'

%%
%initial guess
V1=zeros(na,1);
V2=zeros(na,1);

% Return matrices
% consumption
cmat1 = w*y(1,1)+(1+r)*amat-apmat;
cmat1(cmat1<0)=[-99999999];
cmat2 = w*y(2,1)+(1+r)*amat-apmat;
cmat2(cmat2<0)=[-99999999];
    
%utility
%U1mat_1=u1(cmat1);
%U1mat_2=u1(cmat2);

U2mat_1=u2(cmat1);
U2mat_2=u2(cmat2);

%% Iterations
% for storing
V= repmat(0,na,ny);
ga= repmat(0,na,ny);

maxits = 400;
tol = 0.001;
dif = .5;
its=0;
%tme = cputime;
tic

while dif > tol & its<maxits

  [V1_its,a1_its]=max(U2mat_1 + beta*repmat(V*Pi(1,:)',1,na));
  [V2_its,a2_its]=max(U2mat_2 + beta*repmat(V*Pi(2,:)',1,na));
  
  a_its=[a1_its' a2_its'];
  V_its=[V1_its' V2_its'];
  
  dif=max(abs(V_its-V));
  V=V_its;
  ga=agrid(a_its);
  its = its+1;
end;

%cputime-tme
fprintf('Elapsed Time = %4.2f Seconds\n',toc);
disp('Iterations until convergence (value function fixed point)')
disp(its);
%toc

% consumption policy function
a=agrid;
gc1 = w*y(1,1)+(1+r)*agrid-ga(:,1);
gc2 = w*y(2,1)+(1+r)*agrid-ga(:,2);


%% Policy function and solution: plots
figure(1)
plot(agrid,V(:,1),'r');
hold on
plot(agrid,V(:,2),'b');
title('Value function')
xlabel('a')
ylabel('V')
legend('High-income', 'Low-income')
hold off
%print -dpdf q1_1.eps

figure(2)
plot(a,a,'--y');
hold on
plot(agrid,ga(:,1),'r');
hold on
plot(agrid,ga(:,2),'b');
title('Policy function: assets')
xlabel('a')
ylabel('ga')
legend('45-line','High-income', 'Low-income')
hold off

%print -dpdf q1_2.eps

figure(3)
plot(a,a,'--y');
hold on
plot(a,gc1,'r');
hold on
plot(a,gc2,'b');
title('Policy function: consumption')
xlabel('a')
ylabel('gc')
legend('45-line','High-income', 'Low-income')
hold off

% -------------------------------------------------------------------------

%% Setting up the model
clear
close all;

% parameters of the model
r = .04;
rho =.06;
beta = 1/(1+rho);
cbar = 100;
sigma =2;
sigma_y = .1; 
gamma =0;

% normalization
w=1;

% Markov Process (need to make it more general)
y=[1-sigma_y;1+sigma_y];
Pi=[(1+gamma)/2,(1-gamma)/2; (1-gamma)/2,(1+gamma)/2];
ny=length(y);

% functional forms
%u1  = @(c) c.^(1-sigma)/(1-sigma)-1/(1-sigma);   % preferences
u2  = @(c) .5.*(c-cbar).^2 ;  % preferences

%% B) PRECAUTIONARY SAVINGS
% State-space bounds
a_max=30;
% lower bound: depends on the borrowing constraint
% 1= natural limit (precautionary savings - prudence);
% else = no borrowing at all (liquidity constraints);
[a_min]=borrowing_const(1,y,r);

% Discretizing state space
na=100; %number of grids
agrid=zeros(na,1);
for i=1:na
    agrid(i)=a_min+(i-1).*(a_max-a_min)/(na-1);
    agrid(i,1)=agrid(i)';
end

% converting grid into matrices: amat varies along columns, apmat varies along rows.
[amat,apmat] = ndgrid(agrid,agrid);     % ap denotes a'

%%
%initial guess
V1=zeros(na,1);
V2=zeros(na,1);

% Return matrices
% consumption
cmat1 = w*y(1,1)+(1+r)*amat-apmat;
cmat1(cmat1<0)=[-99999999];
cmat2 = w*y(2,1)+(1+r)*amat-apmat;
cmat2(cmat2<0)=[-99999999];
    
%utility
%U1mat_1=u1(cmat1);
%U1mat_2=u1(cmat2);

U2mat_1=u2(cmat1);
U2mat_2=u2(cmat2);

%% Iterations
% for storing
V= repmat(0,na,ny);
ga= repmat(0,na,ny);

maxits = 400;
tol = 0.001;
dif = .5;
its=0;
%tme = cputime;
tic

while dif > tol & its<maxits

  [V1_its,a1_its]=max(U2mat_1 + beta*repmat(V*Pi(1,:)',1,na));
  [V2_its,a2_its]=max(U2mat_2 + beta*repmat(V*Pi(2,:)',1,na));
  
  a_its=[a1_its' a2_its'];
  V_its=[V1_its' V2_its'];
  
  dif=max(abs(V_its-V));
  V=V_its;
  ga=agrid(a_its);
  its = its+1;
end;

%cputime-tme
fprintf('Elapsed Time = %4.2f Seconds\n',toc);
disp('Iterations until convergence (value function fixed point)')
disp(its);
%toc

% consumption policy function
a=agrid;
gc1 = w*y(1,1)+(1+r)*agrid-ga(:,1);
gc2 = w*y(2,1)+(1+r)*agrid-ga(:,2);


%% Policy function and solution: plots
figure(1)
plot(agrid,V(:,1),'r');
hold on
plot(agrid,V(:,2),'b');
title('Value function')
xlabel('a')
ylabel('V')
legend('High-income', 'Low-income')
hold off
print -dpdf pe2_1.eps

figure(2)
plot(a,a,'--y');
hold on
plot(agrid,ga(:,1),'r');
hold on
plot(agrid,ga(:,2),'b');
title('Policy function: assets')
xlabel('a')
ylabel('ga')
legend('45-line','High-income', 'Low-income')
hold off

print -dpdf pe2_2.eps

figure(3)
plot(a,a,'--y');
hold on
plot(a,gc1,'r');
hold on
plot(a,gc2,'b');
title('Policy function: consumption')
xlabel('a')
ylabel('gc')
legend('45-line','High-income', 'Low-income')
hold off

print -dpdf pe2_3.eps
% -------------------------------------------------------------------------