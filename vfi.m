%% PART 1 - VFI
%% A) Brute-force VFI
% The model is:
% ------- Bellman Equation
% V(k)= max_k' ln[k^(1-theta)+(1-delta)k-k']-kappa(nu/1+nu)+ beta V(k')
% ------- Euler Equation
% 1/c=1/c' beta [(1-theta)k^(-theta)+1-delta]
clear
close all;
%parameters of the model
theta=.679;
beta=.988;
delta=.013;
kappa= 5.24; 
nu= 2;

%% STEP 1: Discretize state space
% first, since we know the steady-state,
% let's use it to get a sense of the bounds for k
% k_ss=1;
k_ss=(((1-theta).*beta)./(1-beta+beta.*delta)).^(1./theta); 
k_max=1.2*k_ss;
k_min=.2*k_ss;
p=100; % number of grid-points
eta=(k_max-k_min)/(p-1); %'step'
k=zeros(p,1);
for i=1:p
    k(i)=k_min+(i-1).*eta;
    k(i,1)=k(i)';
end

%% STEP 3: return matrix M
M=zeros(p,p);
c=zeros(p,p);
for i=1:p
    for j=1:p
        c(i,j)=k(i,1).^(1-theta)+(1-delta).*k(i,1)-k(j,1);
        if c(i,j)<0
            M(i,j)=-9999999;
        else
            M(i,j)=log(c(i,j));
        end
    end
end

%% STEP 4:
maxits = 300;
tol = 0.001;
dif = .5;

V=zeros(p,1); %initial guess for the value function
T=M+beta*V';
[V_max,ind]=max(T,[],2);

tStart1 = tic;
its = 0;
while dif > tol & its<maxits
    V=V_max;
    for i=1:p
        for j=1:p
            T(i,j)=M(i,j)+beta*V(j,1);
        end
    end
    [V_max,ind]=max(T,[],2);
    gk1=k(ind);
    gc1=k.^(1-theta)+(1-delta)*k-gk1;
    %
    dif = max(abs(V_max-V));
    its=its+1;
end
tElapsed1 = toc(tStart1);
disp(its)
disp(tElapsed1)

%% Policy function and solution: plots
figure(1)
plot(k,V_max);
title('Value function (100 iterations)')
xlabel('k')
ylabel('V')

print -dpdf q1_1a.eps

figure(2)
plot(k,gk1,'r');
hold on
hline=refline(1,0)
hline.Color = 'b';
hline.LineStyle=':';
title('Capital policy function: g_k')
xlabel('k')
ylabel('g_k')

print -dpdf q1_1b.eps

figure(3)
plot(k,gc1)
title('Consumption policy function: g_c')
xlabel('k')
ylabel('g_c')

print -dpdf q1_1c.eps
% -------------------------------------------------------------------------