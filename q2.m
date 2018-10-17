%% QNT MACRO - PS4
%-------------QUESTION 2 - Business Cycle Fluctuations--------------------%
%% I) Stochastic VFI
clear; close all; clc;

%parameters of the model
theta=.679;
beta=.988;
delta=.013;
kappa= 5.24; 
nu= 2;

%% Stochastic productivity shocs: z
% transition matrix
z=[1/1.01; 1.01]; % it can take 2 possible values
syms x
prob=double(solve(x*.9901+(1-x)*1.0100-1));
Pi=[prob, 1-prob; 1-prob, prob];
q=2;

%% STEP 1: Discretize state space
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

%% STEP 3: return matrix M (it does not change)
M=zeros(p,p);
c=zeros(p,p);
for l=1:q
    for i=1:p
        for j=1:p
            c((p*(l-1)+i),j)=z(l,1).*k(i,1).^(1-theta)+(1-delta).*k(i,1)-k(j,1);
            if c((p*(l-1)+i),j)<0
                M((p*(l-1)+i),j)=-9999999;
            else
                M((p*(l-1)+i),j)=log(c((p*(l-1)+i),j));
            end
        end
    end
end

%% STEP 4: guess V0
% for all possible combinations of z and k
V0=zeros(p*q,1);

%% STEP 5.1: Expectation Matrix (at each iteration)
W1=Pi(1,1)*V0(1:p*q/2,1)'+Pi(1,2)*V0(p*q/2+1:p*q,1)';
W2=Pi(2,1)*V0(1:p*q/2,1)'+Pi(2,2)*V0(p*q/2+1:p*q,1)';
W=[W1;W2];

%% 5.2: Chi-matrix (for each iteration)
for l=1:q
    for i=1:p
        for j=1:p
            chi(p*(l-1)+i,j)=M(p*(l-1)+i,j)+beta*W(l,j);
        end
    end
end

[V_max,ind]=max(chi,[],2);


%% STEP 4:
maxits = 400;
tol = 0.001;
dif = .5;

tStart = tic;
its = 0;
while dif > tol & its<maxits
    V=V_max;
    W1=Pi(1,1)*V(1:p*q/2,1)'+Pi(1,2)*V(p*q/2+1:p*q,1)';
    W2=Pi(2,1)*V(1:p*q/2,1)'+Pi(2,2)*V(p*q/2+1:p*q,1)';
    W=[W1;W2];
    for l=1:q
        for i=1:p
            for j=1:p
                chi(p*(l-1)+i,j)=M(p*(l-1)+i,j)+beta*W(l,j);
            end
        end
    end

    [V_max,ind]=max(chi,[],2);
    gk=k(ind);
    gc=repmat(k.^(1-theta)+(1-delta)*k,2,1)-gk;

    dif = max(abs(V_max-V));
    its=its+1;
end
tElapsed = toc(tStart);
disp(its)
disp(tElapsed)

%% Plots
% for V for the first value shock z1
 plot(V_max(1:100,1))
hold on
% for V for the second value shock z2
plot(V_max(101:200,1))

print -dpdf q2_1.eps

%% Simulation

t=500;
[chain,state] = simulate_markov(z,Pi,Pi(1,:),t);
kt=k_ss*z(1,1);
for i=2:t
    r=find(abs(gk-kt)<.1);
    if state(1,i)==1 & length(r)>=1
        r=r(1,1);
    elseif length(r)>1
        r=r(2,1);
    end
    kp(i,1)=gk(r);
    yt(i,1)=kt^(1-theta)*chain(1,i-1);
    it(i,1)=kp(i,1)-(1-delta)*kt;
    ct(i,1)=yt(i,1)-it(i,1);
    kt=kp(i,1)*chain(1,i);
end

%%
%plot(yt(100:500));
%plot(kp(100:500));

c_l=log(ct(100:500));
y_l=log(yt(100:500));
i_l=log(it(100:500));
k_l=log(kp(100:500));

yvar=var(y_l);
cvar=var(c_l);
ivar=var(i_l);
kvar=var(k_l);

y_hp=hpfilter(yt(100:500),1600);
yvar_hp=var(y_hp);
i_hp=hpfilter(it(100:500),1600);
ivar_hp=var(i_hp);
c_hp=hpfilter(ct(100:500),1600);
cvar_hp=var(c_hp);
k_hp=hpfilter(kp(100:500),1600);
kvar_hp=var(k_hp);

disp(yvar)
disp(yvar_hp)

disp(cvar)
disp(cvar_hp);

disp(ivar)
disp(ivar_hp);

disp(kvar)
disp(kvar_hp);

%% III) IRF
%deviation from ss

y_ss=k_ss^(1-theta);
i_ss=delta*k_ss;
c_ss=y_ss-i_ss;


shock=z(1,1);
kk(1)=k_ss;
iv(1)=i_ss;
y(1)=shock*(k(1)^(1-theta));
c(1)=y(1)-i(1);
t=50;
for i=2:t
    kk(i)=y(i-1)-c(i-1)+(1-delta)*kk(i-1);
    y(i)=kk(i)^(1-theta);
    c(i)=c(i-1)*beta*[(1-theta)*kk(i)^(-theta)+1-delta]; %EE
    iv(i)=y(i)-c(i);
end

for i=1:t
    kd(i)=log(kk(i))-log(k_ss);
    cd(i)=log(c(i))-log(c_ss);
    yd(i)=log(y(i))-log(y_ss);
    ivd(i)=log(iv(i))-log(i_ss);
end




