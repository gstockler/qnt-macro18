%% HOWARD'S POLICY ITERATION

clear
close all;
%parameters of the model
theta=.679;
beta=.988;
delta=.013;
kappa= 5.24; 
nu= 2;

%% STEP 1: Discretize state space
k_ss=1;
k_max=1.2;
k_min=.2;
p=100; % number of grid-points
eta=(k_max-k_min)/(p-1); %'step'
k=zeros(p,1);
for i=1:p
    k(i)=k_min+(i-1).*eta;
    k(i,1)=k(i)';
end

%% STEP 2: return matrix M
M=zeros(p,p);
C = (k.^(1-theta))-k' +(1-delta).*k;
for i=1:p
    for j=1:p
        if C(i,j)<0
            M(i,j)=-9999999;
        else
            M(i,j)=log(C(i,j));
        end
    end
end

%% Howards' 1
%maxits = 100;
tol = 0.001;
dif = .5;
V=zeros(p,1); %initial guess for the value function
T=M+beta*V';
[V_hmax]=max(T,[],2);

%new
ph=50; 
tStart1 = tic;
its = 0;
while dif > tol 
    V=V_hmax;
    for i=1:p
        for j=1:p
            T(i,j)=M(i,j)+beta*V(j,1);
        end
    end
    [V_max,ind]=max(T,[],2);
    gk1=k(ind);
    gc1=k.^(1-theta)+(1-delta)*k-gk1;
    % Howards' additional step: ph is the new integer-loop
    for h=2:ph
        Vh(:,1)=V_max;
        for i=1:p
            Vh(i,h)= log(k(i,1)^(1-theta)+(1-delta)*k(i,1)-gk1(i,1))+beta*Vh(i,h-1);
        end  
    end
    V_hmax=Vh(:,ph);
    dif = max(abs(V_hmax-V));
    its=its+1;
end
tElapsed1 = toc(tStart1);
disp(its)
disp(tElapsed1)


%%
%new
ph=10;
tStart2 = tic;
its = 0;
while dif > tol 
    V=V_hmax;
    for i=1:p
        for j=1:p
            T(i,j)=M(i,j)+beta*V(j,1);
        end
    end
    [V_max,ind]=max(T,[],2);
    gk1=k(ind);
    gc1=k.^(1-theta)+(1-delta)*k-gk1;
    % Howards' additional step: ph is the new integer-loop
    for h=2:ph
        Vh(:,1)=V_max;
        for i=1:p
            Vh(i,h)= log(k(i,1)^(1-theta)+(1-delta)*k(i,1)-gk1(i,1))+beta*Vh(i,h-1);
        end
    end
    V_hmax=Vh(:,ph);
    dif = max(abs(V_hmax-V));
    its=its+1;
end
tElapsed2 = toc(tStart2);
disp(its)
disp(tElapsed2)

%%
ph=99;
tStart3 = tic;
its = 0;
while dif > tol 
    V=V_hmax;
    for i=1:p
        for j=1:p
            T(i,j)=M(i,j)+beta*V(j,1);
        end
    end
    [V_max,ind]=max(T,[],2);
    gk1=k(ind);
    gc1=k.^(1-theta)+(1-delta)*k-gk1;
    % Howards' additional step: ph is the new integer-loop
    for h=2:ph
        Vh(:,1)=V_max;
        for i=1:p
            Vh(i,h)= log(k(i,1)^(1-theta)+(1-delta)*k(i,1)-gk1(i,1))+beta*Vh(i,h-1);
        end
    end
    V_hmax=Vh(:,ph);
    dif = max(abs(V_hmax-V));
    its=its+1;
end
tElapsed3 = toc(tStart3);
disp(its)
disp(tElapsed3)


