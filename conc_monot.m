%% CONCAVITY+MONOTONICITY

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
%% STEP 3: correct
% one g_low for each i
% first iteration, wth V0=0, has no lower bound since k-grid is ascending
V=zeros(p,1);
T0=M+beta*V';
[V_max,ind]=max(T0,[],2);
gk4=k(ind);

maxits = 300;
tol = 0.001;
dif = tol+1000;
tStart4 = tic;
its=0;
%for its=1:maxits
    while dif > tol & its<maxits
        V=V_max;
        gk=gk4;
         for j=2:p
                if k(j,1)>=k(1,1)
                    if M(1,j-1)+beta*V(j-1,1)> M(1,j)+beta*V(j,1)
                        T(1,j)=M(1,j)+beta*V(j,1);
                    end
                end
         end
            
        for i=2:p
            for j=2:p
                if k(j,1)>=gk(i-1,1)
                     if M(i,j-1)+beta*V(j-1,1)> M(i,j)+beta*V(j,1)
                        T(i,j)= M(i,j-1)+beta*V(j-1,1);
                     end
                    %T(i,j)=M(i,j)+beta*V(j,1);
                end
            end
        end
        [V_max,ind]=max(T,[],2);
        gk4=k(ind,1);
        gc4 =k.^(1-theta)+(1-delta)*k-gk4; % policy function
        dif =abs(V_max-V);
        its=its+1;
        %g_low=gk;
   % end
end
 tElapsed4 = toc(tStart4);
 disp(tElapsed4)
 disp(its)

