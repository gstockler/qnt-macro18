%% CONCAVITY

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

%% STEP 3: Iterations

maxits = 200;
tol = 0.001;
dif = tol+1000;
its=0;

V=zeros(p,1);
T0=M+beta*V';
[V_max]=max(T0,[],2);

tStart3 = tic;
%for its=1:maxits
 while dif > tol
        V=V_max;
        for i=1:p
            for j=1:p-1
                if M(i,j)+beta*V(j,1)> M(i,j+1)+beta*V(j+1,1)
                    T(i,j+1)= M(i,j)+beta*V(j,1);
                    %K(i,1)=j-1;
                end
            end
        end
        %for i=1:p
         %   for j=2:p
          %      if M(i,j-1)+beta*V(j-1,1)> M(i,j)+beta*V(j,1)
           %         T(i,j)= M(i,j-1)+beta*V(j-1,1);
                    %K(i,1)=j-1;
            %    end
            %end
        %end
        [V_max,ind]=max(T,[],2);
        gk3=k(ind,1);
        gc3 =k.^(1-theta)+(1-delta).*k-gk3; % policy function
        dif =abs(V_max-V);
        its=its+1;
 end
%end
 tElapsed3 = toc(tStart3);
 disp(its)
 disp( tElapsed3)
 
 %% TRY CONC+MONT
 %% STEP 3
% one g_low for each i
g_low=k;
for i=2:p
    for j=1:p
        if k(j,1)>=g_low(i-1,1)
            T(i,j)=M(i,j)+beta*V(j,1);
        end
    end
end
i=2
V=zeros(p,1);
for j=1:p
        if k(j,1)>=g_low(i-1,1)
            T(i,j)=M(i,j)+beta*V(j,1);
        end
    end
        
    


%% STEP 3:
maxits = 300;
tol = 0.01;
dif = tol+1000;

g_low=repmat(k(1,1),100,1);
[ind_low] = find(k-g_low > 0);
V=zeros(p,1);
T=M(:,ind_low)+beta*V(ind_low,1)';
[V_max]=max(T,[],2);


tStart4 = tic;
for its=1:maxits
    while dif > tol
        V=V_max;
       [ind_low] = find(k-g_low > 0);
       k=k(ind_low,1);
       C = (k.^(1-theta))-k' +(1-delta).*k;
       for i=1:length(ind_low) 
           for j=1:length(ind_low) 
               if C(i,j)<0
                 M(i,j)=-9999999;
               else
                M(i,j)=log(C(i,j));
               end
           end
       end

       for i=1:length(ind_low)
           for j=2:length(ind_low)
                if M(i,j-1)+beta*V(j-1,1)> M(i,j)+beta*V(j,1)
                    T(i,1)= M(i,j-1)+beta*V(j-1,1);
                    K(i,1)=j-1;
                end
            end
       end
        V_max=T;
        gk4=k(K,1);
        gc4 =k.^(1-theta)+(1-delta).*k-gk4; % policy function
        dif =abs(V_max-V);
    end
end
 tElapsed4 = toc(tStart4);