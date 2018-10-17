%% D) Speeding up: local search on the decision rule

clear
close all;
%parameters of the model
theta=.679;
beta=.988;
delta=.013;
kappa= 5.24; 
nu= 2;

%% STEP 1: Discretize state space
k_ss=(((1-theta).*beta)./(1-beta+beta.*delta)).^(1./theta); 
k_max=1.2;
k_min=.2;
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
%%
% STEP 2.5 (new): characterizing local search
s_min=.5;
s_max=.5;

maxits = 3;
tol = 0.01;
dif = tol+1000;


V=zeros(p,1); %initial guess for the value function
[V_max,ind]=max(M,[],2);
gk4=k;

for i=2:p
    k_low(1,1)=k(1);
    k_low(i,1)=gk4(i-1,1)-s_min;
    k_high(1,1)=k(p);
    k_high(i,1)=gk4(i-1,1)+s_max;

end

k_fine=find(k_low(1,1)<k<k_high(1,1));
for i=1:p
    k_fine=find(k>k_low(i,1) & k<k_high(i,1));
    l_grid=k(k_fine,1);
    for j=1:length(l_grid)
        c(i,j)=l_grid(i,1).^(1-theta)+(1-delta).*l_grid(i,1)-l_grid(j,1);
        if c(i,j)<0
            M(i,j)=-9999999;
        else
            M(i,j)=log(c(i,j));
        end
        T=M(i,j)+beta*V(j,1);
        V_max(i,1)=max(T,[],2);
    
    end
    %V_max(i,1)=max(T(i,:),[],2);
    
end


% starting iterations
tStart4 = tic;
for its=1:maxits
    while dif>tol
      V=V_max;
    
      for i=2:p
            j_low(1,1)=1;
            j_high(1,1)=p;
             for j=j_low(1,1):j_high(1,1)
          
              T(1,j)=M(1,j)+beta*V(j,1);
          end
          j_low(i,1)=max([1,gk4(i-1,1)-s_min]);
          j_high(i,1)=min([1,gk4(i-1,1)+s_max]);
          for i=2:p
          for j=j_low(i,1):j_high(i,1)
              T(i,j)=M(i,j)+beta*V(j,1);
          end
          end
      end
      [V_max,ind]=max(T,[],2);
      gk4 =k(ind,1); %policy function
      gc4 =k.^(1-theta)+(1-delta)*k-gk4; % policy function
      dif =abs(V_max-V);
   
    end  
end
 tElapsed4 = toc(tStart4);