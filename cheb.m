%% Chebyshev regression algorithm to approximate the value function
% cheb poly order: 10
% nodes: 20

clear
close all;

%parameters of the model
theta=.679;
beta=.988;
delta=.013;
kappa= 5.24; 
nu= 2;

tol = 0.001; %tolerance VFI
maxits=100;

% cheb regression definitios
n=10; %cheb poly order
m= 20; % number of nodes
p=m;

%% STEP 1.1: Nodes/Knots
% call function interpolation_cheb: for cheb nodes
a=0;
b=50;
[x]=interpolation_cheb(0,50,m);

%% STEP 1.2: Discretize state space
%k_ss=1;
k_ss=(((1-theta).*beta)./(1-beta+beta.*delta)).^(1./theta); 
k_max=1.25*k_ss;
k_min=.25*k_ss;
%p=100; % number of grid-points
eta=(k_max-k_min)/(p-1); %'step'
k=zeros(p,1);
for i=1:p
    k(i)=k_min+(i-1).*eta;
    k(i,1)=k(i)';
end

% converting grid into [-1,1];
kx(:,1)=(k(:,1)-k_min)*(2/(k_max-k_min))-1;

%% STEP 3: Initial guess for the Cheb Poly
% intial guess for the value function on the nodes x
% let it be ki=kj 
% V0=1/(1-beta)*log(x.^(1-theta)-delta.*x); % number of nodes (m)

% initial guess for the coefs
% c0(1,1)=1/m*(sum(V0)); % n=0
% for i=1:n
%    psi(:,i)=cos(i*acos(2*x));
% end
c0(1,2:n+1)=(2/m)*V0*psi; % n=1..10

% then, the implied value function is (using cheb algorithm)
syms K 
psik=cos((0:n).*acos(2*((K-a)/(b-a))-1));
s=c0(1,(1:n+1)).*psik(1,(1:n+1));
poly=sum(s);

%% STEP 4: policy function - argmax
U=k.^(1-theta)+(1-delta)*k-K;
T=U+beta*poly;
for i=1:p
    
     g(i,1)=fminbnd(matlabFunction(T(i,1)),k_min,k_max);
     
end

%% STEP 5: Updating value function
K=g;
Vg=eval(poly);
V=real(log(k.^(1-theta)+(1-delta)*k-g)+beta*Vg); % new "V0"

%% STEP 6: obtain updated coef.
c(1,1)=1/m*(sum(V)); % n=0
for i=1:n
    psi(:,i)=cos(i*acos(2*x));
end
c(1,2:n+1)=(2/m)*V'*psi; % n=1..10

%% STEP 7: Comparing coef.
dif=max(abs(c-c0));
% if dif< tol: done!
% otherwise, redo from step 3 with the updated coef.

%% LOOP
% initial guesses
V0=1/(1-beta)*log(x.^(1-theta)-delta.*x); % number of nodes (m)

c(1,1)=1/m*(sum(V0)); % n=0
for i=1:n
    psi(:,i)=cos(i*acos(2*x));
end
c(1,2:n+1)=(2/m)*V0*psi; % n=1..10

% then, the implied value function is (using cheb algorithm)
syms K 
psik=cos((0:n).*acos(2*((K-a)/(b-a))-1));
s=c(1,(1:n+1)).*psik(1,(1:n+1));
poly=sum(s);
its=0;
dif=.5;

tStart_cheb = tic;
while dif > tol & its<maxits
    syms K 
    c0=c;
    psik=cos((0:n).*acos(2*((K-a)/(b-a))-1));
    s=c0(1,(1:n+1)).*psik(1,(1:n+1));
    poly=sum(s);
    
    U=k.^(1-theta)+(1-delta)*k-K;
    T=U+beta*poly;
    for i=1:p
        g(i,1)=fminbnd(matlabFunction(T(i,1)),k_min,k_max);
    end
    K=g;
    Vg=eval(poly);
    V=real(log(k.^(1-theta)+(1-delta)*k-g)+beta*Vg); % new "V0"
    
    c(1,1)=1/m*(sum(V)); % n=0
    for i=1:n
        psi(:,i)=cos(i*acos(2*x));
    end
    c(1,2:n+1)=(2/m)*V'*psi; % n=1..10
    
    dif=max(abs(c-c0));
    its=its+1;
end
tElapsed_cheb = toc(tStart_cheb);
disp(its)
disp(tElapsed_cheb)
V_opt=V;
g_opt=g;% are used to constructoptimal decision rule function
    

 





