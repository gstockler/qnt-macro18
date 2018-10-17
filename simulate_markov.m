function [chain,state] = simulate_markov(x,P,pi0,T);

%x     = states
%P     = Markov transition matrix
%  pi0   = probability distribution over initial state
%  T     = number of periods to simulate

n = length(x); %% what is the size of the state vector?
E = rand(T,1); %% T-vector of draws from independent uniform [0,1]  

cumsumP = P*triu(ones(size(P)));

E0   = rand(1,1);
ppi0 = [0,cumsum(pi0)];
s0   = ((E0<=ppi0(2:n+1)).*(E0>ppi0(1:n)))';
s    = s0; 

%%%%% ITERATE ON THE CHAIN

for t=1:T,
    state(:,t) = s;
    ppi        = [0,s'*cumsumP];
    s          = ((E(t)<=ppi(2:n+1)).*(E(t)>ppi(1:n)))';
end

chain = x'*state;    


