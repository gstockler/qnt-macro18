%% Sovereign Debt & Default - Quant Macro II - PS3
% Gabriela Barbosa
%
%   This is based on Prof. Luis Rojas Code for the model in the framework
% of Aguiar and Amador (2014)
%__________________________________________________________________________
%% %%%%%%%%%%%%%%%%% Bond-Price Schedule Computation %%%%%%%%%%%%%%%%%%%%%%

%% 0. Model Setup & Parameters 
clear;
clc;

%%%%%%%%%%%%%%%%%% Discount factor 
betta=0.95;

%%%%%%%%%%%%%%%%%% Possible values for GDP
y_grid=[0.9, 1, 1.05];
%y_grid=[0.6, 1, 1.5];

%%%%%%%%%%%%%%%%%% Transition matrix for GDP
piy=[0.5, 0.3, 0.2;...
     0.1, 0.6, 0.3;...
     0.2, 0.4, 0.4];

%piy=[0.5, 0.3, 0.2;...
%     0.05, 0.65, 0.3;...
%     0.02, 0.55, 0.43];

%%%%%%%%%%%%%%%%%% Risk aversion
sig=2;
%sig=1.5;

%%%%%%%%%%%%%%%%%% Utility function $u(c)=\frac{c^{1-\sigma}}{1-\sigma}$
u = @(c) c.^(1-sig)./(1-sig);

%%%%%%%%%%%%%%%%%% Risk-free interest rate
R=1;

%%%%%%%%%%%%%%%%%% Probability of regaining access to capital markets
%%%%%%%%%%%%%%%%%% next period
lamda=0.3;

%%%%%%%%%%%%%%%%%% GDP loss during default (not varying with the level of
%%%%%%%%%%%%%%%%%% GDP.)
tau=0.2;
%tau=[0.1,0.4,0.5];

%__________________________________________________________________________ 
%% 1. Discretization of the state space
% Possible levesl of debt issuance $b\in B= \{0,0.05,0.1,...,0.5\}$

B=0:0.05:0.5;
%B=0:0.05:2.8;

%__________________________________________________________________________
%% 2. Initial values

%   We set the inital guess by assuming that GDP and the level of debt do not 
% change over time.
%   Notice that we only need initial guess for the Value functions with
% default and the before the default decision (Vd and V - not for Vnd)
%   Initial guess of bonds price is 0
 
%%%%%%%%%%%%%%%%%% Value function of no default 
V= ones(size(B,2),1)*u(y_grid)/(1-betta);

%%%%%%%%%%%%%%%%%% Value function of default
Vd = u((1-tau).*y_grid)/(1-betta);

%__________________________________________________________________________
%% 3. Equilibrium Computation

%%%%%%%%%%%%%%%%%% Initial guess for bond-price schedule
% assume that Gov. never defaults
q=ones(size(B,2),size(y_grid,2));

%%%%%%%%%%%%%%%%%% Iterations on the price schedule q(b',y)

% Iteration parameters:
max_iter=1000;
max_iterq=1000;
kapa = 0.3;

iter=1;
iter_q=1;

while iter_q<max_iterq
    
iter_q 
 
%%%%%%%%%%%%%%%%%% Value function iterations
while iter<max_iter
    
    for iy=1:size(y_grid,2)  
        for ib=1:size(B,2)
            
            Vnd_iter(ib,iy) = max(u(max(y_grid(iy)+q(:,iy).*B'-B(ib),0))+betta*V*piy(iy,:)');
            Vd_iter(iy)=u((1-tau)*y_grid(iy))+(1-lamda)*betta*piy(iy,:)*Vd'+lamda*betta*V(1,:)*piy(iy,:)';
            
            V_iter(ib,iy)=max(Vnd_iter(ib,iy),Vd_iter(iy));
        end   
    end
    
    dev = max(max(abs([V_iter-V;Vd_iter-Vd])));
    
        if dev<=0.000000000001
            break
        end
        
     Vd= Vd_iter;
     Vnd = Vnd_iter;
     V = V_iter;
        
    iter=iter+1;
        
end

%%%%%%%%%%%%%%%%%%  Updating the bond price menu
for iy=1:size(y_grid,2)
    for ib=1:size(B,2)
        
        % computing q for every possible combination of (b',y)
        q_iter(ib,iy)=1-piy(iy,:)*double(Vnd(ib,:)<=Vd)';
        
        % updatind q
        %q_up(ib,iy)=kapa.*q_iter(ib,iy)+(1-kapa).*q(ib,iy);
        
    end
end
 
  q=q_iter;
  %dev_q = max(max(abs(q_up-q)));
  %dev_q = max(max(abs(q_iter-q)));
 
  %if dev_q <=0.000000000001
      break
  %else
      q_iter = q;
  %end
  
  iter_q=iter_q+1;
  
end
fprintf('Done computing Bond-price menu, iteration: %3i, q distance: %2.6f \n',[iter_q]);

%% 4. Results

%%%%%%%%%%%%%%%%%% recovering the policy function

 for iy=1:size(y_grid,2)
       for ib=1:size(B,2)
           
       % No default value and debt issuance (conditional on no deault)
        [Vnd(ib,iy),bp(ib,iy)]=max(u(max(y_grid(iy)+q(:,iy).*B'-B(ib),0))+betta*V*piy(iy,:)')  ;
        
        % Policy function for defaulting: default decision
        gD(ib,iy)=double(Vnd(ib,iy)<=Vd(iy)) ; 
        
        %q_eq(ib,iy)=1-piy(iy,:)*gD(ib,:)'; %checking
       
    end
 end
 
% Policy function for bonds
gB=B(bp);

%%%%%%%%%%%%%%%%%% Graphs
figure(1)
subplot(2,2,1);
plot_q(B,q);
title('Bond-price menu')

figure(2)
subplot(1,2,1);
plot(B,gB)
title('Policy function for Bonds: if no default')

subplot(1,2,2);
plot(B,gD)
title('Policy function for Defaulting')

%__________________________________________________________________________ 

%% 5. Simulated sequence of GDP (exogenous)
T=500;

%  starting value (index)
yt=1;

for t=2:T
    draw_t=rand;
    yt(t)=1+(draw_t>=piy(yt(t-1),1))+(draw_t>=sum(piy(yt(t-1),1:2)));
end

%% 6. Initial values

%%%%%%%%%%%%%%%%%% Debt: index in B
bt=ones(T,1);

%%%%%%%%%%%%%%%%%% Default decision: index =1 and the default state
Def_b=nan(1,T);
Def_state(T)=0;

%% 7. Simulation

for t=2:T
    
%%%%%%%%%%%%%%%%%% Decisions Path
   
   if Def_state(t-1)==0 % had not defaulted
       
       Def_b(t)=gD(bt(t-1),yt(t)); % default decision (decided at t)     
       
       if Def_b(t)==0 % if decides not to default  
           
           bt(t)=bp(bt(t-1),yt(t)); % Debt issuance decision (decided at t) 
           Def_state(t)=0;
           
       else % it decides to default
           
           bt(t)=1; % no debt issuance
           Def_state(t)=1;
           
       end
   
   elseif rand<=lamda % if had defaulted & get to return to financial markets
   
       Def_b(t)=gD(bt(t-1),yt(t));
       
       if Def_b(t)==0 % if does not default again
           
           bt(t)=bp(bt(t-1),yt(t)); % Debt issuance decision (decided at t) 
           Def_state(t)=0;
           
       else % it deaults again
           
           bt(t)=1;
           Def_state(t)=1;
       
       end
   
   else % if had defaulted & do not get to return to financial markets
       
       Def_state(t)=1; % defaults again
       bt(t)=1;
       
   end
   
%%%%%%%%%%%%%%%%%%   Observed risk spread (1/q-1)

   r_spread(t)=1/q(bt(t),yt(t))-1;
 
%%%%%%%%%%%%%%%%%%    Default probability

   p_model(t)=1-q(bt(t),yt(t));
   
end

%%%%%%%%%%%%%%%%%% Graphs
figure(1)
subplot(1,2,1);
plot(Def_b);
title('Deafult decision path')

subplot(1,2,2);
plot(B(bt))
title('Debt issuance path')
ylim([0 1])

figure(2)
risk_premia_graph(Def_state, r_spread)
title('Risk preamia')

disp('Number of period in default');
disp(size(find(Def_b==1),2));
disp('Proportion on time in default:');
disp(size(find(Def_b==1),2)/T);

%__________________________________________________________________________

%% %%%%%%%%%%%%%%%%% Calibration %%%%%%%%%%%%%%%%%%%%%%
%   Adjust the parameters to have that the country is 3% of the time in
%   default
%   Results above show that the country never defaults.
%   As mentioned by the authors, to replicate levels of debt high enough to
%   induce default calibration usually relies on betta*R beign below 1 and
%   on state-contingent punishment, that is, tau varying with output level
%   realization
%   With this in mind, the following code is the same as the one above
%   execpt that now betta is lower, tau can take 3 values (same size as
%   the output grid), sigma is lower, lambda is higher
%__________________________________________________________________________ 

%% 0. Model Setup & Parameters 
clear;
clc;

%%%%%%%%%%%%%%%%%% Discount factor 
betta=0.9; %CHANGED

%%%%%%%%%%%%%%%%%% Possible values for GDP
%y_grid=[0.9, 1, 1.05];
y_grid=[0.6, 1, 1.5]; %CHANGED

%%%%%%%%%%%%%%%%%% Transition matrix for GDP
%piy=[0.5, 0.3, 0.2;...
%     0.1, 0.6, 0.3;...
%     0.2, 0.4, 0.4];

piy=[0.5, 0.3, 0.2;...
     0.05, 0.65, 0.3;...
     0.02, 0.55, 0.43];

%%%%%%%%%%%%%%%%%% Risk aversion
%sig=2;
sig=1.5;

%%%%%%%%%%%%%%%%%% Utility function $u(c)=\frac{c^{1-\sigma}}{1-\sigma}$
u = @(c) c.^(1-sig)./(1-sig);

%%%%%%%%%%%%%%%%%% Risk-free interest rate
R=1;

%%%%%%%%%%%%%%%%%% Probability of regaining access to capital markets
%%%%%%%%%%%%%%%%%% next period
lamda=0.5;

%%%%%%%%%%%%%%%%%% GDP loss during default (not varying with the level of
%%%%%%%%%%%%%%%%%% GDP.)
%tau=0.2;
tau=[0.1,0.4,0.5];

%__________________________________________________________________________ 
%% 1. Discretization of the state space
% Possible levesl of debt issuance $b\in B= \{0,0.05,0.1,...,0.5\}$

%B=0:0.05:0.5;
B=0:0.05:2.8;

%__________________________________________________________________________
%% 2. Initial values

%   We set the inital guess by assuming that GDP and the level of debt do not 
% change over time.
%   Notice that we only need initial guess for the Value functions with
% default and the before the default decision (Vd and V - not for Vnd)
%   Initial guess of bonds price is 0
 
%%%%%%%%%%%%%%%%%% Value function of no default 
V= ones(size(B,2),1)*u(y_grid)/(1-betta);

%%%%%%%%%%%%%%%%%% Value function of default
Vd = u((1-tau).*y_grid)/(1-betta);

%__________________________________________________________________________
%% 3. Equilibrium Computation

%%%%%%%%%%%%%%%%%% Initial guess for bond-price schedule
% assume that Gov. never defaults
q=ones(size(B,2),size(y_grid,2));

%%%%%%%%%%%%%%%%%% Iterations on the price schedule q(b',y)

% Iteration parameters:
max_iter=1000;
max_iterq=1000;
kapa = 0.3;

iter=1;
iter_q=1;

while iter_q<max_iterq
    
iter_q 
 
%%%%%%%%%%%%%%%%%% Value function iterations
while iter<max_iter
    
    for iy=1:size(y_grid,2)  
        for ib=1:size(B,2)
            
            Vnd_iter(ib,iy) = max(u(max(y_grid(iy)+q(:,iy).*B'-B(ib),0))+betta*V*piy(iy,:)');
            Vd_iter(iy)=u((1-tau(iy))*y_grid(iy))+(1-lamda)*betta*piy(iy,:)*Vd'+lamda*betta*V(1,:)*piy(iy,:)';
            
            V_iter(ib,iy)=max(Vnd_iter(ib,iy),Vd_iter(iy));
        end   
    end
    
    dev = max(max(abs([V_iter-V;Vd_iter-Vd])));
    
        if dev<=0.000000000001
            break
        end
        
     Vd= Vd_iter;
     Vnd = Vnd_iter;
     V = V_iter;
        
    iter=iter+1;
        
end

%%%%%%%%%%%%%%%%%%  Updating the bond price menu
for iy=1:size(y_grid,2)
    for ib=1:size(B,2)
        
        % computing q for every possible combination of (b',y)
        q_iter(ib,iy)=1-piy(iy,:)*double(Vnd(ib,:)<=Vd)';
        
        % updatind q
        %q_up(ib,iy)=kapa.*q_iter(ib,iy)+(1-kapa).*q(ib,iy);
        
    end
end
 
  q=q_iter;
  %dev_q = max(max(abs(q_up-q)));
  %dev_q = max(max(abs(q_iter-q)));
 
  %if dev_q <=0.000000000001
      break
  %else
      q_iter = q;
  %end
  
  iter_q=iter_q+1;
  
end
%fprintf('Done computing Bond-price menu, iteration: %3i, q distance: %2.6f \n',[iter_q]);

%% 4. Results

%%%%%%%%%%%%%%%%%% recovering the policy function

 for iy=1:size(y_grid,2)
       for ib=1:size(B,2)
           
       % No default value and debt issuance (conditional on no deault)
        [Vnd(ib,iy),bp(ib,iy)]=max(u(max(y_grid(iy)+q(:,iy).*B'-B(ib),0))+betta*V*piy(iy,:)')  ;
        
        % Policy function for defaulting: default decision
        gD(ib,iy)=double(Vnd(ib,iy)<=Vd(iy)) ; 
        
        %q_eq(ib,iy)=1-piy(iy,:)*gD(ib,:)'; %checking
       
    end
 end
 
% Policy function for bonds
gB=B(bp);

%%%%%%%%%%%%%%%%%% Graphs
figure(1)
plot_q(B,q);
title('Bond-price menu')

figure(2)
subplot(1,2,2);
plot(B,gB)
title('Policy function for Bonds: if no default')

subplot(2,2,3);
plot(B,gD)
title('Policy function for Defaulting')


%__________________________________________________________________________

%% 5. Simulated sequence of GDP (exogenous)
T=500;

%  starting value (index)
yt=1;

for t=2:T
    draw_t=rand;
    yt(t)=1+(draw_t>=piy(yt(t-1),1))+(draw_t>=sum(piy(yt(t-1),1:2)));
end

%% 6. Initial values

%%%%%%%%%%%%%%%%%% Debt: index in B
bt=ones(T,1);

%%%%%%%%%%%%%%%%%% Default decision: index =1 and the default state
Def_b=nan(1,T);
Def_state(T)=0;

%% 7. Simulation

for t=2:T
    
%%%%%%%%%%%%%%%%%% Decisions Path
   
   if Def_state(t-1)==0 % had not defaulted
       
       Def_b(t)=gD(bt(t-1),yt(t)); % default decision (decided at t)     
       
       if Def_b(t)==0 % if decides not to default  
           
           bt(t)=bp(bt(t-1),yt(t)); % Debt issuance decision (decided at t) 
           Def_state(t)=0;
           
       else % it decides to default
           
           bt(t)=1; % no debt issuance
           Def_state(t)=1;
           
       end
   
   elseif rand<=lamda % if had defaulted & get to return to financial markets
   
       Def_b(t)=gD(bt(t-1),yt(t));
       
       if Def_b(t)==0 % if does not default again
           
           bt(t)=bp(bt(t-1),yt(t)); % Debt issuance decision (decided at t) 
           Def_state(t)=0;
           
       else % it deaults again
           
           bt(t)=1;
           Def_state(t)=1;
       
       end
   
   else % if had defaulted & do not get to return to financial markets
       
       Def_state(t)=1; % defaults again
       bt(t)=1;
       
   end
   
%%%%%%%%%%%%%%%%%%   Observed risk spread (1/q-1)

   r_spread(t)=1/q(bt(t),yt(t))-1;
 
%%%%%%%%%%%%%%%%%%    Default probability

   p_model(t)=1-q(bt(t),yt(t));
   
end

%%%%%%%%%%%%%%%%%% Graphs
figure(1)
subplot(1,2,1);
plot(Def_b);
title('Deafult decision path')

subplot(1,2,2);
plot(B(bt))
title('Debt issuance path')
ylim([0 1])

figure(2)
risk_premia_graph(Def_state, r_spread)
title('Risk preamia')

disp('Number of period in default');
disp(size(find(Def_b==1),2));
disp('Proportion on time in default:');
disp(size(find(Def_b==1),2)/T);

% Debt-GDP ratio
ratio=B(bt)./yt;
plot(ratio,p_model);

%__________________________________________________________________________
%% %%%%%%%%%%%%%%%%% Data: Logit Estimation %%%%%%%%%%%%%%%%%%%%%%
%  Arranging the data. 
% 
%  We have to be sure that if a default spell lasts for more than one period, 
% we only include the first time default was declared (this is why we
% 
%  distinguish between the default decision and the default state)
% 
% 
% 
% % Number of observations 

N=sum(1-(isnan(Def_b)));
% Regressors. Constant and debt/GDP
X=[B(bt)'./y_grid(yt)'];

% Default decision
Y=Def_b';

% Logit regression. Binomial outcome (0, 1)
[par_est,dev,stats]=glmfit(X(1:end-1),Y(2:end),'binomial');


% Estimated Default probability
X_grid=0:0.1:2;
p_est=1./(1+exp(-par_est(1)-par_est(2)*X_grid));

Fit_model_graph(X(1:end-1), Y(2:end), X_grid, p_est)


% default probability in the simulation
p_est_sim=1./(1+exp(-par_est(1)-par_est(2)*X));

sim_and_model_graph([p_est_sim,p_model'])