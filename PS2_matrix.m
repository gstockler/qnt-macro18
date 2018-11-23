%% K&S with Endogenous Labor Choice - Quant Macro II - PS2
% Gabriela Barbosa
%__________________________________________________________________________

%% %%%%%%%%%%%%%%%%%% Model without Aggregate Risk %%%%%%%%%%%%%%%%%%%%%%%%
clear;
clc;

%% 0. SET UP

%%%%%%%%%%%%%%%%%% Parameters

beta=0.95;
delta=0.0025;
z=[1.01 0.99]; %good/bad
e=[0;1];
alfa=0.36;
Gamma=8; % labor utility constant
gamma = 8; % labor utility parameter 
hp=.00001; % home production

%__________________________________________________________________________

%%%%%%%%%%%%%%%%%%  Transtition Matrices:

%%% Endogenous matrix of individual state (z,e) 

% The system of equations 
A= [ 1  0  0  0  1  0  0  0  0  0  0  0  0  0  0  0 ; ...
     0  1  0  0  0  1  0  0  0  0  0  0  0  0  0  0 ; ...
     0  0  0  0  0  0  0  0  0  0  1  0  0  0  1  0 ; ...
     0  0  0  0  0  0  0  0  0  0  0  1  0  0  0  1 ; ...
     0  0  0  0  0  0  0  0  1  0  0  0  1  0  0  0 ; ...
     0  0  0  0  0  0  0  0  0  1  0  0  0  1  0  0 ; ...
     0  0  1  0  0  0  1  0  0  0  0  0  0  0  0  0 ; ...
     0  0  0  1  0  0  0  1  0  0  0  0  0  0  0  0 ; ...
     1  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 ; ...
     0  0  0  0  0  0  0  0  0  0  1  0  0  0  0  0 ; ...
     0  0  0  0  0  0  0  0 5.6 0 -1  0  0  0  0  0 ; ...
    -1 0 28/3 0  0  0  0  0  0  0  0  0  0  0  0  0 ; ...
  .02 .48 .05 .45 0 0  0  0  0  0  0  0  0  0  0  0 ; ...
     0  1  0  0  0  0  0  0  0  0  0  0  0  0  0  0 ; ...
     0  0  0  0  0  0  0 0 .02 .48 .05 .45 0 0 0  0 ; ...
     0  0  0  0  0  0  0  0  0  0  0  1  0  0  0  0 ];
  
 
b= [7/8; 7/8; 7/8; 7/8; 1/8; 1/8; 1/8; 1/8; 7/24; 21/40; 0; 0; 0.02; 0.005; 0.05; 0.02];

Pize = reshape(A^-1*b,4,4); %good+unemp / good+emp / bad+unemp / bad+emp

Pize_g = Pize(1:2,1:2); % good: unemp/emp
for i=1:2
    Pize_gn(1,i)=Pize_g(1,i)./sum(Pize_g(1,:));
    Pize_gn(2,i)=Pize_g(2,i)./sum(Pize_g(2,:));
end

Pize_b = Pize(3:4,3:4); % bad: unemp/emp
for i=1:2
    Pize_bn(1,i)=Pize_b(1,i)./sum(Pize_b(1,:));
    Pize_bn(2,i)=Pize_b(2,i)./sum(Pize_b(2,:));
end

%%% Exogenous matrix of aggregate State z

Piz= [7/8  1/8; 1/8  7/8]; % good/bad

%__________________________________________________________________________

%%%%%%%%%%%%%%%%%% Discretezation: capital grid
ne=size(e,1);

%agrid=[0:0.1:5,5.3:0.3:50]; % for individual a
%na=size(agrid,2); % number of points in agrid

a_max=45;
a_min=0;

% Discretizing state space
na=250; %number of grids
agrid=zeros(na,1);
for i=1:na
    agrid(i)=a_min+(i-1).*(a_max-a_min)/(na-1);
    agrid(i,1)=agrid(i)';
end

agrid=agrid';
% vectorize the grid in two dimensions
amat = repmat(agrid',1,ne);            % values of a change vertically
emat = repmat(e',na,1);                % values of e change horizontally

%__________________________________________________________________________

%% 1. FUNCTIONS

%   We need write the problem of the HH in terms of only one unkown: asset
% level next period. 
%   To do so, we will use the functions below (more details in the pdf
%   document)

%%%%%%%%%%%%%%%%%%  (inverse) marginal utility functions

up    = @(c) 1./c;        % marginal utility of consumption
invup = @(x) x;      % inverse of marginal utility of consumption
vp    = @(n) Gamma.*n.^gamma;      % marginal disutility of labor supply
invvp = @(x) Gamma.^(-1/gamma) .*x.^(1./gamma);        % inverse marginal disutility of labor supply  

%%%%%%%%%%%%%%%%%%  endogenous functions

% optimal labor supply
F1  = @(c,s,w) invvp(up(c).*s.*w);

% current consumption level, cp0(anext,ynext) is the guess
F2 = @(cp0,r) invup(beta*(1+r)*up(cp0)*pi');
                
% current asset level, c0 = C0(cp0(anext,ynext))
F3 = @(anext,s,c0,r,w) 1/(1+r)...
                    *(c0+anext-F1(c0,s,w).*s.*w-hp);
%__________________________________________________________________________

%% 2A. Economy at BAD TIMES

%%%%%%%%%%%%%%%%%% Constant TFP
%   Let's fix TFP shock at bad state

zbar= z(1);
 pi=Pize_b; % picking bad times transtion matrix


%% 3A. Price Setting  

%   For a given value for the interest rate, the following gives the solution
% to the HH optimatization problem

r0=0.05; % guess: all results are conditional on it
R0=1+r0-delta;

%%%%%%%%%%%%%%%%%% Equilibrium prices and Aggregate Capital/Labor

%   Given the guess, we need to find the aggregate K to L ratio that
% satisfies firm's optimal decision. 
%   Then, we can compute equilibrium wage

f_KL = @(r) ((r-delta)./(zbar.*alfa)).^(1./(alfa-1));
KL=f_KL(r0);
w0= zbar.*(1-alfa).*KL.^alfa;

%% 4A. CONSUMPTION FUNCTION

% initial guess on future consumption (consume asset income plus labor
% income from working n=1 (e=1).

cp0 = r0*amat+emat*w0+hp;

c0 = F2(cp0,r0);

% derive current assets
a0 = F3(amat,emat,c0,r0,w0);

cp  = zeros(na,ne);
options = optimset('Display','off');

% for unemployed: e(1)
cp(:,1) = fsolve(@(c) (1+r0)*amat(:,1)+hp-c,cp0(:,1),options);
% for employed: e(2)
cp(:,2) = fsolve(@(c) (1+r0)*amat(:,2)+F1(c,e(2),w0).*e(2).*w0+hp-c,cp0(:,2),options);

% Recovering labor consistent with consumption function

cf=cp;

nf=zeros(na,ne);
for nk=1:na
    nf(nk,2)=F1(cf(nk,2),2,w0);
end

%n2=  Gamma.^(-1/gamma).*(e(2).*w0).^(1./gamma).*cf(:,2).^(-1./gamma);
%n1=zeros(na,1);
%n=[n1 n2];

%u=log(cf)-Gamma .* (1+gamma).^(-1).*n.^(1+gamma);

%% Checking consumption and labor functions
figure(1)
subplot(2,2,1)
plot(agrid, nf(:,2)) 
title('Labor vs Assets holdings');
hold off

subplot(2,2,2)
plot(agrid, cf(:,1))
hold on
plot(agrid, cf(:,2))
legend('for unemployed','for employed')
title('Consumption vs Assets holdings')
hold off;

subplot(2,2,3)
plot(cf(:,2),nf(:,2))
title('Consumption vs Labor for employed')
hold off;

subplot(2,2,4)
plot(nf(:,1),cf(:,1))
title('Consumption vs Labor for unemployed')
hold off;

suptitle('Consumption and Labor Endogenous Functions at bad times')
print -dpdf fig1.eps

%plot(agrid, [u(:,1), u(:,2)])

%% 5A. HOUSEHOLD PROBLEM

tol = 10^(-6);
dif_v    = tol+1;
maxits = 10^(3); 
its=0;

% initial guess for value function

V1b=zeros(na,1); % for employed
V0b=zeros(na,1); % for unemployed

%%%%%%%%%%%%%%%%%%%%%%%%%%%% CONSUMER FUNCTIONS

% Now we solve the problem taking into account the functions for
% consumption and labor; that is, the possible values for each variable at
% each state

cv = @(nk,s) max((1+r0).*agrid(nk)+w0.*nf(nk,s)-agrid+hp,0);
uv= @(nk,s) log(  max((1+r0).*agrid(nk)+w0.*nf(nk,s)-agrid+hp,0) ) -Gamma .* (1+gamma).^(-1).*nf(nk,s).^(1+gamma);

%%%%%%%%%%%%%%%%%%%%%%%%%%%% VFI

fprintf('VFI for bad times loop, running... \n');

while (dif_v>tol&&its<maxits)
      
    for nk=1:na
        
        V0b_its(nk,1) = max( uv(nk,1)'+ beta .* ( pi(1,:) * ( [V0b,V1b]' ) )' );
        V1b_its(nk,1) = max( uv(nk,2)'+ beta .* ( pi(2,:) * ( [V0b,V1b]' ) )' );
        
    end
        
    dif_v= max(max(abs( [V0b_its-V0b,V1b_its-V1b])));
         
    V0b=V0b_its;
    V1b=V1b_its;
         
    its=its+1;
    
end

fprintf('VFI loop for bad times done, iteration: %3i, VF distance: %2.6f \n',[its,dif_v]);

%% 6A. Results

%%%%%%%%%%%%%%%%%%%%%%%%%%%% Policy Functions

% Recovering asset policy function

for nk=1:na
    
    [V0b_its(nk,1), a0b_its(nk,1)] = max( uv(nk,1)+ beta .* ( pi(1,:) * ( [V0b,V1b]' ) ) );
    [V1b_its(nk,1), a1b_its(nk,1)] = max( uv(nk,2)+ beta .* ( pi(1,:) * ( [V0b,V1b]' ) ) );
    
end
       
ga_b = NaN(na,ne);
ga_b(:,1) = agrid(a0b_its);
ga_b(:,2)= agrid(a1b_its);

% Labor policy function

gn_b = NaN(na,ne);
gn(:,1) = nf(a0b_its,1);
gn_b(:,2)= nf(a1b_its,2);

% Consumption policy function

gc_b = NaN(na,ne);
gc_b(:,1) = cf(a0b_its,1);
gc_b(:,2)= cf(a1b_its,2);


%%%%%%%%%%%%%%%%%%%%%%%%%%%% Plots

figure(2)
subplot(2,2,1)
plot (agrid, V0b, 'b');
hold on
plot (agrid, V1b, 'r');
legend ('for unemployed', 'for employed');
title('Value Function');
hold off

subplot(2,2,2)
plot (agrid, ga_b(:,1), 'b');
hold on
plot (agrid, ga_b(:,2), 'r');
legend ('for unemployed', 'for employed');
title('Asset Policy Function');
hold off

subplot(2,2,3)
plot (agrid, gn_b(:,1), 'b');
hold on
plot (agrid, gn_b(:,2), 'r');
legend ('for unemployed', 'for employed');
title('Labor Policy Function');
hold off

subplot(2,2,4)
plot (agrid, gc_b(:,1), 'b');
hold on
plot (agrid, gc_b(:,2), 'r');
legend ('for unemployed', 'for employed');
title('Consumption Policy Function');
hold off
suptitle('Household Solution for bad times')

print -dpdf fig2.eps
%__________________________________________________________________________

%% 2B. Economy at GOOD TIMES

%%%%%%%%%%%%%%%%%% Constant TFP
%   Let's fix TFP shock at bad state

zbar= z(2);
pi=Pize_g; % picking bad times transtion matrix


%% 3B. Price Setting  

%   For a given value for the interest rate, the following gives the solution
% to the HH optimatization problem

r0=0.05; % guess: all results are conditional on it
R0=1+r0-delta;

%%%%%%%%%%%%%%%%%% Equilibrium prices and Aggregate Capital/Labor

%   Given the guess, we need to find the aggregate K to L ratio that
% satisfies firm's optimal decision. 
%   Then, we can compute equilibrium wage

f_KL = @(r) ((r-delta)./(zbar.*alfa)).^(1./(alfa-1));
KL=f_KL(r0);
w0= zbar.*(1-alfa).*KL.^alfa;

%% 4B. CONSUMPTION FUNCTION

% initial guess on future consumption (consume asset income plus labor
% income from working n=1 (e=1).

cp0 = r0*amat+emat*w0+hp;

% derive current consumption
    c0 = F2(cp0,r0);

% derive current assets
    a0 = F3(amat,emat,c0,r0,w0);

    cp  = zeros(na,ne);
    options = optimset('Display','off');

% for unemployed: e(1)
    cp(:,1) = fsolve(@(c) (1+r0)*amat(:,1)+hp-c,cp0(:,1),options);
% for employed: e(2)
    cp(:,2) = fsolve(@(c) (1+r0)*amat(:,2)+F1(c,e(2),w0).*e(2).*w0+hp-c,cp0(:,2),options);

% Recovering labor consistent with consumption function

cf=cp;

nf=zeros(na,ne);
for nk=1:na
    nf(nk,2)=F1(cf(nk,2),2,w0);
end


% Checking consumption and labor functions
figure(3)
subplot(2,2,1)
plot(agrid, n2) 
title('Labor vs Assets holdings');
hold off

subplot(2,2,2)
plot(agrid, cp(:,1))
hold on
plot(agrid, cp(:,2))
legend('for unemployed','for employed')
title('Consumption vs Assets holdings')
hold off;

subplot(2,2,3)
plot(nf(:,2),cp(:,2))
title('Consumption vs Labor for employed')
hold off;

subplot(2,2,4)
plot(nf(:,1),cp(:,1))
title('Consumption vs Labor for unemployed')
hold off;

suptitle('Consumption and Labor Endogenous Functions at good times')
print -dpdf fig3.eps

%% 5A. HOUSEHOLD PROBLEM

tol = 10^(-6);
dif_v    = tol+1;
maxits = 10^(3); 
its=0;

% initial guess for value function

V1g=zeros(na,1); % for employed
V0g=zeros(na,1); % for unemployed

%%%%%%%%%%%%%%%%%%%%%%%%%%%% CONSUMER FUNCTIONS

% Now we solve the problem taking into account the functions for
% consumption and labor; that is, the possible values for each variable at
% each state

cv = @(nk,s) max((1+r0).*agrid(nk)+w0.*nf(nk,s)-agrid+hp,0);
uv= @(nk,s) log(  max((1+r0).*agrid(nk)+w0.*nf(nk,s)-agrid+hp,0) ) -Gamma .* (1+gamma).^(-1).*nf(nk,s).^(1+gamma);

%%%%%%%%%%%%%%%%%%%%%%%%%%%% VFI

fprintf('VFI for good times loop, running... \n');

while (dif_v>tol&&its<maxits)
      
    for nk=1:na
        
        V0g_its(nk,1) = max( uv(nk,1)'+ beta .* ( pi(1,:) * ( [V0g,V1g]' ) )' );
        V1g_its(nk,1) = max( uv(nk,2)'+ beta .* ( pi(2,:) * ( [V0g,V1g]' ) )' );
        
    end
        
    dif_v= max(max(abs( [V0g_its-V0g,V1g_its-V1g])));
         
    V0g=V0g_its;
    V1g=V1g_its;
         
    its=its+1;
    
end

fprintf('VFI for good times loop done, iteration: %3i, VF distance: %2.6f \n',[its,dif_v]);

%% 6B. Results

%%%%%%%%%%%%%%%%%%%%%%%%%%%% Policy Functions

% Recovering asset policy function

for nk=1:na
    
    [V0g_its(nk,1), a0g_its(nk,1)] = max( uv(nk,1)+ beta .* ( pi(1,:) * ( [V0g,V1g]' ) ) );
    [V1g_its(nk,1), a1g_its(nk,1)] = max( uv(nk,2)+ beta .* ( pi(1,:) * ( [V0g,V1g]' ) ) );
    
end
       
ga_g = NaN(na,ne);
ga_g(:,1) = agrid(a0g_its);
ga_g(:,2)= agrid(a1g_its);

% Labor policy function

gn_g = NaN(na,ne);
gn_g(:,1) = nf(a0b_its,1);
gn_g(:,2)= nf(a1b_its,2);

% Consumption policy function

gc_g = NaN(na,ne);
gc_g(:,1) = cf(a0g_its,1);
gc_g(:,2)= cf(a1g_its,2);


%%%%%%%%%%%%%%%%%%%%%%%%%%%% Plots

figure(4)
subplot(2,2,1)
plot (agrid, V0g, 'b');
hold on
plot (agrid, V1g, 'r');
legend ('for unemployed', 'for employed');
title('Value Function');
hold off

subplot(2,2,2)
plot (agrid, ga_g(:,1), 'b');
hold on
plot (agrid, ga_g(:,2), 'r');
legend ('for unemployed', 'for employed');
title('Asset Policy Function');
hold off

subplot(2,2,3)
plot (agrid, gn_g(:,1), 'b');
hold on
plot (agrid, gn_g(:,2), 'r');
legend ('for unemployed', 'for employed');
title('Labor Policy Function');
hold off

subplot(2,2,4)
plot (agrid, gc_g(:,1), 'b');
hold on
plot (agrid, gc_g(:,2), 'r');
legend ('for unemployed', 'for employed');
title('Consumption Policy Function');
hold off
suptitle('Household Solution for good times')
print -dpdf fig4.eps
%__________________________________________________________________________
 
%% 7. GENERAL EQUILIBRIUM

% bounds for the interest-rate (which will solve for)
r_min=delta;
r_max=1/beta-1;

% initial guess for interest rate
rng(0,'twister');
r_guess = (r_max-r_min).*rand(1,1) + r_min;
r=r_guess;

% for the big loop
ED=1; %excess demand for capital
tolg=10^2;
maxitg=500;
itg=1;
lambda=.5; %weight for the convex combination of current r and next-step r

% for HH inside loop (same as above)
tol = 10^(-6);
dif_v    = tol+1;
maxits = 10^(3); 
its=1;

% LOOP
while (abs(ED)>tolg&&itg<maxitg)
    its=1;
    r=r_guess;
    
    % Firms
    KL=f_KL(r0);
    w= zbar.*(1-alfa).*KL.^alfa;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%% HOUSEHOLD
    
    V1=zeros(na,1); % for employed
    V0=zeros(na,1); % for unemployed

    cv = @(nk,s) max((1+r0).*agrid(nk)+w0.*nf(nk,s)-agrid+hp,0);
    uv= @(nk,s) log(  max((1+r0).*agrid(nk)+w0.*nf(nk,s)-agrid+hp,0) ) -Gamma .* (1+gamma).^(-1).*nf(nk,s).^(1+gamma);

    fprintf('VFI for good times loop, running... \n');

    while (dif_v>tol&&its<maxits)
      
        for nk=1:na
        
            V0_its(nk,1) = max( uv(nk,1)'+ beta .* ( pi(1,:) * ( [V0,V1]' ) )' );
            V1_its(nk,1) = max( uv(nk,2)'+ beta .* ( pi(2,:) * ( [V0,V1]' ) )' );
        
        end
        
    dif_v= max(max(abs( [V0_its-V0,V1_its-V1])));
         
    V0=V0_its;
    V1=V1_its;
         
    its=its+1;
    
    end

    fprintf('VFI for good times loop done, iteration: %3i, VF distance: %2.6f \n',[its,dif_v]);
    
    % Recovering asset policy function

    for nk=1:na
    
        [V0_its(nk,1), a0_its(nk,1)] = max( uv(nk,1)+ beta .* ( pi(1,:) * ( [V0,V1]' ) ) );
        [V1_its(nk,1), a1_its(nk,1)] = max( uv(nk,2)+ beta .* ( pi(1,:) * ( [V0,V1]' ) ) );
    
    end
    
     a_its=[a0_its'; a1_its'];
     V_its=[V0_its' V1_its'];
    
    ga = NaN(na,ne);
    ga(:,1) = agrid(a0_its);
    ga(:,2)= agrid(a1_its);

    gn = NaN(na,ne);
    gn(:,1) = nf(a0_its,1);
    gn(:,2)= nf(a1_its,2);

    gc = NaN(na,ne);
    gc(:,1) = cf(a0_its,1);
    gc(:,2)= cf(a1_its,2);

 %%%%%%%%%%%%%%%%%%%%%%%%%%%% INVARIANT DISTRIBUTION
 
    % Transition matrix (a,e) pairs: this is an endogenous process
    gt_1=zeros(na,na); % indicator functions
    gt_2=zeros(na,na);
    
    for i=1:na
        gt_1(i,a0_its(i))=1; %if e=0 and a=ai, then gt_1=1 is the best response -> ga1
        gt_2(i,a1_its(i))=1;  %if e=1 and a=ai, then gt_2=1 is the best response -> ga2
    end
    
    % 1-1 = Pi(1,1)*gt_1
    % 1-2 = Pi(1,2)*gt_1
    % 2-1 = Pi(2,1)*gt_2
    % 2-2 = Pi(2,2)*gt_2
    Trans=[ pi(1,1)*gt_1 pi(1,2)*gt_1; pi(2,1)*gt_2 pi(2,2)*gt_2];
    
    % Strationary Distribution: iteration until difference is very small
    Trans= Trans';
    probst = (1/(2*na))*ones(2*na,1);
    test = 1;
    error_trans=10^-7;

    while (test > error_trans)
        
        probst1 = Trans*probst;
     
       test=max(max(abs(probst1-probst)));
       
       probst = probst1;
    end;
    
    Stat_dist = probst;
    
    % AGG. SUPPLY OF CAPITAL: Assets
    aa=ga(:);
    meanA=Stat_dist'*aa;
    
    nn=gn(:);
    meanN=Stat_dist'*nn(:);
    
    % Measure of (a,e) pairs
    phi=zeros(na,ne);
    phi(:)=Stat_dist; %each column is the stat.dist of assets for a given y
   
    Dist_a= sum(phi');     %  stationary distribution of K: sums up to 1
    %row-vector: each cell is the prob. of having that level of K
    %(accounting for both poor and rich)
    
    Dist_a=Dist_a';

    % EXCESS DEMAND OF CAPITAL/LABOR RATIO (Demand - Supply)
    ED=KL-meanA/meanN;
    
    % Updating guess:
    r_old=r;
    r_new=alfa.*(meanA./meanN).^(alfa-1)-delta;
    r_guess=lambda*r_old+(1-lambda)*r_new;

    %if PHI<0
    %    r_guess=lambda*r_old+(1-lambda)*r_min;
    %else
    %   r_guess=lambda*r_old+(1-lambda)*r_max; 
    %end
    
    r=r_guess;
    
    iter=iter+1;
end

disp('Equilibrium rental rate:')
disp(r)
disp('Iterations until GE')
disp(itg);

%% Computing GE
r_ge=r;

% HHs optmize
V_ge=V_its;
ga1_ge=ga(:,1);
ga2_ge=ga(:,2);
gc1_ge=gc(:,1);
gc2_ge=gc(:,2);

% Firms optimize
KL_ge=f_KL(r_ge);
w_ge= zbar.*(1-alfa).*KL_ge.^alfa;

 %% Plots: endogenous distributions

% Capital
figure(1)
plot(agrid,Dist_a,repmat(K_ge,na,1),Dist_k);
xlabel('a')
ylabel('Probability mass')
legend('Distribution','Equilibrium Agg. Capital')
title('Stationary Distribution')

figure(2)
plot(agrid, phi(:,1))
hold on
plot(agrid,phi(:,2))
hold off
xlabel('a')
ylabel('Probability mass')
legend('conditional on y1','conditional on y2')
title('Stationary Distribution of Capital')

figure(3)
plot(agrid, cumsum(Dist_k),agrid,ones(na,1));
xlabel('a')
ylabel('Cumulative probability')
title('Stationary Cumulative Distribution of Capital')

% Wealth
wealth = [ (r_ge*a + y(1,1)*w_ge)  (r_ge*agrid + y(2,1)*w_ge) ]  ; 
[ Dist_w, ind_w ] = sort(wealth(:)); %ascend-ordering income values while keep track of location
phi_vec = phi(:); % transforming phi(na,ny) into vector(na*ny) - to match Dist_inc size

figure(4)
plot(Dist_w,phi_vec(ind_w));
xlabel('wealth')
ylabel('Probability mass')
title('Stationary Distribution of Wealth')

% Income
income = [ repmat(y(1,1)*w_ge,na,1),  repmat(y(2,1)*w_ge,na,1) ]  ; 
[ Dist_inc, ind_inc ] = sort(income(:)); %ascend-ordering income values while keep track of location

figure(4)
plot(Dist_inc(1:100,1),phi_vec(ind_inc(1:100,1)));
hold on
plot(Dist_inc(101:200,1),phi_vec(ind_inc(101:200,1)));
hold off
xlabel('income')
ylabel('Probability mass')
legend('for y1','for y2')
title('Stationary Distribution of Income') 

% Consumption
consumption = [wealth(:,1)+ga1, wealth(:,1)+ga2 ];
%c2=[gc1_ge, gc2_ge]; ?
[Dist_c, ind_c ] = sort(consumption(:));

figure(5)
plot(Dist_c,phi_vec(ind_c));
xlabel('consumption')
ylabel('Probability mass')
title('Stationary Distribution of Consumption')

figure(6)
plot(Dist_c(1:100,1),phi_vec(ind_c(1:100,1)));
hold on
plot(Dist_c(101:200,1),phi_vec(ind_c(101:200,1)));
hold off
xlabel('consumption')
ylabel('Probability mass')
legend('for y1','for y2')
title('Stationary Distribution of Consumption') 

% Demand and Supply of capital
KL_s=ga1_ge+ga2_ge; %supply schedule
r_grid1=linspace(r_min,r_ge,na/2)';
r_grid2=linspace(r_ge,r_max,na/2)';
r_grid=[r_grid1;r_grid2];
KL_d=f_KL(r_grid);
KL_d=KL_d'; %demand schedule

%figure(7)
plot(KL_d,r_grid)
hold on
plot(KL_s,r_grid)
hold on
plot(repmat(KL_ge,na,1),r_grid);
plot(KL_d,repmat(r_ge,na,1))
hold off

% _______________________________________________________________________ %   
