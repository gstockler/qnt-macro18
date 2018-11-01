%% PART II.5 - General Equilibrium with Liquidity Constraint
% The following is using discrete VFI and using PE results from II.2
% Quadratic Utility

%% Setting up the model
clear
close all;

% parameters of the model
rho =.06;
beta = 1/(1+rho);
cbar = 100;
sigma =2;
sigma_y = .1; 
gamma =0;
delta=.02;
theta=.36; %capital income share

% Markov Process (need to make it more general)
y=[1-sigma_y;1+sigma_y];
Pi=[(1+gamma)/2,(1-gamma)/2; (1-gamma)/2,(1+gamma)/2];
ny=length(y);

%% functional forms
%u1  = @(c) c.^(1-sigma)/(1-sigma)-1/(1-sigma);   % preferences
u2  = @(c) .5.*(c-cbar).^2 ;  % preferences

% bounds for the interest-rate (which will solve for)
r_min=delta;
r_max=rho;

% firms CRS production (standard)
% L=1 (exogenous)
f= @(K,L) K.^(theta).*L.^(1-theta);


%% Aggregate Labor and EXOGENOUS STAT. DISTRIBUTION
% we can back-out agg labor from the stationary distribution of shocks
% (Markov chain stationary distribution)

[eig_vec,eig_val]=eig(Pi'); % getting eigenvectors and values of transition matrix

% choose the the eigenvector corresponding to the column where the
% eigenvalue=1 is located
[e1,e1_loc]=max(diag(eig_vec));
Pi_st=eig_vec(:,e1_loc)/sum(eig_vec(:,e1_loc)); % stationary dist. (need to normalize so it sums up to 1)

% then, aggregate labor is given by sum_i(yi*Pi_st(i))
L=y'*Pi_st;

%% Solving HHs' problem given r,w above
% State-space bounds
a_max=30;
a_min=0;
%[a_min]=borrowing_const(2,y,r); want to make it outside the loop

% Discretizing state space
na=100; %number of grids
agrid=zeros(na,1);
for i=1:na
    agrid(i)=a_min+(i-1).*(a_max-a_min)/(na-1);
    agrid(i,1)=agrid(i)';
end

% converting grid into matrices: amat varies along columns, apmat varies along rows.
[amat,apmat] = ndgrid(agrid,agrid);     % ap denotes a'

% _______________________________________________________________________ %


%% GE-Loop
% initial guess for interest rate
rng(0,'twister');
r_guess = (r_max-r_min).*rand(1,1) + r_min;

%for the big loop
ED=1; %excess demand for capital
maxiter=500;
iter=0;
lambda=.8; %weight for the convex combination of current r and next-step r

% for HH inside loop
maxits = 400;
tol = 0.001;
dif = .5;
its=0;

% LOOP
while abs(ED)>.001 & iter<maxiter
    r=r_guess;
    
    % Firms
    K=((r+delta)./theta).^(1/(theta-1));
    w=(1-theta).*K.^theta.*L.^(-theta);

    % HH
    V1=zeros(na,1);
    V2=zeros(na,1);

    cmat1 = w*y(1,1)+(1+r)*amat-apmat;
    cmat1(cmat1<0)=[-99999999];
    cmat2 = w*y(2,1)+(1+r)*amat-apmat;
    cmat2(cmat2<0)=[-99999999];
    
    %U1mat_1=u1(cmat1);
    %U1mat_2=u1(cmat2);
    U2mat_1=u2(cmat1);
    U2mat_2=u2(cmat2);

    V= zeros(na,ny);
    while dif > tol & its<maxits

        [V1_its,a1_its]=max(U2mat_1 + beta*repmat(V*Pi(1,:)',1,na));
        [V2_its,a2_its]=max(U2mat_2 + beta*repmat(V*Pi(2,:)',1,na));
  
        a_its=[a1_its' a2_its'];
        V_its=[V1_its' V2_its'];
  
        ga=agrid(a_its);
        ga1=ga(:,1);
        ga2=ga(:,2);
        
        dif=max(max(abs(V_its-V)));
        V=V_its;
        its = its+1;
        
    end;

    a=agrid;
    gc1 = w*y(1,1)+(1+r)*a-ga1;
    gc2 = w*y(2,1)+(1+r)*a-ga2;

    % TRANSITION MATRIX (a,y) pairs: this is an endogenous process
    gt_1=zeros(na,na); % indicator functions
    gt_2=zeros(na,na);
    
    for i=1:na
        gt_1(i,a1_its(i))=1; %if y=y1 and a=ai, then gt_1=1 is the best response -> ga1
        gt_2(i,a2_its(i))=1;  %if y=y2 and a=ai, then gt_2=1 is the best response -> ga2
    end
    
    % 1-1 = Pi(1,1)*gt_1
    % 1-2 = Pi(1,2)*gt_1
    % 2-1 = Pi(2,1)*gt_2
    % 2-2 = Pi(2,2)*gt_2
    Trans=[ Pi(1,1)*gt_1 Pi(1,2)*gt_1; Pi(2,1)*gt_2 Pi(2,2)*gt_2];
    
    % STATIONARY DISTRIBUTION: iteration until difference is very small
    Trans= Trans';
    probst = (1/(2*na))*ones(2*na,1);
    test = 1;
    error_trans=.00000001;

    while test > error_trans
       probst1 = Trans*probst;
       test=max(abs(probst1-probst));
       probst = probst1;
    end;
    Stat_dist = probst;
    
    % AGG. SUPPLY OF CAPITAL: Assets
    aa=ga(:);
    meanA=Stat_dist'*aa;

    %[VV,DD] = eig(Trans);
    %dd= diag(DD); % check how many unit eigenvalues there are!!
                    % if kmin = 0, there are 2 (one trivial, only relevant if initial capital is zero)
                    % if kmin > 0, there is  1 
    
    %[smallest,index] = min(abs(dd-1)); %% find a unit eigenvalue
    %vv      = VV(:,index);

    %Qbar  = (vv/sum(vv))';
    
    % Measure of (a,y) pairs
    phi=zeros(na,ny);
    phi(:)=Stat_dist; %each column is the stat.dist of assets for a given y
   
    Dist_k= sum(phi');     %  stationary distribution of K: sums up to 1
    %row-vector: each cell is the prob. of having that level of K
    %(accounting for both poor and rich)
    
    Dist_k=Dist_k';

    % EXCESS DEMAND OF CAPITAL (Demand - Supply)
    ED=K-meanA;
    
    % Updating guess:
    r_old=r;
    r_new=theta.*meanA.^(theta-1)-delta;
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
disp(iter);

%% Computing GE
r_ge=r;

% HHs optmize
V_ge=V_its;
ga1_ge=ga1;
ga2_ge=ga2;
gc1_ge=gc1;
gc2_ge=gc2;

% Firms optimize
L_ge=L;
K_ge=((r_ge+delta)./theta).^(1/(theta-1));
w_ge=(1-theta).*K_ge.^theta.*L.^(-theta);

%% Plots: endogenous distributions

% Capital
figure(1)
plot(a,Dist_k,repmat(K_ge,na,1),Dist_k);
xlabel('a')
ylabel('Probability mass')
legend('Distribution','Equilibrium Agg. Capital')
title('Stationary Distribution')
%print -dpdf ge_1.eps

figure(2)
plot(a, phi(:,1))
hold on
plot(a,phi(:,2))
hold off
xlabel('a')
ylabel('Probability mass')
legend('conditional on y1','conditional on y2')
title('Stationary Distribution of Capital')
%print -dpdf ge_2.eps

figure(3)
plot(a, cumsum(Dist_k),a,ones(na,1));
xlabel('a')
ylabel('Cumulative probability')
title('Stationary Cumulative Distribution of Capital')
%print -dpdf ge_3.eps

% Wealth
wealth = [ (r_ge*a + y(1,1)*w_ge)  (r_ge*agrid + y(2,1)*w_ge) ]  ; 
[ Dist_w, ind_w ] = sort(wealth(:)); %ascend-ordering income values while keep track of location
phi_vec = phi(:); % transforming phi(na,ny) into vector(na*ny) - to match Dist_inc size

figure(4)
plot(Dist_w,phi_vec(ind_w));
xlabel('wealth')
ylabel('Probability mass')
title('Stationary Distribution of Wealth')
%print -dpdf ge_4.eps

% Income
income = [ repmat(y(1,1)*w_ge,na,1),  repmat(y(2,1)*w_ge,na,1) ]  ; 
[ Dist_inc, ind_inc ] = sort(income(:)); %ascend-ordering income values while keep track of location

figure(5)
plot(Dist_inc(1:100,1),phi_vec(ind_inc(1:100,1)));
hold on
plot(Dist_inc(101:200,1),phi_vec(ind_inc(101:200,1)));
hold off
xlabel('income')
ylabel('Probability mass')
legend('for y1','for y2')
title('Stationary Distribution of Income') 
%print -dpdf ge_5.eps

% Consumption
consumption = [wealth(:,1)+ga1, wealth(:,1)+ga2 ];
%c2=[gc1_ge, gc2_ge]; ?
[Dist_c, ind_c ] = sort(consumption(:));

figure(6)
plot(Dist_c,phi_vec(ind_c));
xlabel('consumption')
ylabel('Probability mass')
title('Stationary Distribution of Consumption')
%print -dpdf ge_6.eps

figure(6)
plot(Dist_c(1:100,1),phi_vec(ind_c(1:100,1)));
hold on
plot(Dist_c(101:200,1),phi_vec(ind_c(101:200,1)));
hold off
xlabel('consumption')
ylabel('Probability mass')
legend('for y1','for y2')
title('Stationary Distribution of Consumption') 
%print -dpdf ge_7.eps

% Demand and Supply of capital
%K_s=ga1_ge+ga2_ge; %supply schedule
%r_grid1=linspace(r_min,r_ge,na/2)';
%r_grid2=linspace(r_ge,r_max,na/2)';
%r_grid=[r_grid1;r_grid2];
%K_d=((r_grid+delta)./theta).^(1/(theta-1));
%K_d=K_d'; %demand schedule

%figure(7)
%plot(K_d,r_grid)
%hold on
%plot(K_s,r_grid)
%hold on
%plot(repmat(K_ge,na,1),r_grid);
%plot(K_d,repmat(r_ge,na,1))
%hold off

% _______________________________________________________________________ %