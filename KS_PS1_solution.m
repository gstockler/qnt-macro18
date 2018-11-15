%% K&S (1998) - Quant Macro II - PS1

%% %%%%%%%%%%%%%%%%%%%%%% MODEL FRAMEWORK %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear;

%%%%%%%%%%%%%%%%%%  Transtition Matrices

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

Piz= [7/8  1/8; 1/8  7/8]; % good/bad

%__________________________________________________________________________

%%%%%%%%%%%%%%%%%% Parameters

betta=0.95;
delta=0.0025;
z=[1.01 0.99]; %good/bad
alfa=0.36;
L=[0.96, 0.9];

%%%%%%%%%%%%%%%%%% Capital grids

kgrid=[0:0.1:5,5.3:0.3:50]; % for individual k
Nk=size(kgrid,2); % number of points in kgrid

Kgrid=[16:0.04:18.5]; % for aggregate K
NK= size(Kgrid,2); % number of points in Kgrid

%__________________________________________________________________________

%%%%%%%%%%%%%%%%%% Parametrized expectations: initial/guessed values

b0g=0;
b1g=1;
b0b=0;
b1b=1;

%__________________________________________________________________________

%%%%%%%%%%%%%%%%%% Functions

% Utility and consumption
%   nK,nk: indicate which value on the grids
%   e: indicator for employed -> e=1,0
%   g: indicator for good z -> z=1==g,0==b

c = @(nk,nK,g,e) max( (alfa.*z(g).*(Kgrid(nK)/L(g)).^(alfa-1)+1-delta).*kgrid(nk) + ...
                (1-alfa).*z(g).*(Kgrid(nK)./L(g)).^(alfa).*e - kgrid,0 );
            
u = @(nk,nK,g,e) log( max( (alfa.*z(g).*(Kgrid(nK)/L(g)).^(alfa-1)+1-delta).*kgrid(nk) + ...
                (1-alfa).*z(g).*(Kgrid(nK)./L(g)).^(alfa).*e - kgrid,0 ) ); 
            
% VF guess

v1g = @(k,K) log( (alfa.*z(1).*(K/L(1)).^(alfa-1)-delta+1).*k + (1-alfa).*z(1).*(K./L(1)).^(alfa) )./(1-betta);
v1b = @(k,K) log( (alfa.*z(2).*(K./L(2)).^(alfa-1)-delta+1).*k + (1-alfa).*z(2).*(K./L(1)).^(alfa) )./(1-betta);

v0g = @(k,K) log( (alfa.*z(1).*(K./L(1)).^(alfa-1)-delta+1).*k )./(1-betta);
v0b = @(k,K) log( (alfa.*z(2).*(K./L(2)).^(alfa-1)-delta+1).*k ) ./(1-betta);

% Parametrized expectations:  g is an indicador(z=g): g=1 == zg

H=@(K,g) exp( ( b0g+b1g*log(K) ).*g + ( b0b+b1b*log(K) ).*(1-g) );

%__________________________________________________________________________

%%%%%%%%%%%%%%%%%% Initial guess evaluation for VFI

V1g = v1g(kgrid,Kgrid');
V1b = v1b(kgrid,Kgrid');
V0g = v0g(kgrid,Kgrid');
V0b = v0b(kgrid,Kgrid');

V1g=V1g';
V1b=V1b';
V0g=V0g';
V0b=V0b';

%__________________________________________________________________________

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%  LOOP  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% for "big"/complete loop
ITS_max=300; %1000
ITS = 0;

% for consumer's problem
maxits = 300; 
%tol = 0.00001;
dif_V = .5;
its=0;

% for simulation
N = 1000;
T = 2000;

for ITS=1:ITS_max
ITS
%%%%%%%%%%%%%%%%%%%%%%%%%%%  Functions
    
c = @(nk,nK,g,e) max( (alfa.*z(g).*(Kgrid(nK)/L(g)).^(alfa-1)+1-delta).*kgrid(nk) + ...
                (1-alfa).*z(g).*(Kgrid(nK)./L(g)).^(alfa).*e - kgrid,0 );
            
u = @(nk,nK,g,e) log( max( (alfa.*z(g).*(Kgrid(nK)/L(g)).^(alfa-1)+1-delta).*kgrid(nk) + ...
                (1-alfa).*z(g).*(Kgrid(nK)./L(g)).^(alfa).*e - kgrid,0 ) );  


H=@(K,g) exp( ( b0g+b1g*log(K) ).*g + ( b0b+b1b*log(K) ).*(1-g) );
%___________________________

%%%%%%%%%%%%%%%%%%%%%%%%%%%% CONSUMER
    while dif_V > 0.00001 & its<maxits
      
        for nk=1:Nk
            for nK=1:NK 
                [dif_K,i_K]=min( abs( Kgrid-H(Kgrid(nK),1) ) );
 
                V0g_its(nk,nK) = max( u(nk,nK,1,0)'+ betta .* ( Pize(1,:) * ([V0g(:,i_K),V1g(:,i_K),V0b(:,i_K),V1b(:,i_K)]'))');
                V1g_its(nk,nK)= max( u(nk,nK,1,1)'+ betta .* ( Pize(2,:) * ([V0g(:,i_K),V1g(:,i_K),V0b(:,i_K),V1b(:,i_K)]'))');  
                
                
                [dif_K,i_K]=min( abs( Kgrid-H(Kgrid(nK),0) ) );
 
                V0b_its(nk,nK) = max( u(nk,nK,2,0)'+ betta .* ( Pize(3,:) * ([V0g(:,i_K),V1g(:,i_K),V0b(:,i_K),V1b(:,i_K)]'))');
                V1b_its(nk,nK) = max( u(nk,nK,2,1)'+ betta .* ( Pize(4,:) * ([V0g(:,i_K),V1g(:,i_K),V0b(:,i_K),V1b(:,i_K)]'))');
            end
        end
        
         dif_V= max(max(abs( [V0g_its-V0g,V1g_its-V1g,V0b_its-V0b,V1b_its-V1b])));
         
         V0g=V0g_its;
         V1g=V1g_its;
         V0b=V0b_its;
         V1b=V1b_its;
         
         its=its+1;
    end
    
    % Policy function
    for nk=1:Nk
            for nK=1:NK 
                [dif_K,i_K]=min( abs( Kgrid-H(Kgrid(nK),1) ) );
 
                [V0g_its(nk,nK), a0g_its(nk,nK)] = max( u(nk,nK,1,0)'+ betta .* ( Pize(1,:) * ([V0g(:,i_K),V1g(:,i_K),V0b(:,i_K),V1b(:,i_K)]'))');
                [V1g_its(nk,nK), a1g_its(nk,nK)] = max( u(nk,nK,1,1)'+ betta .* ( Pize(2,:) * ([V0g(:,i_K),V1g(:,i_K),V0b(:,i_K),V1b(:,i_K)]'))');  
                
                
                [dif_K,i_K]=min( abs( Kgrid-H(Kgrid(nK),0) ) );
 
                [V0b_its(nk,nK), a0b_its(nk,nK)] = max( u(nk,nK,2,0)'+ betta .* ( Pize(3,:) * ([V0g(:,i_K),V1g(:,i_K),V0b(:,i_K),V1b(:,i_K)]'))');
                [V1b_its(nk,nK), a1b_its(nk,nK)] = max( u(nk,nK,2,1)'+ betta .* ( Pize(4,:) * ([V0g(:,i_K),V1g(:,i_K),V0b(:,i_K),V1b(:,i_K)]'))');
            end
        end
    
    % Store policy function as a(k,K,e,z)
    a = NaN(Nk,NK,2,2);
    a(:,:,2,1) = a0g_its;
    a(:,:,1,1)= a1g_its;
    a(:,:,2,2) = a0b_its;
    a(:,:,1,2)= a1b_its;

%___________________________

%%%%%%%%%%%%%%%%%%%%%%%%%%%%% SIMULATION
% Let the first simulated period define the employment status and shock
% received for each agent
% That is, for all periods, agents will have the same employment status and
% shock realization. Thus what is only changing is the asset distribution
% -> our focus

if ITS == 1
    
    % Shocks series
    zt(1)=1;  

    for t=2:T
        draw=rand;
        zt(t)= 1+(rand>=Piz(zt(t-1),1));
    end
    
    nzt=size(zt,2);

    zt_g=find(zt==1);
    nzt_g=size(zt_g,2);
    zt_b=find(zt==2);
    nzt_b=size(zt_b,2);
    
    if zt(1) ==1
        u0=.04*N;
    else
        u0=.1*N;
    end

    k0=find(kgrid == 17); 

    N_state(1:N-u0,:,1)=ones(N-u0,1)*[k0,1];
    N_state(N-u0+1:N,:,1)=ones(u0,1)*[k0,2];

    K_sim(1)=find(Kgrid == 17);

%%%%%%%%%%%%%%%%%% Simulated distributions

    for t=2:T
    for n=1:N
    
        N_state(n,1,t)=a(N_state(n,1,t-1),K_sim(t-1),N_state(n,2,t-1),zt(t-1));     
        N_state(n,2,t)= 2-(rand>=Pize(1 + zt(t-1)*2 - N_state(n,2,t-1),zt(t)*2-1)/Piz(zt(t-1),zt(t)));  
   
    end

    [dif_Ks, K_sim(t)]=min(abs(kgrid(round(mean(N_state(:,1,t))))-Kgrid));

    end
    
else
    
    for t=2:T
    for n=1:N
    
        N_state(n,1,t)=a( N_state(n,1,t-1),K_sim(t-1),N_state(n,2,t-1),zt(t-1) );
   
    end

    [dif_Ks, K_sim(t)]=min(abs(kgrid(round(mean(N_state(:,1,t))))-Kgrid));

    end

end
%___________________________

%%%%%%%%%%%%%%%%%%%%%%%%%%%%% REGRESSION: PARAMETRIZED EXP.

    burn = 20;
    ztr_g = zt_g(burn:end);
    ztr_b = zt_b(burn:end);

    % OLS for good times
    Yg=log(Kgrid(K_sim(ztr_g))');  
    Xg=[ones(size(ztr_g,2),1),log(Kgrid(K_sim(ztr_g-1))')] ;  
    
    Bg = Xg\Yg; 
    b0g_hat= Bg(1);
    b1g_hat = Bg(2);

    % OLS for bad times
    Yb=log(Kgrid(K_sim(ztr_b))');
    Xb=[ones(size(ztr_b,2),1),log(Kgrid(K_sim(ztr_b-1))')]  ;

    Bb = Xb\Yb;
    b0b_hat=Bb(1);
    b1b_hat=Bb(2);

    % Comparing new parameter estimates and the guess

    dif_b = max(abs([b0g-b0g_hat b1g-b1g_hat b0b-b0b_hat b1b-b1b_hat]));

    if dif_b<=0.01
        break 
    end

    % update guess
    b0g=0.1*b0g_hat+0.9*b0g;
    b1g=0.1*b1g_hat+0.9*b1g;
    b0b=0.1*b0b_hat+0.9*b0b;
    b1b=0.1*b1b_hat+0.9*b1b;
%___________________________

ITS=ITS+1;

end

%__________________________________________________________________________

%% %%%%%%% SOLVED MODEL ANALYSIS
disp('Iterations until convergence of Household problem (VFI)')
disp(its)

disp('Iterations to solve the model')
disp(ITS)

mean_sim = mean(Kgrid(K_sim));
disp('Simulated series Aggregate/Average Capital is')
disp(mean_sim)

for t=1:T
    if zt(t) == 1
        g=1;
    else
        g=0;
    end
    
    K_prime(t) = H(Kgrid(K_sim(t)),g);
end

mean_approx = mean(K_prime)
disp('Approximated Aggregate/Average Capital is')
disp(mean_approx)


%%%%%%%%%%%%%%%%%% Policy Function Plots
% graph is for a given K, e, z
% consider K=Kgrid(1)

% recovering policy function separately to plot
    % for unemployed at good state
gk_0g = kgrid(a(:,1,2,1)) ;
    % for employed at good state
gk_1g = kgrid(a(:,1,1,1)) ;

    % for unemployed at bad state
gk_0b = kgrid(a(:,1,2,2)) ;
     % for employed at bad state
gk_1b = kgrid(a(:,1,1,2)) ;

figure(1)
plot(kgrid, gk_0g )
hold on
plot(kgrid, gk_1g )
hold on
plot(kgrid,kgrid , 'y--')
title ('Asset policy function at good state')
legend('for unemployed', 'for employed', '45^0')
xlabel('k today ')
ylabel('k tomorrow')
print -dpdf ag_fig1.eps
hold off

figure(2)
plot(kgrid, gk_0b )
hold on
plot(kgrid, gk_1b )
plot(kgrid,kgrid , 'y--')
title ('Asset policy function at bad state')
legend('for unemployed', 'for employed','45^0')
xlabel('k today ')
ylabel('k tomorrow')
print -dpdf ab_fig2.eps
hold off

figure(3)
plot(kgrid, gk_0b )
hold on
plot(kgrid, gk_1b )
hold on
plot(kgrid, gk_0g )
hold on
plot(kgrid, gk_1g )
title ('Asset policy function')
legend('bad/unemployed', 'bad/employed', 'good/unemployed', 'good/employed')
xlabel('k today ')
ylabel('k tomorrow')
print -dpdf a_fig3.eps
hold off

%%%%%%%%%%%%%%%%%% Regression Analysis

mdlg = fitlm(Xg(:,2),Yg);
R2g = mdlg.Rsquared.Ordinary;
stdg= mdlg.RMSE; 

mdlb = fitlm(Xb(:,2),Yb);
R2b = mdlb.Rsquared.Ordinary;
stdb= mdlb.RMSE;

% For GOOD times: aggregate capital tomorrow
for t=1:nzt_g    
    K_pg(t) = H(Kgrid(K_sim(zt_g(t))),1);
end

% For BAD times: aggregate capital tomorrow
for t=1:nzt_b    
    K_pb(t) = H(Kgrid(K_sim(zt_b(t))),0);
end

figure(4)
plot(Kgrid(K_sim(zt_g)), K_pg);
hold on
plot(Kgrid(K_sim(zt_b)), K_pb);
hold on
plot(Kgrid(K_sim),Kgrid(K_sim));
title('Aggregate Capital today vs. tomorrow (Figure 1)' )
xlabel(' today' )
ylabel('tomorrow')
legend(strcat('at good times: R^2 =', num2str(R2g), 'MSE=', num2str(stdg)),strcat('at bad times: R^2 =',num2str(R2b), 'MSE=', num2str(stdb)) , '45^o')
print -dpdf aggk_fig4.eps

%%%%%%%%%%%%%%%%%% Capital/weath distribution

for t_ind=1:100  
    
    hist(kgrid(reshape(N_state(:,1,t_ind),1,1000)),40)
    title('Capital Distribution for 1000 individuals after 100 periods ')

end
print -dpdf histk_fig5.eps

%%%%%%%%%%%%%%%%%% Employment distribution

for t_ind=1:100  
    hist(reshape(N_state(:,2,t_ind),1,1000),40)
    title('Employment Distribution for 1000 individuals after 100 periods ')
    xticks([1 2])
    xticklabels({'Employed','Unemployed'})
end
print -dpdf histe_fig6.eps

%%%%%%%%%%%%%%%%%% Asset distribution after 7 periods in bad vs good states
% for consecutive BAD times 
grouped_b = mat2cell( zt_b, 1, diff( [0, find(diff(zt_b) ~= 1), length(zt_b)] )) ;
    
for i=1:size(grouped_b,2)
    if size(grouped_b{i},2)==7
     g7b = grouped_b{i};
     break
    end
end

for t=g7b(1):g7b(7)
hist(kgrid(reshape(N_state(:,1,t),1,1000)),40);
title('Capital Distribution for 1000 individuals after 7 consecutive bad periods ')
end
print -dpdf hist7b_fig7.eps

% for good 
grouped_g = mat2cell( zt_g, 1, diff( [0, find(diff(zt_g) ~= 1), length(zt_g)] )) ;
    
for i=1:size(grouped_g,2)
    if size(grouped_g{i},2)==7
     g7g = grouped_g{i};
     break
    end
end

for t=g7g(1):g7g(7)
hist(kgrid(reshape(N_state(:,1,t),1,1000)),40);
title('Capital Distribution for 1000 individuals after 7 consecutive bad periods ')
end
print -dpdf hist7g_fig8.eps




