%% K&S (1998) - Quant Macro II - PS1 
%
% %%%%%%%%%%%%%%%%%%%%%%%%%  SOLUTION II  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Now only 50% of agents update their PLM 
% The ones that do not update have coefficients equal to the initial guess
% always
% 
%__________________________________________________________________________

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

% NEW: for those who do not update -> notice that it is independent of z
% state now: 
%   H_const =@(K) exp( log(K) ) = K

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
% NEW
V1gc = v1g(kgrid,Kgrid')';
V1bc = v1b(kgrid,Kgrid')';
V0gc = v0g(kgrid,Kgrid')';
V0bc = v0b(kgrid,Kgrid')';



%__________________________________________________________________________

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%  LOOP  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tic; 
% for "big"/complete loop
ITS_max=300; %1000
ITS = 0;

% for consumer's problem
maxits = 300; 
%tol = 0.00001;
dif_V = .5;
dif_Vc = .5;
its=0;
itsc=0;

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

if ITS == 1
%%%%%%%%%%%%%%%%%%%%%%%%%%%% CONSUMER
    while dif_V > 0.00001 & dif_Vc > 0.00001 & its<maxits
      
        for nk=1:Nk
            for nK=1:NK 
                [dif_K,i_K]=min( abs( Kgrid-H(Kgrid(nK),1) ) );
 
                V0g_its(nk,nK) = max( u(nk,nK,1,0)'+ betta .* ( Pize(1,:) * ([V0g(:,i_K),V1g(:,i_K),V0b(:,i_K),V1b(:,i_K)]'))');
                V1g_its(nk,nK)= max( u(nk,nK,1,1)'+ betta .* ( Pize(2,:) * ([V0g(:,i_K),V1g(:,i_K),V0b(:,i_K),V1b(:,i_K)]'))');  
                
                
                [dif_K,i_K]=min( abs( Kgrid-H(Kgrid(nK),0) ) );
 
                V0b_its(nk,nK) = max( u(nk,nK,2,0)'+ betta .* ( Pize(3,:) * ([V0g(:,i_K),V1g(:,i_K),V0b(:,i_K),V1b(:,i_K)]'))');
                V1b_its(nk,nK) = max( u(nk,nK,2,1)'+ betta .* ( Pize(4,:) * ([V0g(:,i_K),V1g(:,i_K),V0b(:,i_K),V1b(:,i_K)]'))');
            
                % For those with constant beliefs
                V0g_itsc(nk,nK) = max( u(nk,nK,1,0)'+ betta .* ( Pize(1,:) * ([V0gc(:,nK),V1gc(:,nK),V0bc(:,nK),V1bc(:,nK)]'))');
                V1g_itsc(nk,nK)= max( u(nk,nK,1,1)'+ betta .* ( Pize(2,:) * ([V0gc(:,nK),V1gc(:,nK),V0bc(:,nK),V1bc(:,nK)]'))');  

                V0b_itsc(nk,nK) = max( u(nk,nK,2,0)'+ betta .* ( Pize(3,:) * ([V0gc(:,nK),V1gc(:,nK),V0bc(:,nK),V1bc(:,nK)]'))');
                V1b_itsc(nk,nK) = max( u(nk,nK,2,1)'+ betta .* ( Pize(4,:) * ([V0gc(:,nK),V1gc(:,nK),V0bc(:,nK),V1bc(:,nK)]'))');
            end
        end
        
         dif_V= max(max(abs( [V0g_its-V0g,V1g_its-V1g,V0b_its-V0b,V1b_its-V1b])));
         dif_Vc= max(max(abs( [V0g_itsc-V0gc,V1g_itsc-V1gc,V0b_itsc-V0bc,V1b_itsc-V1bc])));
         
         V0g=V0g_its;
         V1g=V1g_its;
         V0b=V0b_its;
         V1b=V1b_its;
         
         V0gc=V0g_itsc;
         V1gc=V1g_itsc;
         V0bc=V0b_itsc;
         V1bc=V1b_itsc;
         
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
            
                 % For those with constant beliefs
                [V0g_itsc(nk,nK),a0g_itsc(nk,nK)] = max( u(nk,nK,1,0)'+ betta .* ( Pize(1,:) * ([V0gc(:,nK),V1gc(:,nK),V0bc(:,nK),V1bc(:,nK)]'))');
                [V1g_itsc(nk,nK),a1g_itsc(nk,nK)] = max( u(nk,nK,1,1)'+ betta .* ( Pize(2,:) * ([V0gc(:,nK),V1gc(:,nK),V0bc(:,nK),V1bc(:,nK)]'))');  

                [V0b_itsc(nk,nK),a0b_itsc(nk,nK)] = max( u(nk,nK,2,0)'+ betta .* ( Pize(3,:) * ([V0gc(:,nK),V1gc(:,nK),V0bc(:,nK),V1bc(:,nK)]'))');
                [V1b_itsc(nk,nK),a1b_itsc(nk,nK)] = max( u(nk,nK,2,1)'+ betta .* ( Pize(4,:) * ([V0gc(:,nK),V1gc(:,nK),V0bc(:,nK),V1bc(:,nK)]'))');
           
            end
    end

    % Store policy function as a(k,K,e,z)
    a = NaN(Nk,NK,2,2);
    a(:,:,2,1) = a0g_its;
    a(:,:,1,1)= a1g_its;
    a(:,:,2,2) = a0b_its;
    a(:,:,1,2)= a1b_its;
    
    % NEW: For those who do not update beliefs
    ac = NaN(Nk,NK,2,2);
    ac(:,:,2,1) = a0g_itsc;
    ac(:,:,1,1)= a1g_itsc;
    ac(:,:,2,2) = a0b_itsc;
    ac(:,:,1,2)= a1b_itsc;
         
%%%%%%%%%%%%%%%%%%%%%%%%%%%%% SIMULATION     
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
% NEW: but how to make it random??
    for t=2:T
        
    for n=1:N
    
        N_state(n,1,t)=a(N_state(n,1,t-1),K_sim(t-1),N_state(n,2,t-1),zt(t-1));     
        N_state(n,2,t)= 2-(rand>=Pize(1 + zt(t-1)*2 - N_state(n,2,t-1),zt(t)*2-1)/Piz(zt(t-1),zt(t)));  
   
    end

    [dif_Ks, K_sim(t)]=min(abs(kgrid(round(mean(N_state(:,1,t))))-Kgrid));
    
    end
    
     % randomly picking those who do not update
    ntotal=linspace(1,N,N)';
    const=randperm(1000,500);
    up= ntotal(ismember(ntotal,const));
    
    for i=1:size(const,2)
        cc(i)=const(i);
        N_state_c(i,:,:) = N_state(cc(i),:,:);
    end
    
    for n=1:size(up)
        uu(n)=up(n);
        N_state_up(n,:,:)=N_state(uu(n),:,:);
     
    end
    
else % ONLY FOR THOSE WHO UPDATE: ITS>1
    
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
   for t=2:T
    %for i=1:size(const,2)
    %    %cc(i)=const(i);
    %    N_state(cc(i),1,t)= N_state_c(i,1,t);
  
    %end
    
 
    %for n=1:size(up,2)
    %    uu(n)=up(n);
    %    N_state(uu(n),1,t)=a( N_state(uu(n),1,t-1),K_sim(t-1),N_state(uu(n),2,t-1),zt(t-1) );
   
    %end
     for n=1:500
        N_state_up(n,1,t)=a( N_state_up(n,1,t-1),K_sim(t-1),N_state_up(n,2,t-1),zt(t-1) );
     
     end
    
      for i=1:size(const,2)
        cc(i)=const(i);
        N_state(cc(i),:,:) = N_state_c(i,:,:);
    end
    
    for n=1:size(up)
        uu(n)=up(n);
        N_state(uu(n),:,:)=N_state_up(n,:,:);
 
    end
    
    %[dif_Ks, K_sim(t)]=min(abs(kgrid(round(mean(N_state(:,1,t))))-Kgrid));
    [dif_Ks, K_sim(t)]=min(abs(kgrid(round(mean(N_state(:,1,t))))-Kgrid));

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
end

ITS=ITS+1;

end
toc;
%__________________________________________________________________________

%% %%%%%%%%%%%%%%%%%%%% SOLVED MODEL ANALYSIS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%% Comparing Welfare
% bad/employed
max(abs(max(V1b_its-V1b_itsc)));
% 0.0578

% good/employed
max(abs(max(V1g_its-V1g_itsc)));
% 0.0456

% bad/unemployed
max(abs(max(V0b_its-V0b_itsc)));
%0.0550
    
% bad/employed
max(abs(max(V0g_its-V0g_itsc)));
%0.0446

%%%%%%%%%%%%%%%%%% Policy Function Plots
% graph is for a given K, e, z
% consider K=Kgrid(1)

% Policy function for agents with best fit
    % for unemployed at good state
gk_0g = kgrid(a(:,1,2,1)) ;
    % for employed at good state
gk_1g = kgrid(a(:,1,1,1)) ;

    % for unemployed at bad state
gk_0b = kgrid(a(:,1,2,2)) ;
     % for employed at bad state
gk_1b = kgrid(a(:,1,1,2)) ;

% Policy function for agents with constant expectations
gk_0gc = kgrid(ac(:,1,2,1)) ;
    % for employed at good state
gk_1gc = kgrid(ac(:,1,1,1)) ;

    % for unemployed at bad state
gk_0bc = kgrid(ac(:,1,2,2)) ;
     % for employed at bad state
gk_1bc = kgrid(ac(:,1,1,2)) ;

figure(1)
subplot(2,2,1)
plot(kgrid, gk_0g,kgrid, gk_0gc )
title ('unemployed at good state')
legend('best fit', 'constant')

subplot(2,2,2)
plot(kgrid, gk_0b,kgrid, gk_0bc )
title ('unemployed at bad state')
legend('best fit', 'constant')

subplot(2,2,3)
plot(kgrid, gk_1g,kgrid, gk_1gc )
title ('employed at good state')
legend('best fit', 'constant')

subplot(2,2,4)
plot(kgrid, gk_1b,kgrid, gk_1bc )
title ('employed at bad state')
legend('best fit', 'constant')
suptitle('Asset policy function')
print -dpdf a2_fig1.eps

%%%%%%%%%%%%%%%%%% Asset Distribution
figure(2)
subplot(2,1,1)
for t_ind=1:100 

    hist(kgrid(reshape(N_state_c(:,1,t_ind),1,500)),40)
    title('Capital Distribution for agents with constant beliefs ')

end

subplot(2,1,2)
for t_ind=1:100 
  
    hist(kgrid(reshape(N_state_up(:,1,t_ind),1,500)),40)
    title('Capital Distribution for agents with best fit beliefs ')

end

print -dpdf histk_fig21.eps

