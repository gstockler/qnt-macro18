%% TRANSITION
clear all; clc;

T=100;
%betta=.9805;
thetta=.67;
delta=0.0625;
h=.31;

%% SS
[hs1,betta1,zs1,ks1,cs1,invs1,outts1] = steady_state(thetta,delta);
[hs2,betta2,zs2,ks2,cs2,invs2,outts2] = steady_state2(thetta,delta);
savss1=outts1-invs1;
savss2=outts2-invs2;
%% Persistent shock
% Storage matrices
kmat=zeros(T,1);
cmat=zeros(T,1);
ymat=zeros(T,1);
imat=zeros(T,1);

kmat(1,1)=ks1;
cmat(1,1)=cs1;
ymat(1,1)=outts1;
imat(1,1)=invs1;
kmat(T,1)=ks2;
cmat(T,1)=cs2;
ymat(T,1)=outts2;
imat(T,1)=invs2;

for i=1:T-2
    kmat(i+1,1)=ymat(i,1)-cmat(i,1)+(1-delta)*kmat(i,1);
    ymat(i+1,1)=kmat(i+1,1)^(1-thetta)*(zs2*h)^thetta;
    cmat(i+1,1)=.75*ymat(i+1,1);
    imat(i+1,1)=.25*ymat(i+1,1);
end

figure(1)
%plotting consumption path
subplot(2,3,1)
plot(cmat,'b');
hold on
plot(repmat(cs1,100,1),'r-.');
hold on
plot(repmat(cs2,100,1),'m-.');
ylim([0,2]);
title('consumption path');

subplot(2,3,2)
plot(kmat,'b');
hold on
plot(repmat(ks1,100,1),'r-.');
hold on
plot(repmat(ks2,100,1),'m-.');
ylim([3,9]);
title('capital path');

subplot(2,3,3)
plot(imat,'b');
hold on
plot(repmat(invs1,100,1),'r-.');
hold on
plot(repmat(invs2,100,1),'m-.');
ylim([0,1]);
title('investment path');

subplot(2,3,5)
plot(ymat,'b');
hold on
plot(repmat(outts1,100,1),'r-.');
hold on
plot(repmat(outts2,100,1),'m-.');
ylim([0,3]);
title('output path');

subplot(2,3,4)
plot((ymat-imat),'b');
hold on
plot(repmat(savss1,100,1),'r-.');
hold on
plot(repmat(savss2,100,1),'m-.');
ylim([0,2]);
title('savings path');



print -dpdf q1c_pers.eps


%% Temporary shock
kmat=zeros(T,1);
cmat=zeros(T,1);
ymat=zeros(T,1);
imat=zeros(T,1);

kmat(1,1)=ks1;
cmat(1,1)=cs1;
ymat(1,1)=outts1;
imat(1,1)=invs1;
kmat(T,1)=ks1;
cmat(T,1)=cs1;
ymat(T,1)=outts1;
imat(T,1)=invs1;


for i=1:T-2
    kmat(i+1,1)=ymat(i,1)-cmat(i,1)+(1-delta)*kmat(i,1);
    ymat(i+1,1)=kmat(i+1,1)^(1-thetta)*(zs2*h)^thetta;
    cmat(i+1,1)=.75*ymat(i+1,1);
    imat(i+1,1)=.25*ymat(i+1,1);
    if i>10
        kmat(i+1,1)=ymat(i,1)-cmat(i,1)+(1-delta)*kmat(i,1);
        ymat(i+1,1)=kmat(i+1,1)^(1-thetta)*(zs1*h)^thetta;
        cmat(i+1,1)=.75*ymat(i+1,1);
        imat(i+1,1)=.25*ymat(i+1,1);
    end
    
end

figure(2)
%plotting consumption path
subplot(2,3,1)
plot(cmat,'b');
hold on
plot(repmat(cs1,100,1),'r-.');
hold on
plot(repmat(cs2,100,1),'m-.');
ylim([0,2]);
title('consumption path');

subplot(2,3,2)
plot(kmat,'b');
hold on
plot(repmat(ks1,100,1),'r-.');
hold on
plot(repmat(ks2,100,1),'m-.');
ylim([2,9]);
title('capital path');

subplot(2,3,3)
plot(imat,'b');
hold on
plot(repmat(invs1,100,1),'r-.');
hold on
plot(repmat(invs2,100,1),'m-.');
ylim([0,1]);
title('investment path');

subplot(2,3,5)
plot(ymat,'b');
hold on
plot(repmat(outts1,100,1),'r-.');
hold on
plot(repmat(outts2,100,1),'m-.');
ylim([0,3]);
title('output path');

subplot(2,3,4)
plot((ymat-imat),'b');
hold on
plot(repmat(savss1,100,1),'r-.');
hold on
plot(repmat(savss2,100,1),'m-.');
ylim([0,2]);
title('savings path');

print -dpdf q1c_temp.eps