%% Multivariate CES Function Approximation - sigma=1
clc
clear all;

syms k h alpha sigma
f=((1-alpha)*k^((sigma-1)/sigma)+alpha*h^((sigma-1)/sigma))^(sigma/(sigma-1));
alpha=0.5;
sigma=0.99;
f_e=eval(f);

%% A) Show that sigma is the elasticity of substitution (ES)
n=diff(f_e,k)/diff(f_e,h);
d=k/h;
es=diff(log(n))/diff(log(d));
ES=-(1/es);
if sigma == ES;
    disp(ES,'Sigma is the ES')
end

%% B) Labor share
%It's defined as the marginal productivity of labor
syms h k
f_h=diff(f,h);
h_share=(h*f_h)/f;
he_share=eval(h_share);

%% C) Chebyshec Approximation
syms k h
syms x y
k=x;
h=y;
f=eval(f_e);
%f=eval(f_e) not used since it's constant
[x]=interpolation_cheb(0,10,20);
[y]=interpolation_cheb(0,10,20);
wk=eval(f);
%wk=ones(1,20);
xk=x;
yk=y;
m=20;
%% order 3
n=3;
coef=zeros(1,(n+1));
for i=1:n+1
    c=wk(1:m).*cos((i-1).*acos(xk(1,(1:m)))).*cos((i-1).*acos(yk(1,(1:m))));
    d1=cos((i-1)*acos(xk(1,(1:m)))).* cos((i-1)*acos(xk(1,(1:m))));
    d2=cos((i-1)*acos(yk(1,(1:m)))).* cos((i-1)*acos(yk(1,(1:m))));
    coef(1,i)=sum(c)./(sum(d1)*sum(d2)); 
end
syms x y
psii1=cos((0:n).*acos(2*((x)/(10))-1));
psii2=cos((0:n).*acos(2*((y)/(10))-1));
p=coef(1,(1:n+1)).*psii1(1,(1:n+1)).*psii2(1,(1:n+1));
poly_3=sum(p);

%plotting
[x y]=meshgrid(0:1:10);

figure(1);
ax1=subplot(1,3,2);
z_p=eval(poly_3);
mesh(ax1,x,y,z_p);
title(ax1,'approximation');

z_f=ones(11,11);
ax2=subplot(1,3,1);
mesh(ax2,x,y,z_f);
title(ax2,'function');

z_error=abs(z_f-z_p);
ax3=subplot(1,3,3);
mesh(ax3,x,y,z_error);
title(ax3,'error');

suptitle('3rd-order Chebyshev polynomial')

print -dpdf ex5c3_3.eps

%% order 4
n=4;
coef=zeros(1,(n+1));
for i=1:n+1
    c=wk(1:m).*cos((i-1).*acos(xk(1,(1:m)))).*cos((i-1).*acos(yk(1,(1:m))));
    d1=cos((i-1)*acos(xk(1,(1:m)))).* cos((i-1)*acos(xk(1,(1:m))));
    d2=cos((i-1)*acos(yk(1,(1:m)))).* cos((i-1)*acos(yk(1,(1:m))));
    coef(1,i)=sum(c)./(sum(d1)*sum(d2)); 
end
syms x y
psii1=cos((0:n).*acos(2*((x)/(10))-1));
psii2=cos((0:n).*acos(2*((y)/(10))-1));
p=coef(1,(1:n+1)).*psii1(1,(1:n+1)).*psii2(1,(1:n+1));
poly_4=sum(p);

%% order 5
n=5;
coef=zeros(1,(n+1));
for i=1:n+1
    c=wk(1:m).*cos((i-1).*acos(xk(1,(1:m)))).*cos((i-1).*acos(yk(1,(1:m))));
    d1=cos((i-1)*acos(xk(1,(1:m)))).* cos((i-1)*acos(xk(1,(1:m))));
    d2=cos((i-1)*acos(yk(1,(1:m)))).* cos((i-1)*acos(yk(1,(1:m))));
    coef(1,i)=sum(c)./(sum(d1)*sum(d2)); 
end
syms x y
psii1=cos((0:n).*acos(2*((x)/(10))-1));
psii2=cos((0:n).*acos(2*((y)/(10))-1));
p=coef(1,(1:n+1)).*psii1(1,(1:n+1)).*psii2(1,(1:n+1));
poly_5=sum(p);


%% order 6
n=6;
coef=zeros(1,(n+1));
for i=1:n+1
    c=wk(1:m).*cos((i-1).*acos(xk(1,(1:m)))).*cos((i-1).*acos(yk(1,(1:m))));
    d1=cos((i-1)*acos(xk(1,(1:m)))).* cos((i-1)*acos(xk(1,(1:m))));
    d2=cos((i-1)*acos(yk(1,(1:m)))).* cos((i-1)*acos(yk(1,(1:m))));
    coef(1,i)=sum(c)./(sum(d1)*sum(d2)); 
end
syms x y
psii1=cos((0:n).*acos(2*((x)/(10))-1));
psii2=cos((0:n).*acos(2*((y)/(10))-1));
p=coef(1,(1:n+1)).*psii1(1,(1:n+1)).*psii2(1,(1:n+1));
poly_6=sum(p);
%plotting

[x y]=meshgrid(0:1:10);

figure(4);
ax1=subplot(1,3,2);
z_p=eval(poly_6);
mesh(ax1,x,y,z_p);
title(ax1,'approximation');

z_f=ones(11,11);
ax2=subplot(1,3,1);
mesh(ax2,x,y,z_f);
title(ax2,'function');

z_error=abs(z_f-z_p);
ax3=subplot(1,3,3);
mesh(ax3,x,y,z_error);
title(ax3,'error');

suptitle('6th-order Chebyshev polynomial')

print -dpdf ex5c6_3.eps
%% order 7
n=7;
coef=zeros(1,(n+1));
for i=1:n+1
    c=wk(1:m).*cos((i-1).*acos(xk(1,(1:m)))).*cos((i-1).*acos(yk(1,(1:m))));
    d1=cos((i-1)*acos(xk(1,(1:m)))).* cos((i-1)*acos(xk(1,(1:m))));
    d2=cos((i-1)*acos(yk(1,(1:m)))).* cos((i-1)*acos(yk(1,(1:m))));
    coef(1,i)=sum(c)./(sum(d1)*sum(d2)); 
end
syms x y
psii1=cos((0:n).*acos(2*((x)/(10))-1));
psii2=cos((0:n).*acos(2*((y)/(10))-1));
p=coef(1,(1:n+1)).*psii1(1,(1:n+1)).*psii2(1,(1:n+1));
poly_7=sum(p);

%% order 8
n=8;
coef=zeros(1,(n+1));
for i=1:n+1
    c=wk(1:m).*cos((i-1).*acos(xk(1,(1:m)))).*cos((i-1).*acos(yk(1,(1:m))));
    d1=cos((i-1)*acos(xk(1,(1:m)))).* cos((i-1)*acos(xk(1,(1:m))));
    d2=cos((i-1)*acos(yk(1,(1:m)))).* cos((i-1)*acos(yk(1,(1:m))));
    coef(1,i)=sum(c)./(sum(d1)*sum(d2)); 
end
syms x y
psii1=cos((0:n).*acos(2*((x)/(10))-1));
psii2=cos((0:n).*acos(2*((y)/(10))-1));
p=coef(1,(1:n+1)).*psii1(1,(1:n+1)).*psii2(1,(1:n+1));
poly_8=sum(p);


%% order 9
n=9;
coef=zeros(1,(n+1));
for i=1:n+1
    c=wk(1:m).*cos((i-1).*acos(xk(1,(1:m)))).*cos((i-1).*acos(yk(1,(1:m))));
    d1=cos((i-1)*acos(xk(1,(1:m)))).* cos((i-1)*acos(xk(1,(1:m))));
    d2=cos((i-1)*acos(yk(1,(1:m)))).* cos((i-1)*acos(yk(1,(1:m))));
    coef(1,i)=sum(c)./(sum(d1)*sum(d2)); 
end
syms x y
psii1=cos((0:n).*acos(2*((x)/(10))-1));
psii2=cos((0:n).*acos(2*((y)/(10))-1));
p=coef(1,(1:n+1)).*psii1(1,(1:n+1)).*psii2(1,(1:n+1));
poly_9=sum(p);
%plotting
[x y]=meshgrid(0:1:10);
figure(7)

ax1=subplot(1,3,2);
z_p=eval(poly_9);
mesh(ax1,x,y,z_p);
title(ax1,'approximation');

z_f=ones(11,11);
ax2=subplot(1,3,1);
mesh(ax2,x,y,z_f);
title(ax2,'function');

z_error=abs(z_f-z_p);
ax3=subplot(1,3,3);
mesh(ax3,x,y,z_error);
title(ax3,'error');

suptitle('9th-order Chebyshev polynomial')

print -dpdf ex5c9_3.eps

%% order 10
n=10;
coef=zeros(1,(n+1));
for i=1:n+1
    c=wk(1:m).*cos((i-1).*acos(xk(1,(1:m)))).*cos((i-1).*acos(yk(1,(1:m))));
    d1=cos((i-1)*acos(xk(1,(1:m)))).* cos((i-1)*acos(xk(1,(1:m))));
    d2=cos((i-1)*acos(yk(1,(1:m)))).* cos((i-1)*acos(yk(1,(1:m))));
    coef(1,i)=sum(c)./(sum(d1)*sum(d2)); 
end
syms x y
psii1=cos((0:n).*acos(2*((x)/(10))-1));
psii2=cos((0:n).*acos(2*((y)/(10))-1));
p=coef(1,(1:n+1)).*psii1(1,(1:n+1)).*psii2(1,(1:n+1));
poly_15=sum(p);

%% order 11
n=11;
coef=zeros(1,(n+1));
for i=1:n+1
    c=wk(1:m).*cos((i-1).*acos(xk(1,(1:m)))).*cos((i-1).*acos(yk(1,(1:m))));
    d1=cos((i-1)*acos(xk(1,(1:m)))).* cos((i-1)*acos(xk(1,(1:m))));
    d2=cos((i-1)*acos(yk(1,(1:m)))).* cos((i-1)*acos(yk(1,(1:m))));
    coef(1,i)=sum(c)./(sum(d1)*sum(d2)); 
end
syms x y
psii1=cos((0:n).*acos(2*((x)/(10))-1));
psii2=cos((0:n).*acos(2*((y)/(10))-1));
p=coef(1,(1:n+1)).*psii1(1,(1:n+1)).*psii2(1,(1:n+1));
poly_11=sum(p);


%% order 12
n=12;
coef=zeros(1,(n+1));
for i=1:n+1
    c=wk(1:m).*cos((i-1).*acos(xk(1,(1:m)))).*cos((i-1).*acos(yk(1,(1:m))));
    d1=cos((i-1)*acos(xk(1,(1:m)))).* cos((i-1)*acos(xk(1,(1:m))));
    d2=cos((i-1)*acos(yk(1,(1:m)))).* cos((i-1)*acos(yk(1,(1:m))));
    coef(1,i)=sum(c)./(sum(d1)*sum(d2)); 
end
syms x y
psii1=cos((0:n).*acos(2*((x)/(10))-1));
psii2=cos((0:n).*acos(2*((y)/(10))-1));
p=coef(1,(1:n+1)).*psii1(1,(1:n+1)).*psii2(1,(1:n+1));
poly_12=sum(p);
%plotting
[x y]=meshgrid(0:1:10);
figure(10);

ax1=subplot(1,3,2);
z_p=eval(poly_12);
mesh(ax1,x,y,z_p);
title(ax1,'approximation');

z_f=ones(11,11);
ax2=subplot(1,3,1);
mesh(ax2,x,y,z_f);
title(ax2,'function');

z_error=abs(z_f-z_p);
ax3=subplot(1,3,3);
mesh(ax3,x,y,z_error);
title(ax3,'error');

suptitle('12th-order Chebyshev polynomial');

print -dpdf ex5c12_3.eps
%% order 13
n=13;
coef=zeros(1,(n+1));
for i=1:n+1
    c=wk(1:m).*cos((i-1).*acos(xk(1,(1:m)))).*cos((i-1).*acos(yk(1,(1:m))));
    d1=cos((i-1)*acos(xk(1,(1:m)))).* cos((i-1)*acos(xk(1,(1:m))));
    d2=cos((i-1)*acos(yk(1,(1:m)))).* cos((i-1)*acos(yk(1,(1:m))));
    coef(1,i)=sum(c)./(sum(d1)*sum(d2)); 
end
syms x y
psii1=cos((0:n).*acos(2*((x)/(10))-1));
psii2=cos((0:n).*acos(2*((y)/(10))-1));
p=coef(1,(1:n+1)).*psii1(1,(1:n+1)).*psii2(1,(1:n+1));
poly_13=sum(p);
%plotting
[x y]=meshgrid(0:1:10);
figure(11)
z_p=eval(poly_13);
mesh(x,y,z_p);
hold on
z_f=eval(f);
mesh(x,y,z_f);
legend('Function approx order 13', 'Function')
hold off

%% order 14
n=14;
coef=zeros(1,(n+1));
for i=1:n+1
    c=wk(1:m).*cos((i-1).*acos(xk(1,(1:m)))).*cos((i-1).*acos(yk(1,(1:m))));
    d1=cos((i-1)*acos(xk(1,(1:m)))).* cos((i-1)*acos(xk(1,(1:m))));
    d2=cos((i-1)*acos(yk(1,(1:m)))).* cos((i-1)*acos(yk(1,(1:m))));
    coef(1,i)=sum(c)./(sum(d1)*sum(d2)); 
end
syms x y
psii1=cos((0:n).*acos(2*((x)/(10))-1));
psii2=cos((0:n).*acos(2*((y)/(10))-1));
p=coef(1,(1:n+1)).*psii1(1,(1:n+1)).*psii2(1,(1:n+1));
poly_14=sum(p);
%plotting
[x y]=meshgrid(0:1:10);
figure(12)
z_p=eval(poly_14);
mesh(x,y,z_p);
hold on
z_f=eval(f);
mesh(x,y,z_f);
legend('Function approx order 14', 'Function')
hold off

%% order 15
n=15;
coef=zeros(1,(n+1));
for i=1:n+1
    c=wk(1:m).*cos((i-1).*acos(xk(1,(1:m)))).*cos((i-1).*acos(yk(1,(1:m))));
    d1=cos((i-1)*acos(xk(1,(1:m)))).* cos((i-1)*acos(xk(1,(1:m))));
    d2=cos((i-1)*acos(yk(1,(1:m)))).* cos((i-1)*acos(yk(1,(1:m))));
    coef(1,i)=sum(c)./(sum(d1)*sum(d2)); 
end
syms x y
psii1=cos((0:n).*acos(2*((x)/(10))-1));
psii2=cos((0:n).*acos(2*((y)/(10))-1));
p=coef(1,(1:n+1)).*psii1(1,(1:n+1)).*psii2(1,(1:n+1));
poly_15=sum(p);
%plotting
[x y]=meshgrid(0:1:10);
figure(13)

ax1=subplot(1,3,2);
z_p=eval(poly_15);
mesh(ax1,x,y,z_p);
title(ax1,'approximation');

z_f=ones(11,11);
ax2=subplot(1,3,1);
mesh(ax2,x,y,z_f);
title(ax2,'function');

z_error=abs(z_f-z_p);
ax3=subplot(1,3,3);
mesh(ax3,x,y,z_error);
title(ax3,'error');

suptitle('15th-order Chebyshev polynomial')

print -dpdf ex5c15.eps


%% D) Isoquants
%isoquants with exact function and approximation of order 3
%function contour
[K,H] = meshgrid(0:10, 0:10);
Z=subs(f_e,'h','H');
Z=subs(Z,'k','K');
Z =1./(1./(2.*H.^(1./99)) + 1./(2.*K.^(1/99))).^99;
[K,H] = meshgrid(0:10, 0:10);
[C,h] = contour(K,H,Z);
print -dpdf ex5d1_3.eps

%approx
poly_3s=simplify(poly_3);
poly_3s=subs(poly_3s,'x','K');
poly_3s=subs(poly_3s,'y','H');
Z3=simplify(poly_3s);
Z3 =(H./5 - 1).*(K./450 - 1./90) + (3080777541025067.*((2.*H.^2)./25 - (4.*H)./5 + 1).*((2.*K.^2)./25 - (4.*K)./5 + 1))./73786976294838206464 + (4951920775910657.*((4.*H.^3)./125 - (12.*H.^2)./25 + (9.*H)./5 - 1).*((4.*K.^3)./125 - (12.*K.^2)./25 + (9.*K)./5 - 1))./37778931862957161709568 + 1/4;
[K,H] = meshgrid(0:10, 0:10);
[C,h] = contour(K,H,Z3);

print -dpdf ex5d2_3.eps

%approx vs function contour
[K,H] = meshgrid(0:10, 0:10);
[C,h] = contour(K,H,Z3);
hold on
[C,h] = contour(K,H,Z);
print -dpdf ex5d3_3.eps

%error contour
e=simplify(abs(f_e-poly_3));
e=subs(e,'x','K');
e=subs(e,'y','H');
e=subs(e,'k','K');
e=subs(e,'h','H');
e = abs((H./5 - 1).*(K./450 - 1./90) + (3080777541025067.*((2.*H.^2)./25 - (4.*H)./5 + 1).*((2.*K.^2)./25 - (4.*K)./5 + 1))./73786976294838206464 + (4951920775910657.*((4.*H.^3)./125 - (12.*H.^2)./25 + (9.*H)./5 - 1).*((4.*K.^3)./125 - (12.*K.^2)./25 + (9.*K)./5 - 1))./37778931862957161709568 - (633825300114114700748351602688.*H.*K)./(H.^(1./99) + K.^(1./99)).^99 + 1/4);
[K,H] = meshgrid(0:10, 0:10);
[C,h] = contour(K,H,e);
hold on
[C,h] = contour(K,H,Z3);
print -dpdf ex5d4_3.eps
hold off
