%% Multivariate CES Function Approximation - sigma=0.25
clc
clear all;

syms k h alpha sigma
f=((1-alpha)*k^((sigma-1)/sigma)+alpha*h^((sigma-1)/sigma))^(sigma/(sigma-1));
alpha=0.5;
sigma=0.25;
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
[x]=interpolation_cheb(0,10,20);
[y]=interpolation_cheb(0,10,20);
wk=eval(f);
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

z_f=eval(f);
ax2=subplot(1,3,1);
mesh(ax2,x,y,z_f);
title(ax2,'function');

z_error=abs(z_f-z_p);
ax3=subplot(1,3,3);
mesh(ax3,x,y,z_error);
title(ax3,'error');

suptitle('3rd-order Chebyshev polynomial')

print -dpdf ex5c3.eps

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
%plotting
[x y]=meshgrid(0:1:10);
figure(2)
z_p=eval(poly_4);
mesh(x,y,z_p);
hold on
z_f=eval(f);
mesh(x,y,z_f);
legend('Function approx order 4', 'Function')
hold off

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
%plotting
[x y]=meshgrid(0:1:10);
figure(3)
z_p=eval(poly_5);
mesh(x,y,z_p);
hold on
z_f=eval(f);
mesh(x,y,z_f);
legend('Function approx order 5', 'Function')
hold off

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

z_f=eval(f);
ax2=subplot(1,3,1);
mesh(ax2,x,y,z_f);
title(ax2,'function');

z_error=abs(z_f-z_p);
ax3=subplot(1,3,3);
mesh(ax3,x,y,z_error);
title(ax3,'error');

suptitle('6th-order Chebyshev polynomial')

print -dpdf ex5c6.eps
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
%plotting
[x y]=meshgrid(0:1:10);
figure(5)
z_p=eval(poly_7);
mesh(x,y,z_p);
hold on
z_f=eval(f);
mesh(x,y,z_f);
legend('Function approx order 7', 'Function')
hold off

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
%plotting
[x y]=meshgrid(0:1:10);
figure(6)
z_p=eval(poly_8);
mesh(x,y,z_p);
hold on
z_f=eval(f);
mesh(x,y,z_f);
legend('Function approx order 8', 'Function')
hold off

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

z_f=eval(f);
ax2=subplot(1,3,1);
mesh(ax2,x,y,z_f);
title(ax2,'function');

z_error=abs(z_f-z_p);
ax3=subplot(1,3,3);
mesh(ax3,x,y,z_error);
title(ax3,'error');

suptitle('9th-order Chebyshev polynomial')

print -dpdf ex5c9.eps

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
%plotting
[x y]=meshgrid(0:1:10);
figure(8)
z_p=eval(poly_10);
mesh(x,y,z_p);
hold on
z_f=eval(f);
mesh(x,y,z_f);
legend('Function approx order 10', 'Function')
hold off

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
%plotting
[x y]=meshgrid(0:1:10);
figure(9)
z_p=eval(poly_11);
mesh(x,y,z_p);
hold on
z_f=eval(f);
mesh(x,y,z_f);
legend('Function approx order 11', 'Function')
hold off

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

z_f=eval(f);
ax2=subplot(1,3,1);
mesh(ax2,x,y,z_f);
title(ax2,'function');

z_error=abs(z_f-z_p);
ax3=subplot(1,3,3);
mesh(ax3,x,y,z_error);
title(ax3,'error');

suptitle('12th-order Chebyshev polynomial');

print -dpdf ex5c12.eps
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

z_f=eval(f);
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
Z=subs(f_e,'h','H');
Z=subs(Z,'k','K');
[K,H] = meshgrid(0:10, 0:10);
Z=(0.5.*(K).^((-0.75)/0.25)+(0.5).*(H).^(-0.75/0.25)).^(0.25/-0.75);
[C,h] = contour(K,H,Z);
print -dpdf ex5d1.eps
%approx
poly_3s=simplify(poly_3);
poly_3s=subs(poly_3s,'x','K');
poly_3s=subs(poly_3s,'y','H');
Z3=simplify(poly_3s);
Z3=(H/5 - 1)*(K/450 - 1/90) + (1540388770512533*((2*H^2)/25 - (4*H)/5 + 1)*((2*K^2)/25 - (4*K)/5 + 1))/36893488147419103232 + (19343440530901*((4*H^3)/125 - (12*H^2)/25 + (9*H)/5 - 1)*((4*K^3)/125 - (12*K^2)/25 + (9*K)/5 - 1))/147573952589676412928 + 1/4
[K,H] = meshgrid(0:10, 0:10);
[C,h] = contour(K,H,Z3);
print -dpdf ex5d2.eps

%approx vs function contour
[K,H] = meshgrid(0:10, 0:10);
[C,h] = contour(K,H,Z3);
hold on
[C,h] = contour(K,H,Z);
print -dpdf ex5d3.eps

%error contour
e=simplify(abs(f_e-poly_3));
e=subs(e,'x','K');
e=subs(e,'y','H');
e=subs(e,'k','K');
e=subs(e,'h','H');
e=abs((H./5 - 1).*(K./450 - 1/90) + (1540388770512533.*((2.*H.^2)./25 - (4.*H)./5 + 1).*((2.*K.^2)./25 - (4.*K)./5 + 1))./36893488147419103232 - 2^(1/3)./((H.^3 + K.^3)./(H.^3.*K.^3)).^(1/3) + (19343440530901.*((4.*H.^3)./125 - (12.*H.^2)./25 + (9.*H)./5 - 1).*((4.*K.^3)./125 - (12.*K.^2)./25 + (9.*K)./5 - 1))./147573952589676412928 + 1/4)
[K,H] = meshgrid(0:10, 0:10);
[C,h] = contour(K,H,e);
hold on
[C,h] = contour(K,H,Z3);
print -dpdf ex5d4.eps
hold off