function [poly]=cheb_poly(func,a,b,n,m)
syms x
f=func;
[x]=interpolation_cheb(a,b,m);
wk=eval(f);
xk=x;
coef=zeros(1,(n+1));
for i=1:n+1 %1st is 0-rder
    c=wk(1:m).*cos((i-1).*acos(xk(1,(1:m))));
    d=cos((i-1)*acos(xk(1,(1:m)))).* cos((i-1)*acos(xk(1,(1:m))));
    coef(1,i)=sum(c)./sum(d); 
end
syms x
psii=cos((0:n).*acos(2*((x-a)/(b-a))-1));
p=coef(1,(1:n+1)).*psii(1,(1:n+1));
poly=sum(p);
end
