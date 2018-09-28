function [x]=interpolation_n(a,b,n)
x=zeros(1,n);
 for i=1:n 
    x(1,i)=a+((i-1)/(n-1))*(b-a);
 end

end




