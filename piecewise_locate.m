function [loc]=piecewise_locate(a_nk,a)
nk=length(a_nk);
for ind=1:nk
    if (a_nk(1,ind)<=a) & (a<=a_nk(1,ind+1))
        loc = ind;
        return
    end
end
    
end