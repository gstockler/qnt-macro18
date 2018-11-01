function [a_min]=borrowing_const(i,y,r)
abar=min(y)/r;
if i==1
    a_min=-abar;
else
    a_min=0;
end


end
