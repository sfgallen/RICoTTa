%-----------------------------------------------------------------------
%
%Subroutine NNupdate_dlist - calculates square of distance from
%			     all base points to new axis, assuming
%                                    dlist contains square of all distances
%			     to previous axis dimlast. It also
%			     updates the nearest node to the
%			     point x through which the axes pass.
%
%
%-----------------------------------------------------------------------
%
function [dlist,node,dmin]=NNupdate_dlist(dim,dimlast,dlist,bp,nd,nb,x)

d1 = (x(dimlast)-bp(dimlast,1));
d1 = d1*d1;
dmin = dlist(1)+d1;
node = 1;
d2 = (x(dim)-bp(dim,1));
d2 = d2*d2;
dlist(1) = dmin-d2;
for i=2:nb
    d1 = (x(dimlast)-bp(dimlast,i));
    ds = d1;
    d1 = dlist(i)+d1*d1;
    if (dmin>d1)
        dmin = d1;
        node = i;
    end
    d2 = (x(dim)-bp(dim,i));
    d2 = d2*d2;
    dlist(i) = d1-d2;
end

end