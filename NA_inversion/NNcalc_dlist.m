%
%-----------------------------------------------------------------------
%
%       Subroutine NNcalc_dlist - calculates square of distance from
%                                 all base points to new axis (defined
%                                 by dimension dim through point x.
%                                 It also updates the nearest node and
%                                 distance to the point x.
%
%       This is a full update of dlist, i.e. not using a previous dlist.
%
%-----------------------------------------------------------------------
%
function [dlist,nodex]=NNcalc_dlist(dim,bp,nd,nb,x)

dmin = 0;
for j=1:dim-1
    d = (x(j)-bp(j,1));
    d = d*d;
    dmin = dmin + d;
end
for j=dim+1:nd
    d = (x(j)-bp(j,1));
    d = d*d;
    dmin = dmin + d;
end
dlist(1) = dmin;
d = (x(dim)-bp(dim,1));
d = d*d;
dmin = dmin + d;
nodex = 1;

for i=2:nb
    dsum = 0;
    for j=1:dim-1
        d = (x(j)-bp(j,i));
        d = d*d;
        dsum = dsum + d;
    end
    for j=dim+1:nd
        d = (x(j)-bp(j,i));
        d = d*d;
        dsum = dsum + d;
    end
    dlist(i) = dsum;
    d = (x(dim)-bp(dim,i));
    d = d*d;
    dsum = dsum + d;
    if (dmin>dsum)
        dmin = dsum;
        nodex = i;
    end
end

end
