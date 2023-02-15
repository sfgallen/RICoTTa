%
%-----------------------------------------------------------------------
%
%NNaxis_intersect - find intersections of current Voronoi cell
%		   with current 1-D axis.
%
%       Input:
%      x(nd)		:point on axis
%      dim          :dimension index (defines axis)
%      dlist		:set of distances of base points to axis
%      bp(nd,nb)	:set of base points
%      nd           :number of dimensions
%      nb           :number of base points
%      resetlist	:TRUE if dlist and nodex is to be calculated
%      nodex		:index of base node closest to x
%      dmin_in		:distance of base node closest to x
%      xmin         :start point along axis
%      xmax         :end point along axis

%       Output:
%      x1		:intersection of first Voronoi boundary
%      x2		:intersection of second Voronoi boundary

%       Comment:
%        This method uses a simple formula to exactly calculate
%	the intersections of the Voronoi cells with the 1-D axis.
%	It makes use of the perpendicluar distances of all nodes
%	to the current axis contained in the array dlist.

%        The method involves a loop over ensemble nodes for
%	each new intersection found. For an axis intersected
%	by ni Voronoi cells the run time is proportional to ni*ne.

%	It is assumed that the input point x(nd) lies in
%	the Vcell of nodex, i.e. nodex is the closest node to x(nd).

%	Note: If the intersection points are outside of either
%	      axis range  the axis range is returned, i.e.

%	      		x1 is set to max(x1,xmin) and
%	      		x2 is set to min(x2,xmin) and

%                                       M. Sambridge, RSES, June 1998

%-----------------------------------------------------------------------
%
function [x1,x2]=NNaxis_intersect(x,dim,dlist,bp,nd,nb,nodex,xmin,xmax)

x1 = xmin;
x2 = xmax;
dp0   = dlist(nodex);
x0    = bp(dim,nodex);

% find intersection of current Voronoi cell with 1-D axis
for j=1:nodex-1
    xc    = bp(dim,j);
    dpc   = dlist(j);
    % calculate intersection of interface (between nodes nodex and j) and 1-D axis
    dx = x0 - xc;
    if (dx~=0)
        xi = 0.5*(x0+xc+(dp0-dpc)/dx);
        if (xi>xmin && xi<xmax)
            if (xi>x1 && x0>xc)
                x1 = xi;
            elseif (xi<x2 && x0<xc)
                x2 = xi;
            end
        end
    end
end

for j=nodex+1:nb
    xc    = bp(dim,j);
    dpc   = dlist(j);
    % calculate intersection of interface (between nodes nodex and j) and 1-D axis
    dx = x0 - xc;
    if (dx~=0)
        xi = 0.5*(x0+xc+(dp0-dpc)/dx);
        if (xi>xmin && xi<xmax)
            if (xi>x1 && x0>xc)
                x1 = xi;
            elseif (xi<x2 && x0<xc)
                x2 = xi;
            end
        end
    end
end

end
