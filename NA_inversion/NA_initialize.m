%----------------------------------------------------------------------------

%       NA_initialize - performs minor initialization tasks for NA algorithm.

%       Calls no other routines.

%					M. Sambridge, Oct. 1996

%----------------------------------------------------------------------------

function [restartNA,ranget,x]=NA_initialize(nd,nd_max,range,scales,nsample,ncells)

% initialize variables
range=zeros(2,nd);
ranget=zeros(2,nd);
x=zeros(1,nd);
rval=zeros(1,2);

global verbose debug summary idnext ic sobol sob_seq

if (nd>nd_max)
    dips('Error in subroutine NA_initialize. Number of dimensions greater than maximum.')
    error('NA:ndmax','Increase parameter nd_max')
end

%	set logical switch for first call to NA_sample
%	(ensures distance list is initialized)
restartNA = 1;

%	set initial parameter for NA walk
ic = 1;

if (sobol)
    % create Sobol sequence with nd*nsample numbers
    sob_seq=NA_sobol(nd,itmax*nsample,0);
end
%	Normalize parameter ranges
%	by a-priori model co-variances
if (scales(1)==0)
    %	First option:
    %	No transform (All a priori model co-variances are equal to unity)
    for i=1:nd
        ranget(1,i) = range(1,i);
        ranget(2,i) = range(2,i);
        scales(i+1) = 1;
    end
elseif (scales(1)==-1.0)
    %	Second option:
    %	Use parameter range as a priori model co-variances
    for i=1:nd
        ranget(1,i) = 0;
        ranget(2,i) = 1;
        scales(i+1) = range(2,i)-range(1,i);
    end
else
    %	Third option:
    %	Use scales array as
    %	a priori model co-variances
    for i=1:nd
        if (scales(i+1)==0)
            disp(' Error in subroutine NA_initialize.')
            error('NA:prior_covar','Input a priori model co-variance is equal to zero.')
        end
        ranget(1,i)  = 0;
        ranget(2,i)  = (range(2,i)-range(1,i))/scales(i+1);
    end
end
%	calculate axis increments and initialize current point
%       (used by NA_sample) to mid-point of parameter space

for i=1:nd
    x(i) = (ranget(2,i)+ranget(1,i))/2;
end

end