%----------------------------------------------------------------------------

%       transform2raw - transforms model from scaled to raw units.

%Input:
%      nd		: dimension of parameter space
%      model_sca(nd)	: model in scaled co-ordinates
%      range(2,nd)	: min and max of parameter space
%			  in raw co-ordinates.
%      scales(nd+1)	: range scale factors

%Output:
%      model_raw(nd)	: model in scaled co-ordinates

%Comments:
%         This routine transforms a model in dimensionless scaled
%	 co-ordinates to input (raw) units.

%       Calls no other routines.

%                                               M. Sambridge, March 1998

%----------------------------------------------------------------------------

function model_raw=transform2raw(nd,range,scales,model_sca)

model_raw=zeros(1,nd);

if (scales(1)==0)
    
    for i=1:nd
        model_raw(i) = model_sca(i);
    end
    
elseif (scales(1)==-1)
    
    for i=1:nd
        b = model_sca(i);
        a = 1-b;
        model_raw(i) = a*range(1,i) + b*range(2,i);
    end
    
else
    
    for i=1:nd
        model_raw(i) = range(1,i) + scales(i+1)*model_sca(i);
    end
    
end

end