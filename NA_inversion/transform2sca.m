%----------------------------------------------------------------------------

%       transform2sca - transforms model from raw to scaled units.

%Input:
%      nd		: dimension of parameter space
%      model_raw(nd)	: model in raw co-ordinates
%      range(2,nd)	: min and max of parameter space
%			  in raw co-ordinates.
%      scales(nd+1)	: range scale factors

%Output:
%      model_sca(nd)	: model in scaled co-ordinates

%Comments:
%         This routine transforms a model in raw co-ordinates
%	 to dimensionless units determined by parameter scale factors.

%       Calls no other routines.

%                                               M. Sambridge, March 1998

%----------------------------------------------------------------------------
function [model_sca]=transform2sca(model_raw,nd,range,scales)


if (scales(1)==0)
    for i=1:nd
        model_sca(i) = model_raw(i);
    end
elseif (scales(1)==-1)
    for i=1:nd
        model_sca(i) = (model_raw(i)-range(1,i))/(range(2,i)-range(1,i));
    end
else
    for i=1:nd
        model_sca(i) = (model_raw(i)-range(1,i))/scales(i+1);
    end
    
end
end