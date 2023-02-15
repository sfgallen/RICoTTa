%----------------------------------------------------------------------------

%       NA_random - generates a uniform random sample to be used
%	    instead of NA algorithm.

%       Assumes random number generator ran3 has been initialized

%       Calls ran3.

%					M. Sambridge, Oct. 1996

%----------------------------------------------------------------------------

function [na_models]=NA_random(nd,range,nsample)

% generate samples
if (sobol)
    % use quasi-random number generator
    ran_matrix=NA_sobol(nd,nsample,100);
    for i=1:nsample
        na_models(:,i)=(1-ran_matrix(:,i)).*range(1,:)'+ran_matrix(:,i).*range(2,:)';
    end
else
    % use pseudo-random number generator
    ran_matrix=rand(nd,nrem);
    for i=1:nsample
        na_models(:,i)=(1-ran_matrix(:,i)).*range(1,:)'+ran_matrix(:,i).*range(2,:)';
    end
end

end