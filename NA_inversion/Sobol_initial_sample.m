%----------------------------------------------------------------------------

%       NA_initial_sample - generates initial sample for NA algorithm.

%Comments:

%	Assumes n-dimensional Sobol sequence has been initialized
%	Will generate a minimum of two samples.

%	Assumes ran3 has been initialized.

%       Calls no other routines.

%					M. Sambridge, Oct. 1996
%					(Updated for ran3 Aug. 1997)

%----------------------------------------------------------------------------

function [na_models,misfit]=Sobol_initial_sample(nsample,nd,range)

    %	Generate initial uniform
    %	random sample using a
    %       uniform random distribution
    
    %	Generate initial samples
    %	using quasi random sequences
    if (sobol)
        ran_matrix=sob_seq(:,1:nsample);
        sob_seq(:,1:nsample)=[];
        for i=1:nsample
            na_models(:,i)=(1-ran_matrix(:,i)).*range(1,:)'+ran_matrix(:,i).*range(2,:)';
        end
    else
        %	Use pseudo random number generator rand()
        ran_matrix=rand(nd,nsample);
        for i=1:nsample
            na_models(:,i)=(1-ran_matrix(:,i)).*range(1,:)'+ran_matrix(:,i).*range(2,:)';
        end
    end
    misfit=zeros(1,nsample);

end