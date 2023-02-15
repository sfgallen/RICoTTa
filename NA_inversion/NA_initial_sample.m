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

function [na_models,misfit]=NA_initial_sample(istype,monte,nsample,nd,range,scales)

load na_param_inc.mat

% initialize variables
global sobol sob_seq

%	read in starting models from a NAD file
if (istype==1 && monte==0)
    load na_nad.mat
    % file with variables:
    %             nd               : dimension of parameter space
    %             ne               : number of models in ensemble
    %             misfit(nd)       : array of misfit values for each model (real*4)
    %             na_models(ne,nd) : array of model values  (real*4)
    for i=1:nsample
        [na_models]=transform2sca(models(1,i),nd,range,scales);
    end
elseif (istype==2 && monte==0)
    %	read in models from ascii input file
    load models_in.mat % reads misift, na_models lu_nad=fopen('models.in');
    nsample_r = length(misfit);
    
    %	scale models read in
    for i=1:nsample_r
        [na_models]=transform2sca(models(1,i),nd,range,scales);
    end
    
    if (nsample_r<nsample)
        %	randomly generate the the remaining samples
        nrem = nsample - nsample_r;
        disp([' generating remaining samples ',num2str(nrem)])
        %	Generate initial samples
        %	using quasi random sequences
        if (sobol)
            ran_matrix=NA_sobol(nd,nrem,100);
            for i=nsample_r+1:nsample
                na_models(:,i)=(1-ran_matrix(:,i)).*range(1,:)'+ran_matrix(:,i).*range(2,:)';
            end
        else
            %	se pseufor random number generator ran3
            ran_matrix=rand(nd,nrem);
            for i=nsample_r+1:nsample
                na_models(:,i)=(1-ran_matrix(:,i)).*range(1,:)'+ran_matrix(:,i).*range(2,:)';
            end
        end
    end
else
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

end