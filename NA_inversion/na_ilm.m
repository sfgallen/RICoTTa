%*************************************************************************

%	Main subroutine NA  - sampling a parameter space
%		      	      using a Neighbourhood algorithm

%Input files:
%	    na.in	- Options for Neighbourhood algorithm

%Output files:
%	    na.sum	- summary of results
%	    na.nad	- NAD (NA direct access file of models)
%	    sobol.coeff	- initializing data used for
%			  quasi-random sequences

%Comments:
%	 The NAD file is a direct access compact format file
%	 containing all models generated by the neighbourhood
%	 algorithm. (NAD files can be read in by multi-dimensional
%	 integration program NA-Bayes and plot program S-plot.)

%				M. Sambridge, (RSES, ANU)
%				Last revision Sept. 1999.

%*************************************************************************

% test script for NA
ndin=3;
rangein=zeros(2,ndin);
rangein(1,:)=[0 0 10];
rangein(2,:)=[100 2 50];

% load boundaries of input parameters
load na_param_inc.mat

% Info and Logical unit common blocks used by NA routines
global  verbose debug summary nxsave ndsave ndc nerr ncald nupd cells torder taxis tup tcd tdev tna tres sobol sob_seq

%	User specific setup for forward modelling
nd=ndin;
range(:,1:nd)=rangein(:,1:nd);
%call user_init(nd,range,scales)

% Read in options for Neighbourhood algorithm.
[monte,istype,nsleep,noforward,nclean,nsample,nsamplei,itmax,ncells]=NA_options(nsample_max, nit_max, nmod_max,...
    nsleep_max, nd);


% set other info or debug options
check = 1;
check = 0;
scales(1:nd)=-1; % 0: No transform (All a priori model co-variances are equal to unity); 1: Use parameter range as a priori model co-variances

% Initialize NA routines.
[restartNA,ranget,xcur]=NA_initialize(nd,nd_max,range,scales,nsample,ncells);

% Generate or read in starting models
[new_models,new_misfit]=NA_initial_sample(istype,monte,nsample,nd,range,scales);

%	MAIN OPTIMIZATION LOOP
misfit=[];
models=[];

for it = 1:itmax+1
    %	Calculate misfit values for each model in the current population.
    disp(['Start forward modeling iteration: ',num2str(it)])
    
    %	skip forward modelling for starting models if they have been read in from a NAD file
    if (it>1 || noforward~=1)
        % solve forward problem
        %new_misfit=forward(new_models);
        new_misfit=abs(new_models(1,:)-50);
        misfit = [misfit new_misfit];
        models = [models new_models];
    else
        misfit = [misfit new_misfit];
        models = [models new_models];
    end
    
    %	Calculate properties of current misfit distribution.
    %	(Mean,min,best model etc.)
    [mfitmean,mfitminc,mfitmin,mfitord,mopt]=NA_misfits(misfit);
    
    % Best model
    disp(['Best model parameters: ',num2str(models(:,mopt)')])
    disp(['Lowest misfit: ',num2str(mfitmin)])
    
    ntot=length(misfit);
    calcmovement=0;
    
    % tranform to scale
    for i=1:ntot
        models_sca(:,i)=transform2sca(models(:,i),nd,range,scales);
    end
    
    if (it<itmax+1)
        %	Call main NA routines
        if (monte)
            % Perform Monte Carlo search for comparison to NA.
            [new_models_sca]=NA_random(nd,range,nsample);
        else
            % generate a new sample using Neighbourhood algorithm (resample version)
            [new_models_sca,xcur,restartNA,work_NA1]=NA_sample(models_sca,ntot,nsample,nd,nsleep,ncells,misfit,mfitord,ranget,check,xcur,calcmovement,nclean);
        end
    end
    % transform to raw
    for i=1:nsample
        new_models(:,i)=transform2raw(nd,range,scales,new_models_sca(:,i));
    end

end
