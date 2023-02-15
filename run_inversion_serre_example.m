% clear workpace
clear; clc

% make sure new matlab session use different random numbers
rng shuffle

%% add paths needed to exicute code
addpath([cd,'\river_models']); %add the path to the river models
addpath([cd,'\thermo_models']); %add the path to the thermo models
addpath([cd,'\NA_inversion']); %add the path to the neighborhood algorithum
addpath([cd,'\forward_model_functions']) % add path to function to exicute forward models
addpath([cd,'\Calabria_data']); % add path to the folder with the pre-processed data

% name added to files produced at the end of the inversion
ftag = 'serre_example';

%% load the data
% thermo data
load thermo_serre.mat

% marine terrace data
mt = load('serre_MT.txt');

% load stream data
load stream_param_serre.mat

%% add marine terrace data to parameters
% define marine terrace parameters
param.mt_age = mt(:,1);
param.mt_u = mt(:,2);
param.mt_uc = mt(:,3);

%% User defined parameters

% define time variables
param.start_t = 30;     % model run time in Myr
param.dt = 1e3;         % model time step in years

% stream parameters
param.ee = 15;          % error on observed stream elevations

% thermo
param.T0=20;                  % surface temperature at sea level in C
param.lr=5;                   % atmospheric lapse rate in C/km
param.TD=1E-6;                % thermal diffusivity in m2/s
param.rho_c=2700;             % crustal density in kg/m3
param.hpc=9.6E-10*param.rho_c;% crustal heat production in W/m3
param.cp_c=800;               % specific heat capacity of granite in J/kg*K
param.gg = 25;                % initial geothermal gradient

% other variables
param.kflag = 0;

% define overall parameters
param.nu=4;                   % number of uplift steps
param.nux=1;                  % 1: lateral constant uplift, 

%% additional parameters needed for forward model

% additional thermo parameters
param.thermo_meas = thermo_meas;
param.nr_thermo = length(thermo_meas.x);

% sample = [Time1, Time2, Time3, U1, U2, U3, U4, K, n, m/n]
sample = [20.24,9.36,1.00,1.63,.25,0.036,0.58,1.43e-6,1.22,0.46];
corr_sample_size(sample,param.nu,param.nux,param.kflag)

% Here 'sample' is roughly the MAP solution from pervious model runs
MAP = sample;


%% inversion of data with NA (Neighbourhood Algorithm)
ndin=10;

% define priors
rangein=zeros(2,ndin);
% t1 t2 t3 u1 u2 u3 u4 k n m/n 
rangein(1,:)=[13 5 0.25 0.25 0 0 0 1E-8 1 0.3];
rangein(2,:)=[25 13 5 2 0.5 1.0 2 1E-5 1.5 0.7];
resolution=[0.1 0.1 0.05 0.01 0.002 0.005 0.005 5E-8 0.01 0.01]; % if the range of resampled cells < resolution, than the inverison is starting from a new initial ensemble
inv_method=1; % inversion method: 1 (only resampling of existing ensemble), 0 (resampling and additional random samples are generted)

% load boundaries of input parameters
load na_param_inc.mat

% Info and Logical unit common blocks used by NA routines
global  inv_method verbose debug summary nxsave ndsave ndc nerr ncald nupd cells torder taxis tup tcd tdev tna tres sobol sob_seq

%	User specific setup for forward modelling
nd=ndin;
range(:,1:nd)=rangein(:,1:nd);
%call user_init(nd,range,scales)

% Read in options for Neighbourhood algorithm.
[monte,istype,nsleep,noforward,nclean,nsample,nsamplei,itmax,ncells]=NA_options(nsample_max, nit_max, nmod_max,...
    nsleep_max, nd);


% set other info or debug options
check = 0;
scales(1:nd)=-1; % 0: No transform (All a priori model co-variances are equal to unity); 1: Use parameter range as a priori model co-variances

% Initialize NA routines.
[restartNA,ranget,xcur]=NA_initialize(nd,nd_max,range,scales,nsample,ncells);

% Generate or read in starting models
[new_models,new_misfit]=NA_initial_sample(istype,monte,nsample,nd,range,scales);
new_models(:,1) = sample';

% initialize variables
misfit_all=nan(1,nsample*(itmax+1));
models_all=nan(ndin,nsample*(itmax+1));
logL_topo_all=nan(1,nsample*(itmax+1));
logL_thermo_all=nan(1,nsample*(itmax+1));
logL_MT_all=nan(1,nsample*(itmax+1));


new_logL_topo=nan(1,nsample);
new_logL_thermo=nan(1,nsample);
new_logL_MT = nan(1,nsample);

n = 0;
nn = 1;

for it = 1:itmax+1
    %	Calculate misfit values for each model in the current population.
    
    %	skip forward modelling for starting models if they have been read in from a NAD file
    if (it>1 || noforward~=1)
        hh = waitbar(0,['Running iteration ',num2str(it),' of ', num2str(itmax+1),'...']);
        for i=1:nsample
            misfits=forward_serre(new_models(:,i),param);
            new_misfit(i)= misfits.norm;
            new_logL_topo(i)=misfits.logL_topo;
            new_logL_thermo(i)=misfits.logL_thermo;
            new_logL_MT(i) = misfits.logL_MT;
            waitbar(i/nsample,hh)
        end
        close(hh)
    end
    
    inds = n+1:n+nsample;
    n = n+nsample;
    misfit_all(inds) = new_misfit;
    models_all(:,inds) = new_models;
    logL_topo_all(inds) = new_logL_topo;
    logL_thermo_all(inds) = new_logL_thermo;
    logL_MT_all(inds) = new_logL_MT;
    
    misfit = misfit_all(nn:n);
    models = models_all(:,nn:n);
    
    %	Calculate properties of current misfit distribution.
    %	(Mean,min,best model etc.)
    [mfitmean,mfitminc,mfitmin,mfitord,mopt]=NA_misfits(misfit);
        
    ntot=length(misfit);
    calcmovement=0;
    
    % tranform to scale
    models_sca = nan(size(models));
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
    
    % calculate difference between mest ncells models and either continue
    % inversion or start from a new initial set of samples
    max_d=max(models(:,mfitord(1:ncells))');
    min_d=min(models(:,mfitord(1:ncells))');
    d=max_d-min_d;
   if sum(d<resolution)>0
        misfit=[];
        models=[];
        % Generate or read in starting models
        [new_models,new_misfit]=NA_initial_sample(istype,monte,nsample,nd,range,scales);
        new_models(:,1) = MAP';
        nn = n + 1;
   end 
end

% save relevent outputs
save(['misfit_',ftag,'.mat'],'misfit_all');
save(['models_',ftag,'.mat'],'models_all');
save(['logL_topo_',ftag,'.mat'],'logL_topo_all');
save(['logL_thermo_',ftag,'.mat'],'logL_thermo_all');
save(['logL_MT_',ftag,'.mat'],'logL_MT_all');

% sort the results and pull out the "best" model for plotting
tab = [misfit_all; models_all]';
tab = sortrows(tab);
inds = find(~isnan(tab(:,2)));
tab = tab(inds,:);
param.sample=tab(1,2:end)';

% plot the "best-fit" results
forward_serre_plot(param.sample,param);

% plot the marginal posteriors with the MAP solution as a red dashed line
param_names = {'Time 1 (Ma)','Time 2 (Ma)','Time 3 (Ma)','Uplift rate 1 (mm/yr)',...
    'Uplift rate 2 (mm/yr)','Uplift rate 3 (mm/yr)','Uplift rate 4 (mm/yr)'....
    'K (m^{(1-2m)}/yr)','n','m/n'};
param_dat = tab(:,2:end);
map = param_dat(1,:);
plot_marginals(param_dat,param_names,rangein,map);
