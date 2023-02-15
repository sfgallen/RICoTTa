function [misfit] = forward_serre(sample,param)

%
sample = sample';

%% unpack para structure for forward models
% unpack time parameters
start_t = param.start_t;   % in Myr
dt = param.dt;             % time step in model in yrs

%  unpack MCMC parameters
TI = sample(1:3);
TI = [start_t*1E6 sort(TI,'descend')*1E6 0]; % time range for uplift increments

UI = sample(4:7).*1e-3;     % uplift rates converted from mm/yr to m/yr
K = sample(8);              % erodiblity
n = sample(9);              % slope exponent
m = sample(10).*n;              % drainage area exponent

% unpack steam parameters
S = param.S;                % topotoolbox stream obj
Sd = param.Sd;              % streamwise distance
S_DA = param.S_DA;          % drainage area in m^2
SdA = param.SdA;            % change in drainage area between stream nodes
d = param.d;                % donors
r = param.r;                % recievers
dx = param.dx;              % dx (from fastscape_eroder_prep function)
outlet_nodes = param.outlet_nodes; % stream outlet nodes fastscape_eroder_prep function
obs_elev = param.obs_elev;  % observed elevation
ee = param.ee;              % error on observed elevation

% unpack marine terrace parameters
mt_U = param.mt_u;
mt_UC = param.mt_uc;
mt_age = param.mt_age;

% unpack thermochron parameters
thermo_meas = param.thermo_meas;
T0=param.T0;
lr=param.lr;
TD=param.TD;
rho_c=param.rho_c;
cp_c=param.cp_c;
gg = param.gg;           % initial geothermal gradient

%% river erosion model
% put uplift data into a matrix
U = zeros(length(S_DA),4);
for i = 1:4
    U(:,i) = UI(i);
end
Ui = U(:,1);

% marine terrace model
mod_mt_u = terrace_uplift_forward_model(mt_age,TI(end-1),UI(end-1),UI(end));

% calculate mt misfit
misfit_mt=(sum(log(2*pi/2)+log(mt_UC)+0.5*((mt_U-mod_mt_u)./mt_UC).^2));

% Get the velocity field for the fastscape alg.
A = K.*S_DA.^m; % velocity field

% Calculate Chi for the river network
mn = m/n;

% this model assumes that the river profiles are in steady-state with
% initial uplift field
% calculate initial steady-state river elevations
z0 = calculate_z(Sd,d,r,S_DA,Ui,K,mn,n,obs_elev,outlet_nodes);

% determine time step locations where we change uplift rate and start cosmo
nsteps = round(-diff(TI)/dt);

% % get time loop prepped
steps = 0;
s = 0;

% allocate memory to get thermo variables. These have n-rows = number of
% samples and n-cols = number of timesteps in forward model
topo=zeros(length(thermo_meas.ind),sum(nsteps));
exhumation=zeros(length(thermo_meas.ind),sum(nsteps));

% run the river incision forward model
for l = 1:4
    
    S_U = U(:,l);
    
    for t = 1:nsteps(l)
        
        steps=steps+1;
        s=s+1;
        
        % update uplift field as needed
        z1 = fastscape_eroder_outlets(z0, n, dt, A, d, r, dx, S_U,outlet_nodes);
        
        % calculate and load thermo variables
        topo(:,s)=z1(thermo_meas.ind);   % elevation of the thermo samples through time
        exhumation(:,s)=S_U(thermo_meas.ind)+(z0(thermo_meas.ind)-z1(thermo_meas.ind))./dt; % eroison rate of the thermo samples through time
        
        % update Z0 data with new elevations
        z0 = z1;
        
    end
end

% calculate misfit of river profile
misfit_topo=sum(log(2*pi/2)+log(ee)+0.5*((obs_elev-z1)./ee).^2);

%% Thermochron model(s)
% forward 1D thermal modelling of thermochronological ages
% a finite difference method is used to solve the 1D conductive

nt=length(thermo_meas.x); % number of observations
u=exhumation/(365.25*24*3600); % transform uplift in m/sec
steps=size(u,2); % length of input arrays
time=start_t:-start_t/steps:0; % time vector in Ma
time=time.*(1E6*365.25*24*3600); % time vector in sec
time=time(2:end);
Ts=T0-thermo_meas.elevation*lr/1000;  % surface temperature in C
h=300000/gg;    %<-- Does this represent an initial guess at the rough initial depth of the sample?
hc=h+thermo_meas.elevation;  % crustal thickness in m
nr_nodes=20;
dx=hc/nr_nodes; % spatial resolution in m
hc=dx*nr_nodes;
Tbase=h/1000*gg; % asthenosphere temperature
hpc=9.6E-10*rho_c; % crustal heat production in W/m3
hp_exp=10000; % depth at which hp is 1/e
SC=0.48; % ratio of dt*TD/dx^2
dt_thermo=SC*min(dx)^2/TD;    % determine dt (in seconds)
time_new=floor(time(1)/dt_thermo)*dt_thermo:-dt_thermo:0; % new time vector according to dt in sec
tsteps=length(time_new);
zsteps=hc(1)/dx(1)+1;

% resample topo and u for dt
u_new=zeros(nt,length(time_new));
topo_new=zeros(nt,length(time_new));
for i=1:nt
    u_new(i,:)=interp1(time,u(i,:),time_new);
    topo_new(i,:)=interp1(time,topo(i,:),time_new);
end

% initialize temperature and heat production
T=zeros(tsteps,nt,zsteps);
hp=zeros(1,nt,zsteps);
for i=1:nt
    T(:,i,1)=Ts(i);
    T(:,i,end)=Tbase;
end
for i=1:hc/dx+1
    T(1,:,i)=Ts(:)+i*((Tbase(:)-Ts(:))/(hc/dx+1));
    hp(1,:,i)=hpc*exp(-((i-1)*dx(:))/hp_exp);
end

% calculate stable geotherm
change=1;
t=1;
j=2:round((hc/dx));
while change>0.1
    t=t+1;
    % variable in the lithosphere
    T(t,:,j)=T(t-1,:,j)+SC*(T(t-1,:,j+1)-2*T(t-1,:,j)+T(t-1,:,j-1))+hp(1,:,j)*dt_thermo/(rho_c*cp_c);
    change=sum(abs(T(t-1,1,:)-T(t,1,:)));
end

% time-variable exhumation/burial history
Tvar=zeros(tsteps,nt,zsteps); % variable Geotherm
Tvar(1,:,:)=T(t,:,:); % initial steady-state Geotherm
Tvar(:,:,1)=T(:,:,1);
Tvar(:,:,end)=T(:,:,end);

% calculate total exhumation/burial and original location of samples
u_total=sum(u_new')*dt_thermo; % total exhumation/burial
u_total(u_total<0)=0;
location=u_total; % sample location (change with time)
if tsteps>1E2 % determine step at which the temperature is recorded, min 100
    step=floor(tsteps/1E2);
else
    step=1;
end
tT=zeros(nt,round(tsteps/step)); % Temperature path of sample location
z=zeros(nt,zsteps); % depth vector from sample location to depth of 300? isotherm
for i=1:nt % find temperature at start time
    z(i,:)=0:dx(i):hc(i);
    tT(i,1)=interp1(z(i,:),squeeze(Tvar(1,i,:))',location(i));
end

% model temperature evolution with time and record tT path of samples
j=round(j);
u_matrix=zeros(tsteps,nt,zsteps);
dx_matrix=zeros(tsteps,nt,zsteps);
for o=1:zsteps
    u_matrix(:,:,o)=u_new';
    for i=1:nt
        dx_matrix(:,i,o)=dx(i);
    end
end
i=1;
for t=2:tsteps
    Tvar(t,:,j)=Tvar(t-1,:,j)+SC*(Tvar(t-1,:,j+1)-2*Tvar(t-1,:,j)+Tvar(t-1,:,j-1))+hp(1,:,j)*dt_thermo/(rho_c*cp_c)+(Tvar(t-1,:,j+1)-Tvar(t-1,:,j)).*(u_matrix(t,:,j).*dt_thermo./dx_matrix(t,:,j));
    % find temperature at sample locations
    location=location-(u_new(:,t)'*dt_thermo);
    if mod(t,step)==0
        i=i+1;
        for o=1:nt
            tT(o,i)=interp1(z(o,:),squeeze(Tvar(t,o,:))',location(o));
        end
    end
end
% simplify tT-path
time_i=time_new(1:step:end)/(1E6*365.25*24*3600);
if length(time_i)<length(tT(1,:))
    time_i=[time_i 0];
end
time_i(isnan(sum(tT)))=[];
tT(:,isnan(sum(tT)))=[];
for i=1:nt
    ps = dpsimplify([tT(i,:)' time_i'],0.1);
    thermo{i}.temp=ps(:,1)';
    thermo{i}.time=ps(:,2)';
end

% calculate thermochron ages
for i=1:nt
    
    % calculate ZFT age (Tagami et al. 1998)
    thermo{i}.time(end)=0;
    [modelled_data.zft(i),ftld,ftldmean,ftldsd]=Mad_Zirc(thermo{i}.time,thermo{i}.temp,0,1);
    
    % calculate AFT age and length distribution (Ketcham et al. 2007)
    Dpar=2;
    l0=15.3;
    % reduce TT-path to ~20 nodes from 180?C to surface temperature
    thermo{i}.time=flipdim(thermo{i}.time,2);
    thermo{i}.temp=flipdim(thermo{i}.temp,2);
    nn=find(thermo{i}.temp>180,1,'first');
    dd=1;
    if nn>20
        dd=round(nn*1.5/20);
    else isempty(nn);
        nn=length(thermo{i}.temp);
    end
    [modelled_data.aft(i),modelled_data.aft_length_pdf(i,:),modelled_data.aft_MTL(i)]=ft_dpar_l0(thermo{i}.temp(1:dd:nn),thermo{i}.time(1:dd:nn),Dpar,l0);
    
    % calculate AHe age with RDAAM model (Flowers et al. 2009)
    [modelled_data.ahe(i),~,~]=RDAAM_MEX(thermo{i}.time,thermo{i}.temp,60,20,40,200); % to do, pass SER, U, Th and Sm in from obsverations
    % older model (not recommended), Farley 2000
    %modelled_data.ahe(i)=Mad_He(thermo{i}.time,thermo{i}.temp,length(thermo{i}.temp),1);
    
    % calculate ZHe age with ZRDAAM model (Guenthner et al. 2013)
    [modelled_data.zhe(i),~,~]=ZRDAAM_MEX(thermo{i}.time,thermo{i}.temp,60,400,400,0); % to do, pass SER, U, Th and Sm in from obsverations
    % older model (not recommended), Reiners et al. 2004
    %modelled_data.zhe(i)=Mad_He(thermo{i}.time,thermo{i}.temp,length(thermo{i}.temp),2);
end

% calculate misfit of thermochronological data
misfit_thermo=0;
if isfield(thermo_meas,'ahe')
    misfit_thermo=misfit_thermo+sum(log(2*pi/2)+log(thermo_meas.ahe_error)+0.5*((thermo_meas.ahe-modelled_data.ahe(1:5))./thermo_meas.ahe_error).^2);
end
if isfield(thermo_meas,'aft')
    misfit_thermo=misfit_thermo+sum(log(2*pi/2)+log(thermo_meas.aft_error)+0.5*((thermo_meas.aft-modelled_data.aft(6:8))./thermo_meas.aft_error).^2);
end
if isfield(thermo_meas,'zhe')
    misfit_thermo=misfit_thermo+sum(log(2*pi/2)+log(thermo_meas.zhe_error)+0.5*((thermo_meas.zhe-modelled_data.zhe)./thermo_meas.zhe_error).^2);
end
if isfield(thermo_meas,'zft')
    misfit_thermo=misfit_thermo+sum(log(2*pi/2)+log(thermo_meas.zft_error)+0.5*((thermo_meas.zft-modelled_data.zft)./thermo_meas.zft_error).^2);
end

% if length data provided on AFT, compare.
modelled_data.aft_length_pdf(modelled_data.aft_length_pdf<1E-5)=1E-5;
if isfield(thermo_meas,'aft_length')
    for i=1:length(thermo_meas.aft)
        if length(thermo_meas.aft_length(:,1))>length(thermo_meas.aft(i))
            misfit_thermo=misfit_thermo-sum(log(modelled_data.aft_length_pdf(i+5,(round((thermo_meas.aft_length(i,:)+0.1).*10)))))/201;
        else
            misfit_thermo=misfit_thermo+log(2*pi/2)+log(thermo_meas.aft_length_sd(i))+0.5*((thermo_meas.aft_length(i)-modelled_data.aft_MTL(i+5))./thermo_meas.aft_length_sd(i)).^2;
        end
    end
end

%% pack up all the misfit data into a structure
% individual logL and global normalized logL
misfit.logL_topo=misfit_topo;
misfit.logL_thermo=misfit_thermo;
misfit.logL_MT = misfit_mt;
misfit.norm=misfit_topo/param.nr_topo+misfit_thermo/param.nr_thermo+misfit_mt./length(mt_U);

end