function [misfit] = forward_sila_plot(sample,param)

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
SdAw = SdA;       % drainage area weighting factor
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

% unpack cosmo-related parameters
t_record = param.t_record*1e6;  % time to start doing cosmo model converted from Myr to yrs
rho_c = param.rho_c;
muon = param.muon;
cosmo = param.cosmo;
cosmo_meas = param.cosmo_meas;

% unpack thermochron parameters
thermo_meas = param.thermo_meas;
T0=param.T0;
lr=param.lr;
TD=param.TD;
rho_c=param.rho_c;
cp_c=param.cp_c;
gg = param.gg;           % initial geothermal gradient

%% put uplift data into a matrix
U = zeros(length(S_DA),4);
for i = 1:4
    U(:,i) = UI(i);
end
Ui = U(:,1);

%% marine terrace model
mod_mt_u = terrace_uplift_forward_model(mt_age,TI(end-1),UI(end-1),UI(end));

% calculate mt misfit
misfit_mt=(sum(log(2*pi/2)+log(mt_UC)+0.5*((mt_U-mod_mt_u)./mt_UC).^2));

figure
errorbar(mt_age,mt_U,mt_UC,'k.'); hold on
plot(mt_age,mt_U,'ko','markerfacecolor',[0.5 0.5 0.5]);
plot(mt_age,mod_mt_u,'ko','markerfacecolor',[0.0 0.0 0.9]);
title('Marine terrace results')
xlabel('Terrace age (kyr)');
ylabel('Total uplift (m)');

%% River model 
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
steps_cosmo=t_record/dt;

% % get time loop prepped
% tsteps = round(run_time/dt);
steps = 0;
s = 0;
o = 0;

% allocate memory to get erosion rates for cosmo
erosion_river = zeros(length(z0),round(t_record/dt));
erosion_catchment = zeros(length(z0),round(t_record/dt));
elevation_river = zeros(length(z0),round(t_record/dt));

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
        
        % if within the cosmo timescale
        if sum(nsteps)-steps<steps_cosmo
            % calculate erosion/incision rate along channel
            o=o+1;
            erosion_river(:,o)=(S_U + (z0-z1)./dt).*1E3;
            erosion_river(outlet_nodes,o)=S_U(outlet_nodes).*1E3;
            erosion_catchment(:,o) = erosion_river(:,o).*SdA;
            elevation_river(:,o) = z1;
        end
        
        % update Z0 data with new elevations
        z0 = z1;
        
    end
end

% calculate misfit of river profile
misfit_topo=sum(log(2*pi/2)+log(ee)+0.5*((obs_elev-z1)./ee).^2);

% plot results
figure
plot(Sd,obs_elev,'.','color',[0.5 0.5 0.5]); hold on
plot(Sd,z1,'b.');
title('River results')
xlabel('Distance (m)');
ylabel('Elevation (m)');
legend('Observed','Modeled')
%% cosmo model
% these are the river channel nodes upstream of the
node=unique(cosmo_meas.ind); % find nodes that need to be calculated
node(node==0)=[];

TCN_channel = zeros(length(z1),1);

cosmo.pressure=1;
load al_be_consts_v22
nuclide=10;

for o=1:length(node)
    % prepare variables
    cosmo=rmfield(cosmo,'pressure');
    elevation=elevation_river(node(o),:)';   % all river elevations for a given node for cosmo time
    erosion=erosion_river(node(o),:)';       % all river erosion rates for a given node for cosmo time
    
    % get the lat and lon for the node
    cosmo.lat=repmat(cosmo.lat_river(node(o)),1,length(elevation));
    cosmo.long=repmat(cosmo.long_river(node(o)),1,length(elevation));
    
    % forward modelling of a nuclide concentration of river cosmo
    % input:    cosmo - time and uplift rates
    % output:   modelled cosmogenic 10Be concentration
    
    
    if (~isfield(cosmo,'pressure'))
        if (strcmp(cosmo.aa,'std'))
            % Old code
            % cosmo.pressure = stdatm(cosmo.elevation);
            % New code
            % Negative longitude check is in NCEPatm_2.m
            cosmo.pressure = NCEPatm_2(cosmo.lat',cosmo.long',elevation);
        elseif (strcmp(cosmo.aa,'ant'))
            cosmo.pressure = antatm(elevation);
        end
    end
    
    % Catch confusion with pressure submission. If cosmo.pressure is already
    % set, it should have a submitted value. If zero, something is wrong.
    if isempty(cosmo.pressure)
        error(['Sample.pressure extant but empty on cosmo ' cosmo.sample_name]);
    elseif (cosmo.pressure == 0)
        error(['Sample.pressure = 0 on cosmo ' cosmo.sample_name]);
    end
    
    % Convert cosmo thickness to g/cm2; get thickness SF
    
    cosmo.thickgcm2 = cosmo.thickness.* rho_c/1E3;
    if cosmo.thickness > 0
        cosmo.thickSF = thickness(cosmo.thickness,al_be_consts.Lsp,rho_c/1E3);
    else
        cosmo.thickSF = 1;
    end
    
    % Negative longitude catch
    
    if cosmo.long < 0
        cosmo.long = cosmo.long + 360;
    end
    
    % Initialize the result flags.
    
    results.flags = [];
    
    % ------------------ 2. NUCLIDE-SPECIFIC ASSIGNMENTS ---------------------
    
    if nuclide == 10
        %N = cosmo.N10;
        %delN = cosmo.delN10;
        mc.Natoms = al_be_consts.Natoms10;
        mc.sigma190 = al_be_consts.sigma190_10;
        mc.k_neg = al_be_consts.k_neg10;
        mc.delsigma190 = al_be_consts.delsigma190_10;
        mc.delk_neg = al_be_consts.delk_neg10;
        l = al_be_consts.l10;
        L = al_be_consts.Lsp;
        P_ref_St = al_be_consts.P10_ref_St; delP_ref_St = al_be_consts.delP10_ref_St;
        P_ref_Du = al_be_consts.P10_ref_Du; delP_ref_Du = al_be_consts.delP10_ref_Du;
        P_ref_De = al_be_consts.P10_ref_De; delP_ref_De = al_be_consts.delP10_ref_De;
        P_ref_Li = al_be_consts.P10_ref_Li; delP_ref_Li = al_be_consts.delP10_ref_Li;
        P_ref_Lm = al_be_consts.P10_ref_Lm; delP_ref_Lm = al_be_consts.delP10_ref_Lm;
        nstring='Be-10';
    elseif nuclide == 26
        N = cosmo.N26;
        delN = cosmo.delN26;
        mc.Natoms = al_be_consts.Natoms26;
        mc.sigma190 = al_be_consts.sigma190_26;
        mc.k_neg = al_be_consts.k_neg26;
        mc.delsigma190 = al_be_consts.delsigma190_26;
        mc.delk_neg = al_be_consts.delk_neg26;
        l = al_be_consts.l26;
        L = al_be_consts.Lsp;
        P_ref_St = al_be_consts.P26_ref_St; delP_ref_St = al_be_consts.delP26_ref_St;
        P_ref_Du = al_be_consts.P26_ref_Du; delP_ref_Du = al_be_consts.delP26_ref_Du;
        P_ref_De = al_be_consts.P26_ref_De; delP_ref_De = al_be_consts.delP26_ref_De;
        P_ref_Li = al_be_consts.P26_ref_Li; delP_ref_Li = al_be_consts.delP26_ref_Li;
        P_ref_Lm = al_be_consts.P26_ref_Lm; delP_ref_Lm = al_be_consts.delP26_ref_Lm;
        nstring='Al-26';
    end
    
    % ---------------------- 3. INITIAL GUESS ---------------------------------
    
    % P_mu_0_diag = P_mu_total(0,mean(cosmo.pressure),mc,'yes');
    %
    % P_temp = (P_ref_St.*stone2000(mean(cosmo.lat),mean(cosmo.pressure),1).*cosmo.thickSF.*cosmo.topocorr)...
    %     + P_mu_0_diag.P_fast + P_mu_0_diag.P_neg;
    
    
    % -------------------- 4. GET THE EROSION RATES -----------------------------
    
    % Make P(t). Go out to 10M. See the hard-copy documentation for this
    % function and get_al_be_age.m for more details.
    
    tv = [0:500:6500 6900 7500:1000:11500 12000:1000:800000 logspace(log10(810000),7,200)];
    temp_M = interp1(al_be_consts.t_M(1:end-1),al_be_consts.M(1:end-1),tv(16:end));
    
    % uncomment for calculation of erosion rates
    %     opts = optimset('fzero');
    %     opts = optimset(opts,'tolx',1e-8,'display','off');
    
    if muon==1
        P_St = stone2000(cosmo.lat(i),cosmo.pressure(i),0.9938).*P_ref_St.*cosmo.topocorr; % scalar
        P_mu = braucher2013(cosmo.pressure(i),0.9938);
    else
        P_St = stone2000(cosmo.lat(i),cosmo.pressure(i),1).*P_ref_St.*cosmo.topocorr; % scalar
        % Precompute P_mu(z) to ~200,000 g/cm2
        % This log-spacing setup for the step size has relative accuracy near
        % 1e-3 at 1000 m/Myr erosion rate.
        % start at the mid-depth of the cosmo.
        z_mu = [0 logspace(0,5.3,100)]+(cosmo.thickgcm2./2);
        P_mu_z = zeros(size(z_mu));
        P_mu_z = P_mu_total(z_mu,cosmo.pressure(i),mc);
        c3.z_mu = z_mu-(cosmo.thickgcm2./2); % take 1/2 depth away again so t will match P
        c3.P_mu_z = P_mu_z;
    end
    
    % Yet another constants block for ET_objective.m
    c3.tv = tv;
    c3.l = l;
    c3.tsf = cosmo.thickSF;
    c3.L = L;
    
    c3.P_sp_t = P_St;
    c3.P_mu = P_mu;
    c3.L_muon=4656; 	% attenuation length scale of muons in g/cm-2 according to Braucher et al. 2013
    
    % ------------------------ 5. Forward modelling of surface CN concentration ------------------
    
    E1=erosion./10*rho_c/1E3; % erosion rate through time in g/cm2*yr
    if E1(1)==0
        E1(1)=1E-4;
    end
    % calculate steady state CN profile
    [N_profile,depth_array]=calc_steady_profile(c3,E1(1),muon);
    % find erosion rate that differs from steady state erosion rate
    I=find(abs(E1(:)-E1(1))>1E-10,1,'first');
    if I>0 & I<length(E1)
        depth_vector=cumsum(flipud(erosion(I:end)).*dt*rho_c/1E4);
        depth_vector=[0;depth_vector(1:end-1)];
        % interpolate depth_vector to dt=100 years
        xi=dt/100;
        depth_vector=interp1(xi:xi:length(depth_vector)*xi,depth_vector,xi:1:length(depth_vector)*xi);
        % find CN of node that end up at the surface
        N=interp1(depth_array,N_profile,depth_vector(end));
        if isnan(N)
            N=min(N_profile);
        end
        
        % interpolate muon production for modelled time vector
        if muon==1
            % calculate nuclide production according to spallogenic und muon reactions for modelled time vector (depth)
            P_sp_z_target=c3.P_sp_t*exp(-depth_vector/c3.L);
            P_m_z_target=c3.P_mu*exp(-depth_vector/c3.L_muon);
            P_total=P_sp_z_target+P_m_z_target;
            for j=1:length(depth_vector)
                N=(N+P_total(end-j+1)*dt/xi)*exp(-c3.l*dt/xi);
            end
        else
            P_mu_z_target=interp1(c3.z_mu,c3.P_mu_z,depth_vector);
            % calculate nuclide production for modelled time vector (depth)
            P_sp_z_target=c3.P_sp_t*exp(-depth_vector/c3.L);
            for j=1:length(depth_vector)
                N=(N+P_mu_z_target(end-j+1)*dt/10+P_sp_z_target(end-j+1)*dt/10)*exp(-c3.l*dt/10);
            end
        end
    else
        N = N_profile(1);
    end
    
    % calculate a corresponding erosion rate
    %         [cosmo.E_mod(i),fval_St,exitflag_St,output] = ...
    %             fzero(@(x) ET_objective_new(x,c3,N(i)),E1(end),opts);
    %         cosmo.E_mod(i)=cosmo.E_mod(i)*1E4/rho_c;
    
    TCN_channel(node(o))=N;
end

% calculate the mean tcn concentration; for catchments the in-situ TCN are multiplied
% with the local erosion rate and the sum of all the products is divided by
% the sum of all local erosion rates
tcn_mod = zeros(1,length(cosmo_meas.lat));
riv_ero_final = erosion_river(:,end);
for i=1:length(cosmo_meas.lat)
    % This horrendous line of code calculates the mean cosmo concentration
    % upstream of a given node by weighting for drainage area between nodes
    % and the physical erosion rate for each node
    tcn_mod(i)= sum(TCN_channel(cosmo_meas.ind(i,cosmo_meas.ind(i,:)>1)).*...
        (SdAw(cosmo_meas.ind(i,cosmo_meas.ind(i,:)>1)).*riv_ero_final(cosmo_meas.ind(i,cosmo_meas.ind(i,:)>1)))./...
        (sum(SdAw(cosmo_meas.ind(i,cosmo_meas.ind(i,:)>1)).*riv_ero_final(cosmo_meas.ind(i,cosmo_meas.ind(i,:)>1)))));
    
end

% calculate misfit TCN
if isreal(tcn_mod)
    misfit_TCN=sum(log(2*pi/2)+log(cosmo_meas.tcn_error)+0.5*((tcn_mod-cosmo_meas.tcn)./cosmo_meas.tcn_error).^2);
else
    misfit_TCN=1E5;
end

% plot results
figure
subplot(1,3,1)
errorbar(cosmo_meas.tcn,cosmo_meas.tcn_error,'k.'); hold on
plot(cosmo_meas.tcn,'ko','markerfacecolor',[0.5 0.5 0.5]);
plot(tcn_mod,'ko','markerfacecolor',[0.1 0.1 0.9]);
title('Cosmogenic radionuclide results')
xlabel('Sample')
ylabel('Concentration (atoms/g)')
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
    % constant hp
    %T(t,j)=T(t-1,j)+SC*(T(t-1,j+1)-2*T(t-1,j)+T(t-1,j-1))+hpc*dt_thermo/(rho_c*cp_c);
    % variable in the lithosphere
    T(t,:,j)=T(t-1,:,j)+SC*(T(t-1,:,j+1)-2*T(t-1,:,j)+T(t-1,:,j-1))+hp(1,:,j)*dt_thermo/(rho_c*cp_c);
    %T(t,j)=T(t-1,j)+SC*(T(t-1,j+1)-2*T(t-1,j)+T(t-1,j-1))+hp(j)*dt_thermo/(rho_c*cp_c)+(T(t-1,j+1)-T(t-1,j))*(u*dt_thermo/dx);
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
    misfit_thermo=misfit_thermo+sum(log(2*pi/2)+log(thermo_meas.ahe_error)+0.5*((thermo_meas.ahe-modelled_data.ahe)./thermo_meas.ahe_error).^2);
end
if isfield(thermo_meas,'aft')
    misfit_thermo=misfit_thermo+sum(log(2*pi/2)+log(thermo_meas.aft_error)+0.5*((thermo_meas.aft-modelled_data.aft(1:7))./thermo_meas.aft_error).^2);
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
            misfit_thermo=misfit_thermo-sum(log(modelled_data.aft_length_pdf(i,(round((thermo_meas.aft_length(i,:)+0.1).*10)))))/201;
        else
            misfit_thermo=misfit_thermo+log(2*pi/2)+log(thermo_meas.aft_length_sd(i))+0.5*((thermo_meas.aft_length(i)-modelled_data.aft_MTL(i))./thermo_meas.aft_length_sd(i)).^2;
        end
    end
end


subplot(1,3,2)
errorbar(thermo_meas.aft,thermo_meas.aft_error,'k.'); hold on
plot(thermo_meas.aft,'ko','markerfacecolor',[0.5 0.5 0.5]);
plot(modelled_data.aft,'ko','markerfacecolor',[0.1 0.1 0.9]);
title('AFT results')
xlabel('Sample')
ylabel('Age (Myr)')

subplot(1,3,3)
errorbar(thermo_meas.aft_length,thermo_meas.aft_length_sd,'k.'); hold on
plot(thermo_meas.aft_length,'ko','markerfacecolor',[0.5 0.5 0.5]);
plot(modelled_data.aft_MTL,'ko','markerfacecolor',[0.1 0.1 0.9]);
title('AFT results')
xlabel('Sample')
ylabel('Mean teack length (\mum)')


%% pack up all the misfit data into a structure
% individual logL and global normalized logL
misfit.logL_topo=misfit_topo;
misfit.logL_TCN=misfit_TCN;
misfit.logL_thermo=misfit_thermo;
misfit.logL_MT = misfit_mt;
misfit.norm=misfit_topo/param.nr_topo+misfit_TCN/param.nr_cosmo+misfit_thermo/param.nr_thermo+misfit_mt./length(mt_U);

end