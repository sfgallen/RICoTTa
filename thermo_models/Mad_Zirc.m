function [fta,ftld,ftldmean,ftldsd]=Mad_Zirc(time_i,temp_i,out_flag,param_flag)

% subroutine Mad_Trax to calculate zircon fission track age from a given thermal history

% in input:
% real*4  time_i(n)  :   the time values (in Myr) in descending order
%                        at which the thermal history is given
%                        (ex: 100,50,20,10,0); the last value
%                        should always be 0; the first value
%                        should be smaller than 1000.
% real*4  temp_i(n)  :   the thermal history in degree Celsius
% integer n          :   the number of time-temperature pairs used
%                        to describe the temperature history
% integer out_flag   :   =0 only calculate fission track age
%                        =1 also calculate track length distribution
%                           and statistics
% NB integer out_flag=1 not implemented yet for Zircon !

% integer param_flag :   =1 uses parameters for alpha-damaged zircon
%                        =2 uses parameters for zero-damage zircon
%                        (parameters from Rahn et al. 2004)

% in output:
% real*4  fta        :   fission track age in Myr
% real*4  ftld(17)   :   normalised track length distribution where
%                        ftld(k) is the percentage of track with
%                        length between k-0.5 and k+0.5 microns
% real*4  ftldmean   :   mean track length in microns
% real*4  ftldsd     :   track length standard deviation in microns

% This subroutine is based on the subroutine "ftmod.pas" provided by
% Peter van der Beek in December 1995. The algorithm is explained in
% Peter's PhD thesis and is based on the work by Lutz and Omar (1991)
% This adaptation of the program for zircon fission-track annealing was
% written by Peter in August/September 2006 and is based on algorithms
% given by Galbraith & Laslett (1997), Tagami et al. (1998) and
% Rahn et al. (2004)

% References:

% van der Beek, P., 1995.Tectonic evolution of continental rifts, PhD Thesis,
%       Faculty of Earth Sicences, Free University, Amsterdam.

% Lutz, T.M. and Omar, G.I., 1991. An inverse method of modeling thermal
%       histories from apatite fission-track data. EPSL, 104, 181-195.

% Galbraith, R. F., and G. M. Laslett (1997), Statistical modelling of thermal
%       annealing of fission tracks in zircon, Chemical Geology, 140, 123-135.

% Tagami, T., et al. (1998), Revised annealing kinetics of fission tracks in 
%        zircon and geological implications, in Advances in Fission-Track Geochronology, 
%        edited by P. Van den haute and F. De Corte, pp. 99-112, Kluwer Academic 
%        Publishers, Dordrecht, Netherlands.

% Rahn, M. K., et al. (2004), A zero-damage model for fission track annealing in zircon,
%        American Mineralogist, 89, 473-484.


r=zeros(1,1000);
prob=zeros(1,101);

if param_flag==1
    % alpha-damaged zircon
    a=-10.77;
    b=2.599E-4;
    c=1.026E-2;
else
    % zero-damage zircon
    a=-11.57;
    b=2.755E-4;
    c=1.075E-2;
end

% unannealed fission track length
xind=10.8;

% mean length of spontaneous tracks in standards
xfct=10.8;

% calculate the number of time steps assuming 1My time step length
% if run time > 100 My, else take 100 time steps
nstep=floor(time_i(1));

%     if (nstep>8000) stop 'Your final time is greater than '//
%    &                        'the age of the Universe...'
%     if (nstep>4500) stop 'Your final time is greater than '//
%    &                        'the age of the Earth...'
%     if (nstep>1000) stop 'Fission track does not work very well '//
%    &                        'for time spans greater than 1Byr...'
if (nstep>2000)
    error('Fission track does not work very well for time spans greater than 1Byr...')
elseif (nstep>100)
    nstep=100;
    time_interval=time_i(1)/100;
else
    time_interval=1;
end
deltat=time_interval*1E6*365.25*24*3600;

% calculate final temperature
tempp=interp1(time_i,temp_i,0)+273;
rp=0.5;

% beginning of time stepping
for i=1:nstep;
    time=i*time_interval;
    % calculate temperature by linear interpolation
    temp=interp1(time_i,temp_i,time)+273;

 
    % calculate mean temperature over the time step
    tempm=(temp+tempp)/2;
    
    % calculate the "equivalent time", teq
    if (i==1)
        teq=0;
    else
        teq=exp((log(1-rp)-a-(c*tempm))/(b*tempm));
    end
    
    
    % check if we are not getting too close to r=0
    % in which case r remains 0 for all following time steps
    
    %        if (param_flag<4)
    %           rcrit=(expos(1./b,a)-a*c0-1.)/a/c1*(1./tempm-c3)-c2/c1
    %        else
    %           rcrit=(expos(1./b,a)-a*c0-1.)/a/c1*(alog(1./tempm)-c3)-c2/c1
    %        end
    
    %          if (alog(teq+deltat).ge.rcrit)
    %            do j=i,nstep
    %            r(j)=0.
    %            enddo
    %          nstep=i
    %          goto 90
    
    %          end
    
    % otherwise calculate reduction in length, r, over the time step, dt
    
    dt=teq+deltat;
    
    gr=a+((b*tempm)*log(dt))+(c*tempm);
    
    r(i)=1-exp(gr);
    % stop calculation if r<0.4
%     if (r(i)<=0.4)
%         nstep=i;
%         break
%     end
    
    % update variables for next time step
    tempp=temp;
    rp=r(i);
    
    %       print*,i,time,temp,r(i)
end


% all reduction factors for all time steps have been calculated
% now estimate the fission track age by simple summation
% (here it helps to use 1Myr time steps)

sumdj=0;

for i=1:nstep
    if (r(i)<=0.4)
        dj=0;
    elseif (r(i)<=0.66)
        dj=2.15*r(i)-0.76;
    else
        dj=r(i);
    end
    sumdj=sumdj+dj;
end

fta=(xind/xfct)*sumdj*time_interval;

% now (if out_flag.ne.0) let's do some statistics

if (out_flag==0)
    
    % first, calculate probability density function using Luts and Omar (1991)
    % method and assuming a Gaussian distribution
    
    sumprob=0;
    
    for j=1:101
        
        rr=(j-1)/100;
        
        if (rr<=0.43)
            h=2.53;
        elseif (rr<=0.67)
            h=5.08-5.93*rr;
        else
            h=1.39-.61*rr;
        end
        
        fr=0;
        
        i=1:nstep;
        x=(rr-r(i))*xind/h;
        fr=sum((0.39894228*exp(-(x.^2/2)))/h);

        
        prob(j)=fr/nstep;
        sumprob=sumprob+prob(j);
        
    end
    
    % now let's rescale the track length distribution, its mean and standard
    % deviation
    
    ftld(11)=100;
    imin=1;
    
    for l=1:10
        
        imax=round(l*100/xind);
        ftld(l)=0;
        
        for i=imin:imax
            ftld(l)=ftld(l)+prob(i);
        end
        
        ftld(l)=(ftld(l)*100./sumprob);
        ftld(11)=ftld(11)-ftld(l);
        
        imin=imax+1;
        
    end
    
    sumftld=0;
    
    for l=1:11
        sumftld=sumftld+ftld(l)*(l-0.5);
    end
    
    ftldmean=sumftld/100;
    
    devftld=0;
    
    for l=1:11
        devftld=devftld+ftld(l)*(l-0.5-ftldmean)^2;
    end
    
    ftldsd=sqrt(devftld/100);

else
    ftld=NaN;
    ftldmean=NaN;
    ftldsd=NaN;
end