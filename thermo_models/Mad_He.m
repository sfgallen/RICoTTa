function Apatite_age=Mad_He(time,temperature,ntime,iflag)

% iflag=1 He in Apatite
% iflag=2 He in Zircon
% iflag=3 Ar in K-feldspar
% iflag=4 Ar in Biotite
% iflag=5 Ar in Muscovite
% iflag=6 Ar in Hornblende

%      D0a2=10.**7.7*3600.*24.*365.25e6
%      Ea=36.2e3*4.184

if (iflag==1)
    % He in Ap (Farley 2000)
    D0a2=5.e5*3600.*24.*365.25e6;
    Ea=138.e3;
elseif (iflag==2)
    % He in Zi (Reiners et al. 2004)
    D0a2=4.6e3*3600.*24.*365.25e6;
    Ea=168.e3;
elseif (iflag==3)
    % Ar in Ks
    D0a2=5.6*3600.*24.*365.25e6;
    Ea=120.e3;
elseif (iflag==4)
    % Ar in Bi 500mic
    D0a2=160.*3600.*24.*365.25e6;
    Ea=211.e3;
elseif (iflag==5)
    % Ar in Mu 500mic
    D0a2=13.2*3600.*24.*365.25e6;
    Ea=183.e3;
elseif (iflag==6)
    % Ar in Ho 500mic
    D0a2=24.*3600.*24.*365.25e6;
    Ea=276.e3;
else
    stop 'option not supported in Mad_He'
end

R=8.314;
if time(1)>100
    dt0=1;
else
    dt0=0.1;
end
n=100;
age=zeros(1,n);
diagonal=zeros(1,n);
sup=zeros(1,n-1);
inf=zeros(1,n);
f=zeros(1,n);
fact=ones(1,n);
fact(1)=0.5;
fact(n)=0.5;

for itime=1:ntime-1
    
    dt=dt0;
    %nstep=max(1,int((time(itime)-time(itime+1)+tiny(dt))/dt));
    nstep=max(1,round((time(itime)-time(itime+1))/dt));
    dt=(time(itime)-time(itime+1))/nstep;
    alpha=0.5;
    dr=1./(n-1);
    
    % beta determines the geometry: 2 = spherical
    %                               1 = cylindrical
    beta=2;
    
    temps=temperature(itime);
    tempf=temperature(itime+1);
    
    Da2now=D0a2*exp(-Ea/R/(temps+273));
    
    for istep=1:nstep
        
        fstep=istep/nstep;
        temp=temps+(tempf-temps)*fstep;
        Da2then=Da2now;
        Da2now=D0a2*exp(-Ea/R/(temp+273));
        f1=alpha*dt*Da2now/dr^2;
        f2=(1.-alpha)*dt*Da2then/dr^2;
        
        i=2:n-1;
        diagonal(i)=1+2*f1;
        sup(i)=-f1.*(1+beta./(i-1)./2);
        inf(i)=-f1.*(1-beta./(i-1)./2);
        f(i)=age(i)+f2.*((age(i+1)-2.*age(i)+age(i-1)) +...
                beta.*(age(i+1)-age(i-1))./(i-1)./2)+dt;
%        for i=2:n-1
%            diagonal(i)=1+2*f1;
%            sup(i)=-f1*(1+beta/(i-1)/2);
%            inf(i)=-f1*(1-beta/(i-1)/2);
%            f(i)=age(i)+f2*((age(i+1)-2*age(i)+age(i-1)) +...
%                beta*(age(i+1)-age(i-1))/(i-1)/2)+dt;
%        end
        
        diagonal(1)=1;
        sup(1)=-1;
        f(1)=0;
        diagonal(n)=1;
        inf(n)=0;
        f(n)=0;
        
        % tridiagonal matrix algorithm
%        age = TDMAsolver(inf,diagonal,sup,f);
%        age = thomas(diagonal,sup,inf,f)
        A=diag(diagonal)+diag(sup',1)+diag(inf(2:end),-1);
        age=A\f';


        agei=sum(age'.*fact.*dr.^3.*(0:n-1).*(0:n-1));
%         agei=0;
%         for i=1:n
%             fact=1;
%             if (i==1 || i==n)
%                 fact=0.5;
%             end
%             agei=agei+age(i)*fact*dr^3*(i-1)*(i-1);
%         end
        agei=3*agei;
        age=age';
    end
end
Apatite_age=agei;