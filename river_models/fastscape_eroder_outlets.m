function S_Z = fastscape_eroder_outlets(S_Z, n, dt, A, d, r, dx, U, outlets)

% simplified from Campforts and Schwanghart TTLEM funerosion_implin.m and
% funerosion_impnlin.m to model incision along STREAMobj only using the
% detachment limited stream power incision model. There is a data
% preperation frunction, fastscape_eroder_data_prep.m that should be run
% before one goes into the time evolution forloop where this function is
% applied. Modification also fixes all outlet elevations
%
% Inputs:
%   S_Z     River network elevations in meters sorted as STREAMobj
%   n       Slope exponent in stream power incision model
%   dt      time step in years
%   A       the "velocity field" of the steam power incision model
%   d       STEAMobj donors
%   r       STREAMobj recievers
%   dx      distance between nodes in STREAMobj
%   U       scalar or vector of uplift rate in meters per year
%   outlets outlet positions along the stream network
%
% Outputs
%   S_Z     Updated river network elevations in meters sorted as STREAMobj
%
% Modified by: Sean F. Gallen

time=dt;
dte = dt;

while time>0
    time=time-dte;
    if time<0
        dte=dte+time;
        time=0;
    end
    
    S_Z=S_Z+dte.*U; % add uplift to elevation field
    if isscalar(U)
        S_Z(outlets) = S_Z(outlets)-dte.*U; %except for the outlets
    else
        S_Z(outlets) = S_Z(outlets)-dte.*U(outlets);
    end
    
    if n == 1
        for j = numel(d):-1:1
            tt      = A(d(j))*dte/(dx(j));
            S_Z(d(j)) = (S_Z(d(j)) + S_Z(r(j))*tt)./(1+tt);
        end
    else
        for j = numel(d):-1:1
            
            tt      = A(d(j))*dte/(dx(j));
            % z_t
            zt      = S_Z(d(j));
            % z_(t+dt) of downstream neighbor
            ztp1d   = S_Z(r(j));
            % dx
            dx_n      = dx(j);
            
            % initial value for finding root
            if ztp1d < zt
                ztp1    = newtonraphson(zt,ztp1d,dx_n,tt,n);
            else
                ztp1    = zt;
            end
            
            if ~isreal(ztp1) || isnan(ztp1)
               % disp('Non real solutions converted to real')
                ztp1=real(ztp1);
            end
            S_Z(d(j))=ztp1;
        end
    end
end
    function ztp1 = newtonraphson(zt,ztp1d,dx,tt,n)
        
        tempz   = zt;
        tol = inf;
        
        while tol > 1e-3
            % iteratively approximated value
            ztp1  =  tempz - (tempz-zt + ...
                (tt*dx) * ...
                ((tempz-ztp1d)./dx)^n) / ...
                (1+n*tt*((tempz-ztp1d)./dx)^(n-1));
            tol   = abs(ztp1-tempz);
            tempz = ztp1;
        end
    end
end