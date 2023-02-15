function out = NCEPatm_2(site_lat,site_lon,site_elv);

% Looks up surface pressure and 1000 mb temp from NCEP reanalysis
% and calculates site atmospheric pressures using these as inputs to the
% standard atmosphere equation. 
%
% Syntax: pressure = NCEPatm_2(site_lat,site_lon,site_elv);
% 
% Requires:
%       site_lat: latitude (DD). Southern hemisphere is negative.
%       site_lon: longitude (DD). Western hemisphere is negative.
%           Tries to deal with 0-360 longitudes gracefully.
%       site_elv: elevation (m).
%
% Returns site pressure in hPa.
%
% Vectorized. Send vectors of equal length.
%
% Note: this must load the data file NCEP2.mat whenever called. 
% Repeated calls to this function will be slow for this reason. 
%
% Also: This function is OK but not great for Antarctica.
% Use antatm.m instead. 
%
% Remember: it is always better to estimate the average pressure at your 
% site using a pressure-altitude relation obtained from nearby station
% data.
%
% Written by Greg Balco -- UW Cosmogenic Nuclide Lab
% balcs@u.washington.edu
% October, 2007
% Part of the CRONUS-Earth online calculators: 
%      http://hess.ess.washington.edu/math
%
% Copyright 2001-2007, University of Washington
% All rights reserved
% Developed in part with funding from the National Science Foundation.
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License, version 2,
% as published by the Free Software Foundation (www.fsf.org).


% correct negative longitudes

sz = find(site_lon < 0);
site_lon(sz) = site_lon(sz) + 360;

% data load

load NCEP2;

% Interpolate sea level pressure and 1000-mb temperature
% from global reanalysis data grids. 

site_slp = interp2(NCEPlon,NCEPlat,meanslp,site_lon,site_lat);
site_T = interp2(NCEPlon,NCEPlat,meant1000,site_lon,site_lat);

site_T_degK = site_T + 273.15;

% More parameters

gmr = -0.03417; % Assorted constants
dtdz = 0.0065; % Lapse rate from standard atmosphere

% Variable lapse rate experiments -- attempts to make lapse rate a 
% physically reasonable function of temperature. No guarantees on the 
% correctness of this part. 
%
% Not used in final version. Probably more physically correct but 
% has a very limited effect on the overall results. Code retained here as
% comments for those who are interested. If you are an atmospheric
% scientist, feel free to suggest a better way to do this. 
% 
% Temp assumed for lapse rate is 15 degrees less than real temp.
% Chosen because in stdatm, std T is 15 and LR is MALR for 0 deg C.
% Not sure why this is the case. Ask an atmospheric scientist. 
%
%Tlr = site_T_degK - 15;
%
% Calculate the saturation vapor pressure
%
%esat = 10.*0.6112.*exp(17.67.*(Tlr-273.15)./((Tlr-273.15)+243.5));
%
% Calculate the mixing ratio and thence the lapse rate
%
%rv = 0.622.*esat./(1013.25-esat);
%Lv = 2.501e6; R = 287; E = 0.62; g = 9.8066; cpd = 1005.7;
%dtdz = g.*(1+(Lv*rv./R./Tlr))./(cpd + (Lv.^2.*rv.*E./R./(Tlr.^2)));

% Calculate site pressure using the site-specific SLP and T1000 with the
% standard atmosphere equation.

out = site_slp .* exp( (gmr./dtdz) .* ( log(site_T_degK) - log(site_T_degK - (site_elv.*dtdz)) ) );
