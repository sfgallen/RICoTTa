function stream_data = prep_streams_for_inversion(S,Sz,Sa)

stream_data.S = S;          % Topotoolbox STREAMobj
stream_data.Sd = S.distance; % streamwise distance
stream_data.S_DA = Sa;      % topologically sorted drainage area in m^2
stream_data.Szf = Sz;       % observed elevation of river network


% data useful for fastscape data prep.
[d, r, dx, outlet_nodes] = inversion_data_prep(S);

stream_data.d = d;          % donors
stream_data.r = r;          % recievers
stream_data.dx = dx;        % dx
stream_data.outlet_nodes = outlet_nodes; % outlet nodes
stream_data.obs_elev = Sz;

% calculate dA per stream node for sed flux calculation
SdA = Sa;
%SdA(S.ixc) = SdA(S.ixc) - SdA(S.ix);
% calculating selevation for the entire river network
for lp = numel(r):-1:1
    SdA(r(lp)) = Sa(r(lp))-Sa(d(lp));
end

id2 = streampoi(S,'confluences','logical');
SdA(id2) = S.cellsize^2;

stream_data.SdA = SdA;


end


function [d, r, dx, outlet_nodes] = inversion_data_prep(S)

% simplified from Campforts and Schwanghart TTLEM updateDrainDir.m function
% to handel incision along STREAMobj only
%
% Inputs:
%   S               STREAMobj
%   S_DA            drainage area in meters^2 along the stream network
%
% Outputs:
%   d               donors
%   r               recievers
%   dx              distance between nodes along network
%   outlet_nodes	locations of outlets
%
% Modified by: Sean F. Gallen

% get donor and recievers from STREAMobj
d = S.ix;
r = S.ixc;

% Find outlet nodes in STREAMobj
outlet_nodes = streampoi(S,'outlets','ix');
for i = 1:length(outlet_nodes)
    outlet_nodes(i) = find(S.IXgrid == outlet_nodes(i));
end

% get distance between nodes along profile
Sd = S.distance;
dx = abs(Sd(d) - Sd(r));

end