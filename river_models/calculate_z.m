function z = calculate_z(Sd,d,r,S_DA,S_U,K,mn,n,obs_elev,outlet_nodes)
% function calculate the steady state elevation of the river network
%
% Author: Sean F. Gallen
% Contact: sean.gallen[at]colostate.edu
% Date modified: 09/28/2020

    z = zeros(size(Sd));
    z(outlet_nodes) = obs_elev(outlet_nodes);
    Sa = (S_U./K).^(1/n).*(1./(S_DA)).^mn;       % chi transformation variable
    
    % calculating selevation for the entire river network
    for lp = numel(d):-1:1
        z(d(lp)) = z(r(lp)) + (Sa(r(lp))+(Sa(d(lp))-Sa(r(lp)))/2) *(abs(Sd(r(lp))-Sd(d(lp))));
    end
end