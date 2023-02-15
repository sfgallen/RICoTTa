function [N_profile,depth_array]=calc_steady_profile(c3,E,muon)

dt=100;
t_target=0:dt:1000000;
depth_array=t_target*E; % array in g/cm2

if muon==1
    P_mu_z_target=c3.P_mu*exp(-depth_array/c3.L_muon);
else
    P_mu_z_target=interp1(c3.z_mu,c3.P_mu_z,depth_array); % interpolate Muon production for modelled time vector
    P_mu_z_target(isnan(P_mu_z_target))=0;
end
P_sp_z_target=c3.P_sp_t*exp(-(depth_array)/c3.L); % calculate nuclide production for modelled time vector (depth)
N_profile=zeros(size(P_mu_z_target));
P_total=P_mu_z_target+P_sp_z_target;
P_total2=[P_total(2:end) 0];

%N_profile=P_total/(1-exp(-c3.l));
N_profile(end)=P_total(end)/(1-exp(-c3.l));
test=N_profile(end);

% more accurate calculation
% for i=1:length(P_total)-1
%     N_profile(end-i)=(N_profile(end-i+1)+P_total(end-i)*dt)*exp(-c3.l*dt);
% end
% fast calculation 
N_profile=fliplr(cumsum(fliplr(P_total)*dt*exp(-c3.l*dt)));


% for i=1:length(P_total)-1
%     i
%     N_profile(end-i)=(N_profile(end-i+1)+P_total(end-i)*dt)*exp(-c3.l*dt);
%     for j=1:i
%        N_profile(end-j)=(N_profile(end-j+1)+P_total(end-j)*dt)*exp(-c3.l*dt);
%     end
% end



% figure
% plot(N_profile,t_target*E/2.7)
% set(gca,'YDir','reverse');
% ylim([0 1000])

% opts = optimset('fzero');
% opts = optimset(opts,'tolx',1e-8,'display','off');
% 
% [E_profile,fval_St,exitflag_St,output] = ...
%     fzero(@(x) ET_objective(x,c3,N_profile(1)),E*10,opts);
% % resulting apparent erosion rate
% E_profile=E_profile*10/2.7