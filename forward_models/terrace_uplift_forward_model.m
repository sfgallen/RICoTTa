function [mod_ter_z] = terrace_uplift_forward_model(age,U_change_time,Ui,Uf)
% this function calculates the total uplift experenced by a terrace since
% the time of terrace formation based on a simple two stage uplift history.
% The first stage of uplift is assumed to be due to spatially uniform
% uplift. The second stage is modeled a broken plate flexure and uplift is
% parameterized as terrace Euclidean distance from the broken segement
% (e.g. a normal fault).
%
% Inputs:
% age - terrace ages in kyr
% f_dist - Euclidean distance from broke segment
% lambda - flexural parameter
% fault_init_time - fault initiation time in years
% Ui - initial block uplift rate in m/yr
% Uf - uplift rate at broken segment in m/yr
%
% Outputs:
% ter_z - total rock uplift of terrace
%
% Author: Sean F. Gallen
% Contact: sean.gallen[at]colostate.edu
% Date Modified: 09/28/2020

    
terrace_age = age.*1e3;
tdif = terrace_age - U_change_time;
tdif(tdif<0) = 0;
Uf_time = terrace_age - tdif;
mod_ter_z = (Uf*Uf_time)+Ui.*tdif;

end
