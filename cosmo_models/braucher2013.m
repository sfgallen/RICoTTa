function out = braucher2013(P,Fsp);% Calculates the geographic scaling factor for cosmogenic-muon production as % a function of site atmospheric pressure, according to:% Braucher et al., 2013, Determination of muon attenuation lengths in depth 
% profiles from in situ produced cosmogenic nuclides. Nucl. Instr. Meth. Phys. 
% Res. B 294, 484-490. % Syntax: scalingfactor = braucher2013(pressure)%% Units: % pressure in hPa% Elevation can be converted to pressure with the functions% stdatm.m (general use) and antatm.m (Antarctica). % muon production at sea-level in atoms g-1 a-1MP0=0.028;	
MP0_sd=0.004;
out = MP0.*exp((1013.25 - P)./247);