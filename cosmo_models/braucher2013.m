function out = braucher2013(P,Fsp);
% profiles from in situ produced cosmogenic nuclides. Nucl. Instr. Meth. Phys. 
% Res. B 294, 484-490. 
MP0_sd=0.004;
out = MP0.*exp((1013.25 - P)./247);