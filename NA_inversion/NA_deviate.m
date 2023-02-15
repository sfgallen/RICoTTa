%-------------------------------------------------------------------------€”

%     NA_deviate - generates a random deviate according to
%                  a given distribution using a 1-D SAS sequence.
%	   or a pseudo-random sequence depending of logical
%	   `sobol'

%     Comments:
%	If sobol = 1;: (Quasi-random number)

%	   This routine generates a random number
%	   between x1 and x2.
%	   The parameter i is the sequence number from
%	   which the quasi random devaite is drawn.

%	If sobol = 0; (Pseduo-random number)

%	   ran3 is called to calculate a deviate which is
%	   scaled to the input boundaries x1,x2.

%	This version is for resample mode and simply generates
%	a deviate between input values x1 and x2.

%        calls NA_sobol

%-------------------------------------------------------------------------”

function deviate=NA_deviate(x1,x2,id)

global sobol sob_seq
% Use SAS sequence
if (sobol)
    ran=sob_seq(id,1);
    sob_seq(id,1)=[];
    deviate = x1 + (x2-x1)*ran;
else
    ran = rand();
    deviate = x1 + (x2-x1)*ran;
    
end

end