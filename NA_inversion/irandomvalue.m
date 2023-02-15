%----------------------------------------------------------------------------

%       irandomvalue - generates a random integer between i1 and i2.

%       Comments:
%                Uses either Sobol-Antonov_Saleev quasi-sequence
%                or pseudo-random number generator, depending on the
%                logical variable sobol.

%                Assumes either random sequence has been initialized.

%                                               M. Sambridge, Aug. 1997

%----------------------------------------------------------------------------

function iran=irandomvalue(i1,i2)


global sobol

%                                               Use SAS sequence
if (sobol)
    rval = NA_sobol(1,1,1);
    iran = i1 + round((rval(1)*(i2-i1+1)));
else
    rran = rand();
    iran = i1 + round((rran*(i2-i1+1)));
end

end