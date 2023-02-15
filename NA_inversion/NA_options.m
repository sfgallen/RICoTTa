%----------------------------------------------------------------------------

%       NA_options - reads in all input options for NA algorithm.

%       Calls no other routines.

%Comments:
%                This routine assumes that direct access
%                files are opened with the record length specified
%                in bytes. This is the default for most machines
%                but not Dec machines. (Often a compiler option is
%                available on the DEC/compaq to use bytes rather than
%                4-byte words.

%					M. Sambridge, Oct. 1996

%----------------------------------------------------------------------------

function	[monte,istype,nsleep,noforward,nclean,nsamplei,nsample,itmax,ncells]=NA_options(nsample_max, nit_max, nmod_max,...
    nsleep_max, nd)

global	lu_na lu_out lu_sum lu_det lu_sob lu_dis lu_nad verbose debug timing summary
global	sobol iseed iproc nproc lroot ndin

%       fprintf header to summary file.
if (lroot) 
    fprintf(lu_sum,'Summary information of Neighbourhood algorithm performance') ;
end

%       Question and answer session to read in NA options.
load na_input.mat

%       Note: logical unit for reading in NA options (lu) may be set to the screen (5) or a file.
monte = 0;
if (itype==1)
    monte = 1;
elseif (itype==0)
    %       shrink = 1;
    % elseif (itype==-1)
    %       shrink = 0;
else
    disp('Error - in NA options file')
    error('NA:NAoptions',' Invalid algorithm type entered ')
end

nsleep = 1;

%	removed option from input
sobol = 0;
if (yesorno=='y' || yesorno=='Y') sobol=1; end

%	read initial sample type
noforward = 0;
if (istype<0)
    istype = -istype;
    noforward = 1;
end

%	read in size of initial sample from nad file
if (istype==1)
    load na_nad.mat % na.nad file with len,nd_nad,nsamplei,nh,nhu
    if (nd_nad~=nd)
        disp('Error detected in reading starting models from NAD file')
        disp('')
        disp(' Dimension of parameter space obtained')
        disp(' from user provided subroutine NA')
        disp(' differs from that read in from NAD file')
        disp('')
        disp(' Perhaps this is the wrong NAD file ?')
        disp(' or user subroutine is in error ?')
        disp('')
        disp(' Remedy: correct input NAD file or')
        disp('         user_init and recompile')
        disp('')
        error('NA:NaNadFile','Wrong na.nad file')
    end
end

ntotal = nsamplei + nsample*itmax;
nclean = 500;
verbose = 0;
summary = 0;
if (infolevel==1) summary=1; end
if (infolevel==2) summary=1; end
if (infolevel==2) verbose=1; end
timing = 0;
if (yesornotiming=='y' || yesornotiming=='Y') timing=1; end
debug = 0;
if (yesornodebug=='y' || yesornodebug=='Y') debug=1; end

%	perform parameter checking
nsam = max([nsample nsamplei]);
if (nsam>nsample_max)
    error('NA:nsample_max','Decrease sample size or increase parameter nsample_max')
elseif (ntotal>nmod_max)
    error('NA:nmod_max','Increase maximum number of iterations and recompile')
elseif (itmax>nit_max)
    error('NA:nit_max','Decrease number of iterations or increase parameter nit_max and recompile')
end
if (ncells>nsample || ncells>nsamplei)
    error('NA:re_sample','Decrease number of re-sampled cells or increase sample size')
elseif (nsleep>nsleep_max)
    error('NA:nsleepmax','Decrease size of input variable nsleep or increase parameter nsleep_max and recompile')
end
if (istype>2 || istype<0)
    error('NA:istype','Only values 0, 1 and -1 currently supported')
end

%       Display NA options
if (lroot)
    fprintf(lu_sum,'---------------------------------------------------------');
    fprintf(lu_sum,'Options set: ');
    if (monte)
        fprintf(lu_sum,' Uniform Monte Carlo Algorithm');
    else
        fprintf(lu_sum,' Neighbourhood Algorithm');
    end
    fprintf(lu_sum,[' Initial sample size           : ',num2str(nsamplei)]);
    fprintf(lu_sum,[' Sample size                   : ',num2str(nsample)]);
    fprintf(lu_sum,[' Total number of iterations    : ',num2str(itmax)]);
    fprintf(lu_sum,[' Number of cells re-sampled    : ',num2str(ncells)]);
    fprintf(lu_sum,[' Total number of models        : ',num2str(ntotal)]);
    fprintf(lu_sum,[' Random seed value             : ',num2str(iseed)]);
    if (sobol)
        fprintf(lu_sum,' SAS Quasi-random sequence used');
    else
        fprintf(lu_sum,' Pseudo-random sequence ');
    end
    if (istype==1)
        fprintf(lu_sum,' Starting models read in from NAD file');
    elseif (istype==2)
        fprintf(lu_sum,' Starting models specified by user');
    else
        fprintf(lu_sum,' Starting models generated randomly ');
    end
    fprintf(lu_sum,'---------------------------------------------------------')
    
    if (summary)
        fprintf(lu_out,'---------------------------------------------------------')
        fprintf(lu_out,'Options set: ');
        if (monte)
            fprintf(lu_out,' Process                       : Uniform Monte Carlo');
        else
            fprintf(lu_out,' Process                       :  Neighbourhood Algorithm');
        end
        fprintf(lu_out,[' Initial sample size           : ',num2str(nsamplei)]);
        fprintf(lu_out,[' Sample size                   : ',num2str(nsample)]);
        fprintf(lu_out,[' Total number of iterations    : ',num2str(itmax)]);
        %       fprintf(lu_out,[' Number of axis samples        : ',num2str(naxis)]);
        %       fprintf(lu_out,' Number of samples per reset   ');
        %       fprintf(lu_out,[' of distance database          : ',num2str(nclean)]);
        fprintf(lu_out,[' Number of cells re-sampled    : ',num2str(ncells)]);
        fprintf(lu_out,[' Total number of models        : ',num2str(ntotal)]);
        %       fprintf(lu_out,[' Length of walk in cell        : ',num2str(nsleep)]);
        fprintf(lu_out,[' Random seed value             : ',num2str(iseed)]);
        if (sobol)
            fprintf(lu_out,' SAS Quasi-random sequence used');
        else
            fprintf(lu_out,' Pseudo-random sequence used');
        end
        if (istype==1)
            fprintf(lu_out,' Starting models read in from NAD file');
        elseif (istype==2)
            fprintf(lu_out,' Starting models specified by user');
        else
            fprintf(lu_out,' Starting models generated randomly ');
        end
        fprintf(lu_out,'---------------------------------------------------------')
    end
end
if (verbose && lroot)
    lu_det=fopen('na.det');
    fprintf(lu_det,'---------------------------------------------------------')
    fprintf(lu_det,' Details of each iteration')
    fprintf(lu_det,'---------------------------------------------------------')
end
end