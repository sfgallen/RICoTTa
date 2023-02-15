%----------------------------------------------------------------------------

%       NA_sample - generates a new sample of models using
%                   the Neighbourhood algorithm by distributing
%	    nsample new models in ncells cells.

%    Comments:
%	 If xcur is changed between calls  restartNA
%	 must be set to true. logical restartNA must also
%	 be set to true on the first call.

%       Calls are made to various NA_routines.

%					M. Sambridge
%					Last updated Sept. 1999.

%----------------------------------------------------------------------------
function [new_models,xcur,restartNA,dlist]=NA_sample(na_models,ntot,nsample,nd,nsleep,ncells,misfit,mfitord,range,check,xcur,calcmovement,nclean)


global nxsave ndsave ndc nerr ncald nupd cells taxis tup tcd tdev tna tres
global taxis2 ic verbose debug sobol sob_seq inv_method

new_models=zeros(nd,nsample);

% choose initial axis randomly
idnext=1+round(rand()*(nd-1));

% initialize some variables
ic = ic + 1;
info = 0;
if (mod(ic,nclean)==0)
    resetlist = 1;
end
if (mod(ic,nclean)==0 && verbose)
    info = 1;
end
idiff = 0;
ndc = 0;
cells = 0;
nxsave = 0;
ndsave = 0;
nerr = 0;
taxis = 0;
taxis2 = 0;
tup = 0;
tna = 0;
tcd = 0;
tdev = 0;
tres = 0;
cell = 1;
mopt = mfitord(cell);
ind_cellnext = mopt;
ind_celllast = 0;
dsum = 0;
dcount = 0;
%nsampercell=floor(nsample/ncells);
nsampercell=1;

icount = 0;
if (debug)
    disp([' nsample     = ',num2str(nsample)])
    disp([' nsampercell = ',num2str(nsampercell)])
    disp([' ncells      = ',num2str(ncells)])
end

% loop over samples according to defined inversion method
if inv_method==1
    loop=nsample;
elseif inv_method==0
    loop=ncells;
end
for is = 1:loop
    % choose Voronoi cell for sampling
    ind_cell = ind_cellnext;
    icount = icount + 1;
    if (debug) 
        disp([num2str(cell),' cell = ',num2str(ind_cell)])
    end
    if (ind_cell~=ind_celllast)
        % reset walk to chosen model
        [xcur,restartNA]=NA_restart(na_models,nd,ind_cell,xcur,debug);
    end
    
    if (restartNA)
        resetlist = 1;
        restartNA = 0;
    end
    
    % loop over walk steps
    for il = 1:nsleep
        for iw = 1:nd
            % update dlist and nodex for new axis
            
            if (resetlist==0)
                % incremental update
                [dlist,nodex,dminx]=NNupdate_dlist(idnext,id,dlist,na_models,nd,ntot,xcur);
                nupd = nupd + 1;
            else
                % full update
                [dlist,nodex]=NNcalc_dlist(idnext,na_models,nd,ntot,xcur);
                ncald = ncald + 1;
                resetlist = 0;
            end
            
            id = idnext;

            % Calculate intersection of current Voronoi cell with current 1-D axis
            [x1,x2]=NNaxis_intersect(xcur,id,dlist,na_models,nd,ntot,nodex,range(1,id),range(2,id));
            
            
            % Generate new node in Voronoi cell of input point
            xcur(id)=NA_deviate(x1,x2,id);
            
            % check Voronoi boundaries
%             if (check)
%                 for i=1:nd
%                     xdum(i) = xcur(i);
%                 end
%                 xdum(id) = x1;
%                 call findnearest(xdum,na_models,ntot,nd,nnode)
%                 disp([' Nearest node to x1 = ',num2str(nnode)])
%                 xdum(id) = x2;
%                 call findnearest(xdum,na_models,ntot,nd,nnode)
%                 str=[' Nearest node to x2 = ',num2str(nnode)]
%                 fprintf(lu_sum,str);
%                 disp(str)
%             end
            
            % increment axis
            idnext = idnext + 1;
            if (idnext>nd) 
                idnext=1; 
            end
            
        end
    end
    % put new sample in list
    for i=1:nd
        new_models(i,is) = xcur(i);
    end

    % check nearest node
    
%     if (check)
%         call findnearest(xcur,na_models,ntot,nd,nnode)
%         str=[' Nearest node to new model',num2str(is+ntot),' = ',num2str(nnode)];
%         fprintf(lu_sum,str);
%         disp(str)
%         if (nnode~=nodex)
%             form100=' WARNING node no longer in original Voronoi cell: …
%             original cell = %d, new cell = %d, is = %d, iw = %d, id = %d, nodes :';
%             fprintf(lu_sum,form100,nodex,nnode,is+ntot,iw,id);
%         end
%     end
    
    % find distance moved from Vcell node
    if (calcmovement && ind_cell==mopt)
        dist = 0;
        for i=1:nd
            dd = (na_models(i,nodex) - xcur(i));
            dd = dd*dd;
            dist = dist + dd;
        end
        dsum = dsum + sqrt(dist);
        dcount = dcount + 1;
    end
    
    ind_celllast = ind_cell;
    
    if (icount==nsampercell)
        icount = 0;
        if is<=ncells
            cell = cell + 1;
        else
            cell = 1;
        end
        ind_cellnext = mfitord(cell);
    end

    %					fprintf out average distance
    %					of samples from current Vcell
    %					node

    if (info) 
        cells = idiff/(nd*nsample); 
    end
    
end

% random sampling of n_rand samples 
if inv_method==0
    n_rand=nsample-ncells;
    new_models(:,ncells+1:nsample) = rand(nd,n_rand);
end

end
