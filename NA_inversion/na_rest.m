
      
    
    %----------------------------------------------------------------------------
    
    %       NA_sas_table - uses a pseudo random number generator to
    %	       build the initializing data for the quasi
    %	       random SAS sequence.
    
    %Input:
    %nt			Number of sequences to be generated
    %lu			A vacant logical unit for writing out
    %			SAS initializing data to a file
    %Output:
    %			table written to file sobol.coeff
    %
    %Comments:
    
    %This routine generates initializing data for multiple
    %       Sobol-Antonov-Saleev quasi-random sequences. For each
    %sequence a degree and order of the primitive polynomial
    %are required, and here they are determined by a particular
    %formula (below).
    %
    %The degree, order and initializing random data for each SAS sequence,
    %are written to a file, the free values have been generated
    %randomly (under-constraints) using ran3 from Numerical Recipes.
    %
    %For each degree and order pair (q,p) q initializing integers are
    %required for each sequence, (M1, M2, ..., Mq), where
    %Mi may be any odd integer less than 2**i. So
    %for the i-th term, there are 2**(i-1) possible values. We fprintf,
    %
    %               Nq = 2**(i-1).
    %
    %Since each initializing datum is independent, the total number
    %of possible sequences for degree q is the product,
    
    %               Ntotal =  N1 x N2 x N3 x ... Nq,
    
    %which gives,
    
    %               Ntotal = prod (for i=1,...,q) 2**(i-1),
    
    %               Ntotal = 2**[sum(for i=1,...,q) (i-1))]
    
    %               Ntotal = 2**(q*(q-1)/2)
    
    %Which gives,
    
    %       q           Ntotal      Nq   Number of primitive polynomials (Np)
    %       1                1       1            1
    %       2                2       2            1
    %       3                8       4            2
    %       4               64       8            2
    %       5             1024      16            6
    %       6            32768      32            6
    %       7          2097152      64           18
    %       8        268435456     128           16
    %       9     6.8719476E10     256           48
    %      10     3.5184372E13     512           60
    c
    c
    %Note the number of possible primitive polynomial orders (Np)
    %and their values are defined the degree.
    %All possible values of polynomial order for degrees up to 10 are
    %contained in the array pporder.
    %(A table can also be found on p 302 of Numerical Recipes in
    %Fortran 2d Ed. Press et al 1992)
    %The product of Np and Ntotal is the total number of possible
    %	sequences for that degree.
    
    %When generating a large number of independent sequences using randomly
    %generated initializing data it is prudent to use only higher degrees
    %because for, say degree 4 there are only 64 possible sequences for
    %each of the two polynomial order values, and so
    %if more than 64 are generated some will be duplicates and hence
    %will produce identical (and not independent) sequences.
    
    %Array pporder contains all possible primitive polynomial orders
    %for each degree up to 10.
    
    %The particular formula used here to choose degree and polynomial
    %order for each independent sequence is to cycle through each degree
    %(starting from 5) and take ntotal/10 sequences from that degree.
    %Once the degree is chosen, the polynomial order cycles through
    %its possible values (given by array pporder). The objective here
    %is to minimize the likelihood of repeated trials with different
    %random seeds reproducing the same sequence. Remember there are
    %a finite number of sequences for each degree (see above).
    
    %       Calls ran3 and assumes that this pseudo random number generator
    %	has been initialized.
    
    %					M. Sambridge, Aug. 1999
    
    %----------------------------------------------------------------------------
    
    function NA_sas_table(nt,lu)
    
    parameter(maxdeg=10)
    
    a=zeros(1,maxdeg);
    ntot=zeros(1,10);
    nprim=zeros(1,10);
    pporder=zeros(60,10);
    
    
    global iproc,nproc,lroot
    
    data		ntot/1.,2.,8.,64.,1024.,32768.,2097152.,
    &                       268435456.,6.8719476E10,3.5184372E13/
    data		nprim/1,1,2,2,6,6,18,16,48,60/
    data		pporder/0,59*0,
    &				1,59*0,
    &				1,2,58*0,
    &				1,4,58*0,
    &				2,4,7,11,13,14,54*0,
    &				1,13,16,19,22,25,54*0,
    &				1,4,7,8,14,19,21,28,31,32,
    &				37,41,42,50,55,56,59,62,42*0,
    &				14,21,22,38,47,49,50,52,56,67,
    &				70,84,97,103,115,122,44*0,
    &                          8,13,16,22,25,44,47,52,55,59,
    &				62,67,74,81,82,87,91,94,103,104,
    &				109,122,124,137,138,143,145,152,
    &				157,167,173,176,181,182,185,191,
    &				194,199,218,220,227,229,230,234,
    &				236,241,244,253,12*0,
    & 				4, 13, 19, 22, 50, 55, 64, 69,
    &				98, 107, 115, 121, 127, 134, 140,
    &				145, 152, 158, 161, 171, 181, 194,
    &				199, 203, 208, 227, 242, 251, 253,
    &				265, 266, 274, 283, 289, 295, 301,
    &				316, 319, 324, 346, 352, 361, 367,
    &				382, 395, 398, 400, 412, 419, 422,
    &				426, 428, 433, 446, 454, 457, 472,
    &				493, 505, 508/
    
    global sobol,iseed
    
    a(1) = 1;
    ii = 0;
    for mdeg=5:10
        if (mdeg>maxdeg) error('NA:degree','too many degrees') end
        numt = floor(ntot(mdeg)/100);
        for num=1:numt
            is = mod(num,nprim(mdeg))+1;
            ip = pporder(is,mdeg);
            ii = ii + 1;
            for j=2:mdeg
                rval = ran3(iseed);
                m = 2**(j-1);
                a(j) = 2*(1+floor(m*rval))-1;
            end
            if (lroot)
                save sobol_coeff.mat nt, iseed, mdeg, ip, a
            end
            if (ii==nt) go to 99 end
            end
        end
        99    continue
        
        fclose(lu)
    end
    
    
    
    
    
    
    
    %—————————————————————————————————————
    
    %     NA_gen_deviates - generates uniform quasi random deviates using
    %                       a Multi-dimensional Sobol-Antonov-Saleev sequence.
    
    %       Input:
    %      n 		: Number of independent sequences to
    %			: be generated (one per dimension)
    %      m 		: Number of terms in each sequence
    %			: to be generated.
    %       Output:
    %      dev(n,m))		: Array of deviates.
    %
    %     Comments:
    %	Generates all deviates required for one execution of
    %	routine NA_sample and stores them in array dev. This is
    %        necessary because of the order the deviates are used
    %	in NA_sample is different from that in which they
    %	can be generated. In this way it is possible to use
    %        independent multi-dimensional SAS sequences in each
    %	resampled cell. Without this `in advance' approach the
    %	same SAS sequence would need to be used across more
    %	than one cell which introduces undesirable cyclicity in
    %	the use of the samples from a single sequence.
    
    %—————————————————————————————————————
    
        function NA_gen_deviates(nd,nc,nsp,nsl,nsample,dev)
            
            global sobol,iseed
            
            nrem = mod(nsample,nc);
            
            n = nd*nc;
            m = nsp*nsl;
            n1 = nd*nrem;
            n2 = nd*(nc-nrem);
            m1 = nsp*nsl;
            m2 = (nsp-1)*nsl;
            
            if (sobol)
                
                if (nrem==0)
                    for i=1:m
                        for j=1:n
                            call NA_sobol(ldummy,j,dev(j,i),1,0)
                        end
                    end
                else
                    for i=1:m1
                        for j=1:n1
                            call NA_sobol(ldummy,j,dev(j,i),1,0)
                        end
                    end
                    for i=1:m2
                        for j=n1+1:n
                            call NA_sobol(ldummy,j,dev(j,i),1,0)
                        end
                    end
                end
            end
            
        end
    
    
                        
                    
                
            
            c
            c-----------------------------------------------------------------------
            c
            %Subroutine NNupdate_dlist - calculates square of distance from
            %			     all base points to new axis, assuming
            %                                    dlist contains square of all distances
            %			     to previous axis dimlast. It also
            %			     updates the nearest node to the
            %			     point x through which the axes pass.
            %
            c
            c-----------------------------------------------------------------------
            c
            Subroutine NNupdate_dlist
            &             (dim,dimlast,dlist,bp,nd,nb,x,node,dmin)
            
            real*4		bp(nd,*)
            real*4		x(nd)
            real*4		dlist(*)
            integer		dim,dimlast
            
            d1 = (x(dimlast)-bp(dimlast,1))
            d1 = d1*d1
            dmin = dlist(1)+d1
            node = 1
            d2 = (x(dim)-bp(dim,1))
            d2 = d2*d2
            dlist(1) = dmin-d2
            for i=2,nb
                d1 = (x(dimlast)-bp(dimlast,i))
                ds = d1
                d1 = dlist(i)+d1*d1
                if (dmin>d1)
                    dmin = d1
                    node = i
                end
                d2 = (x(dim)-bp(dim,i))
                d2 = d2*d2
                dlist(i) = d1-d2
                %          if (i==nb)
                %             disp(' NNupdate_dlist: node 20'
                %             disp(' dlist ',dlist(i)
                %             disp(' dim ',dim
                %             disp(' dimlast ',dimlast
                %             disp(' x ',x(1),x(2)
                %             disp(' bp ',bp(1,20),bp(2,20)
                %             disp(' d1 ',ds
                %             disp(' d2 ',d2
                %          end
                
            end
            
            return
        end
    c-----------------------------------------------------------------------
    c
    %Numerical recipes routine
    c
    c-----------------------------------------------------------------------
    c
    SUBROUTINE indexx(n,arr,indx)
    INTEGER n,indx(n),M,NSTACK
    REAL arr(n)
    PARAMETER (M=7,NSTACK=50)
    INTEGER i,indxt,ir,itemp,j,jstack,k,l,istack(NSTACK)
    REAL a
    for 11 j=1,n
        indx(j)=j
        11    continue
        jstack=0
        l=1
        ir=n
        1     if (ir-l<M)
            for 13 j=l+1,ir
                indxt=indx(j)
                a=arr(indxt)
                for 12 i=j-1,1,-1
                    if (arr(indx(i)).le.a)goto 2
                        indx(i+1)=indx(i)
                        12        continue
                        i=0
                        2         indx(i+1)=indxt
                        13      continue
                        if (jstack==0)return
                            ir=istack(jstack)
                            l=istack(jstack-1)
                            jstack=jstack-2
                        else
                            k=(l+ir)/2
                            itemp=indx(k)
                            indx(k)=indx(l+1)
                            indx(l+1)=itemp
                            if (arr(indx(l+1))>arr(indx(ir)))
                                itemp=indx(l+1)
                                indx(l+1)=indx(ir)
                                indx(ir)=itemp
                                endif
                                if (arr(indx(l))>arr(indx(ir)))
                                    itemp=indx(l)
                                    indx(l)=indx(ir)
                                    indx(ir)=itemp
                                    endif
                                    if (arr(indx(l+1))>arr(indx(l)))
                                        itemp=indx(l+1)
                                        indx(l+1)=indx(l)
                                        indx(l)=itemp
                                        endif
                                        i=l+1
                                        j=ir
                                        indxt=indx(l)
                                        a=arr(indxt)
                                        3       continue
                                        i=i+1
                                        if (arr(indx(i))<a)goto 3
                                            4       continue
                                            j=j-1
                                            if (arr(indx(j))>a)goto 4
                                                if (j<i)goto 5
                                                    itemp=indx(i)
                                                    indx(i)=indx(j)
                                                    indx(j)=itemp
                                                    goto 3
                                                    5       indx(l)=indx(j)
                                                    indx(j)=indxt
                                                    jstack=jstack+2
                                                    if (jstack>NSTACK)stop 'NSTACK too small in indexx'
                                                        if (ir-i+1.ge.j-l)
                                                            istack(jstack)=ir
                                                            istack(jstack-1)=i
                                                            ir=j-1
                                                        else
                                                            istack(jstack)=j-1
                                                            istack(jstack-1)=l
                                                            l=i
                                                            endif
                                                            endif
                                                            goto 1
                                                            END
                                                            c
                                                            %----------------------------------------------------------------------------
                                                            %
                                                            %Numerical Recipes random number generator (used by NA_random)
                                                            c
                                                            %----------------------------------------------------------------------------
                                                            FUNCTION ran3(idum)
                                                            INTEGER idum
                                                            INTEGER MBIG,MSEED,MZ
                                                            %     REAL MBIG,MSEED,MZ
                                                            REAL ran3,FAC
                                                            PARAMETER (MBIG=1000000000,MSEED=161803398,MZ=0,FAC=1./MBIG)
                                                            %     PARAMETER (MBIG=4000000.,MSEED=1618033.,MZ=0.,FAC=1./MBIG)
                                                            INTEGER i,iff,ii,inext,inextp,k
                                                            INTEGER mj,mk,ma(55)
                                                            %     REAL mj,mk,ma(55)
                                                            SAVE iff,inext,inextp,ma
                                                            DATA iff /0/
                                                            %MOD Jean 31/10/08
                                                            call random_number (ran3)
                                                            return
                                                            %end MOD Jean
                                                            if (idum<0 || iff==0)
                                                                iff=1
                                                                mj=MSEED-iabs(idum)
                                                                mj=mod(mj,MBIG)
                                                                ma(55)=mj
                                                                mk=1
                                                                for 11 i=1,54
                                                                    ii=mod(21*i,55)
                                                                    ma(ii)=mk
                                                                    mk=mj-mk
                                                                    if (mk<MZ)mk=mk+MBIG
                                                                        mj=ma(ii)
                                                                        11      continue
                                                                        for 13 k=1,4
                                                                            for 12 i=1,55
                                                                                ma(i)=ma(i)-ma(1+mod(i+30,55))
                                                                                if (ma(i)<MZ)ma(i)=ma(i)+MBIG
                                                                                    12        continue
                                                                                    13      continue
                                                                                    inext=0
                                                                                    inextp=31
                                                                                    idum=1
                                                                                    endif
                                                                                    inext=inext+1
                                                                                    if (inext==56)inext=1
                                                                                        inextp=inextp+1
                                                                                        if (inextp==56)inextp=1
                                                                                            mj=ma(inext)-ma(inextp)
                                                                                            if (mj<MZ)mj=mj+MBIG
                                                                                                ma(inext)=mj
                                                                                                ran3=mj*FAC
                                                                                                return
                                                                                                END
                                                                                                c
                                                                                                c-----------------------------------------------------------------------
                                                                                                c
                                                                                                %Numerical recipes routine
                                                                                                c
                                                                                                c-----------------------------------------------------------------------
                                                                                                c
                                                                                                SUBROUTINE sobseq(n,x)
                                                                                                INTEGER n,MAXBIT,MAXDIM
                                                                                                REAL x(*)
                                                                                                PARAMETER (MAXBIT=30,MAXDIM=6)
                                                                                                INTEGER i,im,in,ipp,j,k,l,ip(MAXDIM),iu(MAXDIM,MAXBIT),iv(MAXBIT*
                                                                                                *MAXDIM),ix(MAXDIM),mdeg(MAXDIM)
                                                                                                REAL fac
                                                                                                SAVE ip,mdeg,ix,iv,in,fac
                                                                                                EQUIVALENCE (iv,iu)
                                                                                                DATA ip /0,1,1,2,1,4/, mdeg /1,2,3,3,4,4/, ix /6*0/
                                                                                                DATA iv /6*1,3,1,3,3,1,1,5,7,7,3,3,5,15,11,5,15,13,9,156*0/
                                                                                                if (n<0)
                                                                                                    for 14 k=1,MAXDIM
                                                                                                        for 11 j=1,mdeg(k)
                                                                                                            iu(k,j)=iu(k,j)*2**(MAXBIT-j)
                                                                                                            11        continue
                                                                                                            for 13 j=mdeg(k)+1,MAXBIT
                                                                                                                ipp=ip(k)
                                                                                                                i=iu(k,j-mdeg(k))
                                                                                                                i=ieor(i,i/2**mdeg(k))
                                                                                                                for 12 l=mdeg(k)-1,1,-1
                                                                                                                    if (iand(ipp,1)~=0)i=ieor(i,iu(k,j-l))
                                                                                                                        ipp=ipp/2
                                                                                                                        12          continue
                                                                                                                        iu(k,j)=i
                                                                                                                        13        continue
                                                                                                                        14      continue
                                                                                                                        fac=1./2.**MAXBIT
                                                                                                                        in=0
                                                                                                                    else
                                                                                                                        im=in
                                                                                                                        for 15 j=1,MAXBIT
                                                                                                                            if (iand(im,1)==0)goto 1
                                                                                                                                im=im/2
                                                                                                                                15      continue
                                                                                                                                stop 'MAXBIT too small in sobseq'
                                                                                                                                1       im=(j-1)*MAXDIM
                                                                                                                                for 16 k=1,min(n,MAXDIM)
                                                                                                                                    ix(k)=ieor(ix(k),iv(im+k))
                                                                                                                                    x(k)=ix(k)*fac
                                                                                                                                    16      continue
                                                                                                                                    in=in+1
                                                                                                                                    endif
                                                                                                                                    return
                                                                                                                                    END
                                                                                                                                    c
                                                                                                                                    c-----------------------------------------------------------------------
                                                                                                                                    c
                                                                                                                                    %Numerical recipes routine adapted to give ind and iselect
                                                                                                                                    c
                                                                                                                                    c-----------------------------------------------------------------------
                                                                                                                                    c
                                                                                                                                    FUNCTION select(k,n,arr,ind,iselect)
                                                                                                                                    INTEGER k,n
                                                                                                                                    REAL select,arr(n)
                                                                                                                                    integer ind(n)
                                                                                                                                    INTEGER i,ir,j,l,mid
                                                                                                                                    REAL a,temp
                                                                                                                                    l=1
                                                                                                                                    ir=n
                                                                                                                                    1     if (ir-l.le.1)
                                                                                                                                        if (ir-l==1)
                                                                                                                                            if (arr(ir)<arr(l))
                                                                                                                                                temp=arr(l)
                                                                                                                                                arr(l)=arr(ir)
                                                                                                                                                arr(ir)=temp
                                                                                                                                                itemp=ind(l)
                                                                                                                                                ind(l)=ind(ir)
                                                                                                                                                ind(ir)=itemp
                                                                                                                                                endif
                                                                                                                                                endif
                                                                                                                                                select=arr(k)
                                                                                                                                                iselect=ind(k)
                                                                                                                                                return
                                                                                                                                            else
                                                                                                                                                mid=(l+ir)/2
                                                                                                                                                temp=arr(mid)
                                                                                                                                                arr(mid)=arr(l+1)
                                                                                                                                                arr(l+1)=temp
                                                                                                                                                itemp=ind(mid)
                                                                                                                                                ind(mid)=ind(l+1)
                                                                                                                                                ind(l+1)=itemp
                                                                                                                                                if (arr(l+1)>arr(ir))
                                                                                                                                                    temp=arr(l+1)
                                                                                                                                                    arr(l+1)=arr(ir)
                                                                                                                                                    arr(ir)=temp
                                                                                                                                                    itemp=ind(l+1)
                                                                                                                                                    ind(l+1)=ind(ir)
                                                                                                                                                    ind(ir)=itemp
                                                                                                                                                    endif
                                                                                                                                                    if (arr(l)>arr(ir))
                                                                                                                                                        temp=arr(l)
                                                                                                                                                        arr(l)=arr(ir)
                                                                                                                                                        arr(ir)=temp
                                                                                                                                                        itemp=ind(l)
                                                                                                                                                        ind(l)=ind(ir)
                                                                                                                                                        ind(ir)=itemp
                                                                                                                                                        endif
                                                                                                                                                        if (arr(l+1)>arr(l))
                                                                                                                                                            temp=arr(l+1)
                                                                                                                                                            arr(l+1)=arr(l)
                                                                                                                                                            arr(l)=temp
                                                                                                                                                            itemp=ind(l+1)
                                                                                                                                                            ind(l+1)=ind(l)
                                                                                                                                                            ind(l)=itemp
                                                                                                                                                            endif
                                                                                                                                                            i=l+1
                                                                                                                                                            j=ir
                                                                                                                                                            a=arr(l)
                                                                                                                                                            ia=ind(l)
                                                                                                                                                            3       continue
                                                                                                                                                            i=i+1
                                                                                                                                                            if (arr(i)<a)goto 3
                                                                                                                                                                4       continue
                                                                                                                                                                j=j-1
                                                                                                                                                                if (arr(j)>a)goto 4
                                                                                                                                                                    if (j<i)goto 5
                                                                                                                                                                        temp=arr(i)
                                                                                                                                                                        arr(i)=arr(j)
                                                                                                                                                                        arr(j)=temp
                                                                                                                                                                        itemp=ind(i)
                                                                                                                                                                        ind(i)=ind(j)
                                                                                                                                                                        ind(j)=itemp
                                                                                                                                                                        goto 3
                                                                                                                                                                        5       arr(l)=arr(j)
                                                                                                                                                                        arr(j)=a
                                                                                                                                                                        ind(l)=ind(j)
                                                                                                                                                                        ind(j)=ia
                                                                                                                                                                        if (j.ge.k)ir=j-1
                                                                                                                                                                            if (j.le.k)l=i
                                                                                                                                                                                endif
                                                                                                                                                                                goto 1
                                                                                                                                                                                END
                                                                                                                                                                                c
                                                                                                                                                                                %----------------------------------------------------------------------------
                                                                                                                                                                                c
                                                                                                                                                                                %       findnearest - finds nearest model to input point
                                                                                                                                                                                c
                                                                                                                                                                                %       Calls no other routines.
                                                                                                                                                                                c
                                                                                                                                                                                %                                               M. Sambridge, Oct. 1996
                                                                                                                                                                                c
                                                                                                                                                                                %----------------------------------------------------------------------------
                                                                                                                                                                                c
                                                                                                                                                                                Subroutine findnearest(smodelt,na_models,ntot,moddim,mopts)
                                                                                                                                                                                
                                                                                                                                                                                real		na_models(moddim,*)
                                                                                                                                                                                real		smodelt(*)
                                                                                                                                                                                
                                                                                                                                                                                common		/findnearestcom/dmin
                                                                                                                                                                                
                                                                                                                                                                                dmin = 0
                                                                                                                                                                                for i=1:moddim
                                                                                                                                                                                    b = na_models(i,1)-smodelt(i)
                                                                                                                                                                                    dmin = dmin + b*b
                                                                                                                                                                                end
                                                                                                                                                                                mopts = 1
                                                                                                                                                                                for j=2,ntot
                                                                                                                                                                                    d = 0
                                                                                                                                                                                    for i=1:moddim
                                                                                                                                                                                        b = na_models(i,j)-smodelt(i)
                                                                                                                                                                                        d = d + b*b
                                                                                                                                                                                    end
                                                                                                                                                                                    if (d<dmin)
                                                                                                                                                                                        dmin = d
                                                                                                                                                                                        mopts = j
                                                                                                                                                                                    end
                                                                                                                                                                                end
                                                                                                                                                                                
                                                                                                                                                                                return
                                                                                                                                                                            end
                                                                                                                                                                            
                                                                                                                                                                        
                                                                                                                                                                        
                                                                                                                                                                        c
                                                                                                                                                                        
                                                                                                                                                                        %---------------------------------------------------------------------
                                                                                                                                                                        c
                                                                                                                                                                        %       fprintf_nad - fprintf a direct access file in
                                                                                                                                                                        %                   multi-record NAD format
                                                                                                                                                                        c
                                                                                                                                                                        %       Input:
                                                                                                                                                                        %             lu                : logical unit of file
                                                                                                                                                                        %             fnme              : filename
                                                                                                                                                                        %             nhmax             : maximum size of array header
                                                                                                                                                                        %             ndmax             : maximum size of array data
                                                                                                                                                                        %             nemax             : maximum size of array models
                                                                                                                                                                        %             iform             : =0  single record format
                                                                                                                                                                        %                                 =1  multi-record format (for large nd)
                                                                                                                                                                        c
                                                                                                                                                                        %       Output:
                                                                                                                                                                        %             nh                : length in bytes of file header
                                                                                                                                                                        %             nhu               : length in bytes of user portion of header
                                                                                                                                                                        %             nd                : dimension of parameter space
                                                                                                                                                                        %             ne                : number of models in ensemble
                                                                                                                                                                        %             header            : header character string of length nh (char)
                                                                                                                                                                        %             data(nd)          : array of data values for each model (real*4)
                                                                                                                                                                        %             models(nd,ne)     : array of model values  (real*4)
                                                                                                                                                                        c
                                                                                                                                                                        %       Comments:
                                                                                                                                                                        %                The direct access NAD file format:
                                                                                                                                                                        c
                                                                                                                                                                        %                VARIABLE       TYPE            SIZE IN BYTES
                                                                                                                                                                        %                nd             int             4
                                                                                                                                                                        %                ne             int             4
                                                                                                                                                                        %                nh             int             4
                                                                                                                                                                        %                nhu            int             4
                                                                                                                                                                        %                header         character       nh
                                                                                                                                                                        %                models         real*4          4*nd*ne
                                                                                                                                                                        %                data           real*4          4*nd
                                                                                                                                                                        %                tail           real*4          4
                                                                                                                                                                        c
                                                                                                                                                                        %                In single record mode a direct access file of length
                                                                                                                                                                        %                [4x(4+nd*ne+ne+1) + nh] bytes is produced.
                                                                                                                                                                        c
                                                                                                                                                                        %                In multi record mode a direct access file of length
                                                                                                                                                                        %                [(ne+1)*(max(20+nh,4(nd+1))] bytes is produced.
                                                                                                                                                                        c
                                                                                                                                                                        %               Calls are made to subroutine read_da.
                                                                                                                                                                        c
                                                                                                                                                                        %                This routine assumes that direct access
                                                                                                                                                                        %                files are opened with the record length specified
                                                                                                                                                                        %                in bytes. This is the default for most machines
                                                                                                                                                                        %                but not Dec machines. (Often a compiler option is
                                                                                                                                                                        %                available on the DEC/compaq to use bytes rather than
                                                                                                                                                                        %                4-byte words.
                                                                                                                                                                        c
                                                                                                                                                                        %                                       M. Sambridge, RSES, November 2001
                                                                                                                                                                        c
                                                                                                                                                                        %---------------------------------------------------------------------
                                                                                                                                                                        c
                                                                                                                                                                        Subroutine fprintf_nad
                                                                                                                                                                        &           (lu,fnme,nd,ne,nh,nhu,header,iform,models,data)
                                                                                                                                                                        c
                                                                                                                                                                        real*4            models(nd,ne)
                                                                                                                                                                        real*4            data(ne)
                                                                                                                                                                        character         header(nh)
                                                                                                                                                                        real*4            tail
                                                                                                                                                                        character*256     fnme
                                                                                                                                                                        logical           warn
                                                                                                                                                                        
                                                                                                                                                                        warn = 1;
                                                                                                                                                                        warn = 0;
                                                                                                                                                                        
                                                                                                                                                                        %                                               fprintf new
                                                                                                                                                                        %                                               multi-record format
                                                                                                                                                                        if (iform==1)
                                                                                                                                                                            
                                                                                                                                                                            %                                               calculate length of header
                                                                                                                                                                            len1 = 4*5+nh
                                                                                                                                                                            len2 = 4*(nd+1)
                                                                                                                                                                            mul  = 1 + (len1-1)/len2
                                                                                                                                                                            %        disp(mul
                                                                                                                                                                            lenh = mul*len2
                                                                                                                                                                            num = ne + mul
                                                                                                                                                                            is1 = num*len2
                                                                                                                                                                            is2 = 4*(5+nd*ne+ne)+nh
                                                                                                                                                                            
                                                                                                                                                                            %        disp(' Number of models                         :',ne
                                                                                                                                                                            %        disp(' Number of dimensions                     :',nd
                                                                                                                                                                            %        disp(' Original header length in bytes          :',len1
                                                                                                                                                                            %        disp(' Final header length in bytes             :',lenh
                                                                                                                                                                            %        disp(' Direct access file record length         :',len2
                                                                                                                                                                            %        disp(' Number of records                        :',num
                                                                                                                                                                            %        disp(' Size of nad file in multi-record format  :',is1
                                                                                                                                                                            %        disp(' Size of nad file in single record format :',is2
                                                                                                                                                                            
                                                                                                                                                                            %                                                       fprintf header
                                                                                                                                                                            open(lu,file=fnme,status='unknown',
                                                                                                                                                                            &          form='unformatted',access='direct',recl=lenh)
                                                                                                                                                                            
                                                                                                                                                                            %                                               fprintf out header
                                                                                                                                                                            %                                               for multi-record format
                                                                                                                                                                            
                                                                                                                                                                            fprintf(lu,rec=1)-mul,nd,ne,nh,nhu,header
                                                                                                                                                                            
                                                                                                                                                                            close(lu)
                                                                                                                                                                            %                                                       fprintf models
                                                                                                                                                                            open(lu,file=fnme,status='unknown',
                                                                                                                                                                            &          form='unformatted',access='direct',recl=len2)
                                                                                                                                                                            
                                                                                                                                                                            for i=1:ne
                                                                                                                                                                                call wnad(lu,mul+i,nd,models(1,i),data(i))
                                                                                                                                                                            end
                                                                                                                                                                            close(lu)
                                                                                                                                                                            
                                                                                                                                                                            
                                                                                                                                                                        else
                                                                                                                                                                            %                                               fprintf original
                                                                                                                                                                            %                                               single record format
                                                                                                                                                                            c
                                                                                                                                                                            %                                               set total length of nad file
                                                                                                                                                                            len = 4+nd*ne+ne+1
                                                                                                                                                                            len = 4*len+nh
                                                                                                                                                                            
                                                                                                                                                                            disp(' size of nad file = ',len
                                                                                                                                                                            %                                               open direct access nad file
                                                                                                                                                                            
                                                                                                                                                                            open(lu,file=fnme,status='unknown',
                                                                                                                                                                            &          form='unformatted',access='direct',recl=len)
                                                                                                                                                                            
                                                                                                                                                                            tail = -999.0
                                                                                                                                                                            fprintf(lu,rec=1)nd,ne,nh,nhu,header,models,data,tail
                                                                                                                                                                            
                                                                                                                                                                            close(lu)
                                                                                                                                                                            
                                                                                                                                                                        end
                                                                                                                                                                        
                                                                                                                                                                        return
                                                                                                                                                                    end
                                                                                                                                                                    
                                                                                                                                                                    Subroutine wnad(lu,i,nd,models,data)
                                                                                                                                                                    
                                                                                                                                                                    real*4            models(nd)
                                                                                                                                                                    real*4            data
                                                                                                                                                                    
                                                                                                                                                                    fprintf(lu,rec=i)models,data
                                                                                                                                                                    
                                                                                                                                                                    return
                                                                                                                                                                end
                                                                                                                                                                
                                                                                                                                                                %---------------------------------------------------------------------
                                                                                                                                                                c
                                                                                                                                                                %       read_da - read a direct access file containing
                                                                                                                                                                %                 ne models with dimension nd, ne data
                                                                                                                                                                %                 values and a header of size nh.
                                                                                                                                                                c
                                                                                                                                                                %       Input:
                                                                                                                                                                %             lu                : logical unit of file
                                                                                                                                                                %             fnme              : filename
                                                                                                                                                                %             nh                : length in bytes of file header (minus padding)
                                                                                                                                                                %             mul               : number of records in the header (with padding)
                                                                                                                                                                %             nd                : dimension of parameter space
                                                                                                                                                                %             ne                : number of models in ensemble
                                                                                                                                                                c
                                                                                                                                                                %       Output:
                                                                                                                                                                %             header            : header character string of length nh (char)
                                                                                                                                                                %             data(nd)          : array of data values for each model (real*4)
                                                                                                                                                                %             models(nd,ne)     : array of model values  (real*4)
                                                                                                                                                                c
                                                                                                                                                                %       Comments:
                                                                                                                                                                %                The direct access NAD file format:
                                                                                                                                                                c
                                                                                                                                                                %                VARIABLE       TYPE            SIZE IN BYTES
                                                                                                                                                                %                nd             int             4
                                                                                                                                                                %                ne             int             4
                                                                                                                                                                %                nh             int             4
                                                                                                                                                                %                header         character       nh
                                                                                                                                                                %                models         real*4          4*nd*ne
                                                                                                                                                                %                data           real*4          4*nd
                                                                                                                                                                %                tail           real*4          4
                                                                                                                                                                c
                                                                                                                                                                %                File must contain a single record of length
                                                                                                                                                                %                [4x(3+nd*ne+ne+1) + nh] bytes
                                                                                                                                                                c
                                                                                                                                                                c
                                                                                                                                                                %                                       M. Sambridge, RSES, Nov. 2001
                                                                                                                                                                c
                                                                                                                                                                %---------------------------------------------------------------------
                                                                                                                                                                c
                                                                                                                                                                Subroutine read_da(lu,fnme,nd,ne,nh,mul,header,models,data)
                                                                                                                                                                c
                                                                                                                                                                real*4            models(nd,ne)
                                                                                                                                                                real*4            data(ne)
                                                                                                                                                                character         header(nh)
                                                                                                                                                                %       real*4            tail
                                                                                                                                                                
                                                                                                                                                                character*256     fnme
                                                                                                                                                                c
                                                                                                                                                                %                                              calculate length of header
                                                                                                                                                                iform = 0
                                                                                                                                                                if (mul~=0)iform = 1
                                                                                                                                                                    
                                                                                                                                                                    lenh  = 4*5+nh
                                                                                                                                                                    len   = 4*(nd+1)
                                                                                                                                                                    
                                                                                                                                                                    if (iform==1)
                                                                                                                                                                        
                                                                                                                                                                        open(lu,file=fnme,status='unknown',
                                                                                                                                                                        &          form='unformatted',access='direct',recl=lenh)
                                                                                                                                                                        
                                                                                                                                                                        read(lu,rec=1)idum1,idum2,idum3,idum4,idum5,header
                                                                                                                                                                        
                                                                                                                                                                        close(lu)
                                                                                                                                                                        
                                                                                                                                                                        %          disp(' record length               ',len
                                                                                                                                                                        %          disp(' header length               ',lenh
                                                                                                                                                                        %          disp(' Number of records in header ',mul
                                                                                                                                                                        c
                                                                                                                                                                        open(lu,file=fnme,status='unknown',
                                                                                                                                                                        &          form='unformatted',access='direct',recl=len)
                                                                                                                                                                        
                                                                                                                                                                        for i=1:ne
                                                                                                                                                                            call rnad(lu,i+mul,nd,models(1,i),data(i))
                                                                                                                                                                        end
                                                                                                                                                                        
                                                                                                                                                                        close(lu)
                                                                                                                                                                        
                                                                                                                                                                    else
                                                                                                                                                                        %                                              read in file as
                                                                                                                                                                        %                                              a single record
                                                                                                                                                                        len = 4+nd*ne+ne
                                                                                                                                                                        len = 4*len+nh
                                                                                                                                                                        
                                                                                                                                                                        open(lu,file=fnme,status='unknown',
                                                                                                                                                                        &          form='unformatted',access='direct',recl=len)
                                                                                                                                                                        
                                                                                                                                                                        read(lu,rec=1)i,j,k,kk,header,models,data
                                                                                                                                                                        
                                                                                                                                                                        close(lu)
                                                                                                                                                                        
                                                                                                                                                                    end
                                                                                                                                                                    c
                                                                                                                                                                    return
                                                                                                                                                                end
                                                                                                                                                                
                                                                                                                                                                Subroutine rnad(lu,i,nd,models,data)
                                                                                                                                                                
                                                                                                                                                                real*4            models(nd)
                                                                                                                                                                real*4            data
                                                                                                                                                                
                                                                                                                                                                read(lu,rec=i)models,data
                                                                                                                                                                
                                                                                                                                                                return
                                                                                                                                                            end
                                                                                                                                                            c
                                                                                                                                                            %----------------------------------------------------------------------------
                                                                                                                                                            c
                                                                                                                                                            %       cputime - calls system dependent routine to calculate cputime
                                                                                                                                                            %	  since last call.
                                                                                                                                                            c
                                                                                                                                                            %       Calls dtime.
                                                                                                                                                            %					M. Sambridge, June 2001
                                                                                                                                                            c
                                                                                                                                                            %----------------------------------------------------------------------------
                                                                                                                                                            c
                                                                                                                                                            Function cputime(t1,t2)
                                                                                                                                                            
                                                                                                                                                            real*4 tarray(2)
                                                                                                                                                            
                                                                                                                                                            cputime = dtime(tarray)
                                                                                                                                                            t1 = tarray(1)
                                                                                                                                                            t2 = tarray(2)
                                                                                                                                                            
                                                                                                                                                            return
                                                                                                                                                        end
                                                                                                                                                        c
                                                                                                                                                        %----------------------------------------------------------------------------
                                                                                                                                                        c
                                                                                                                                                        %       NA_header - fprintfs NA-specific information to NAD header.
                                                                                                                                                        c
                                                                                                                                                        %	    This routine adds various NA-header info to
                                                                                                                                                        %	    the header written by the user.
                                                                                                                                                        c
                                                                                                                                                        %       Calls no other routines.
                                                                                                                                                        c
                                                                                                                                                        %					M. Sambridge, June 1999
                                                                                                                                                        c
                                                                                                                                                        %----------------------------------------------------------------------------
                                                                                                                                                        c
                                                                                                                                                        Subroutine NA_header
                                                                                                                                                        &             (lu,fnme,header,nh_max,nh,nd,
                                                                                                                                                        &              range,scales,n1,n2,n3,nhu)
                                                                                                                                                        
                                                                                                                                                        real*4		  range(2,nd)
                                                                                                                                                        real*4		  scales(nd+1)
                                                                                                                                                        character*(*)     header
                                                                                                                                                        character*256 	  fnme
                                                                                                                                                        
                                                                                                                                                        %                                               calculate total header length
                                                                                                                                                        rlen = 3*nd + 4
                                                                                                                                                        len = 4*rlen+nh
                                                                                                                                                        nhu = nh
                                                                                                                                                        nh_na   = 4*rlen
                                                                                                                                                        nh_tot  = len
                                                                                                                                                        %disp(' nh_user = ',nhu
                                                                                                                                                        %disp(' nh_na   = ',nh_na
                                                                                                                                                        %disp(' nh_tot  = ',nh_tot
                                                                                                                                                        
                                                                                                                                                        if (nh_tot>nh_max)
                                                                                                                                                            disp('')
                                                                                                                                                            disp(' Error - header array too small'
                                                                                                                                                            disp('')
                                                                                                                                                            disp('         current size = ',nh_max
                                                                                                                                                            disp('        required size = ',nh_tot
                                                                                                                                                            disp('')
                                                                                                                                                            disp(' Remedy - adjust nh_max in parameter',
                                                                                                                                                            &                        ' file and recompile'
                                                                                                                                                            disp('')
                                                                                                                                                            call na_abort
                                                                                                                                                        end
                                                                                                                                                        
                                                                                                                                                        %					fprintf out header information
                                                                                                                                                        call fprintf_header
                                                                                                                                                        &       (lu,fnme,len,nd,nh,range,scales,n1,n2,n3,header)
                                                                                                                                                        
                                                                                                                                                        nh = nh_tot
                                                                                                                                                        %					read header information
                                                                                                                                                        %					into character string
                                                                                                                                                        call read_header
                                                                                                                                                        &       (lu,fnme,nh,len,header)
                                                                                                                                                        
                                                                                                                                                        
                                                                                                                                                        %fprintf(50,*)header(nh_na+1:nh_tot)
                                                                                                                                                        
                                                                                                                                                        return
                                                                                                                                                    end
                                                                                                                                                    c
                                                                                                                                                    %----------------------------------------------------------------------------
                                                                                                                                                    c
                                                                                                                                                    %       read_header - converts header information into a character
                                                                                                                                                    %	      string by writing it to a direct access file
                                                                                                                                                    %	      and  reading it back as a character string
                                                                                                                                                    c
                                                                                                                                                    %       Calls no other routines.
                                                                                                                                                    c
                                                                                                                                                    %					M. Sambridge, June 1999
                                                                                                                                                    c
                                                                                                                                                    %----------------------------------------------------------------------------
                                                                                                                                                    c
                                                                                                                                                    %                                               open direct access
                                                                                                                                                    %					temporary file
                                                                                                                                                    Subroutine read_header
                                                                                                                                                    &  	   (lu,fnme,nh,len,header)
                                                                                                                                                    
                                                                                                                                                    character         header(nh)
                                                                                                                                                    character*256 	  fnme
                                                                                                                                                    
                                                                                                                                                    open(lu,file=fnme,status='unknown',
                                                                                                                                                    &       form='unformatted',access='direct',recl=len)
                                                                                                                                                    
                                                                                                                                                    read(lu,rec=1)header
                                                                                                                                                    
                                                                                                                                                    close(lu)
                                                                                                                                                    
                                                                                                                                                    return
                                                                                                                                                end
                                                                                                                                                c
                                                                                                                                                %----------------------------------------------------------------------------
                                                                                                                                                c
                                                                                                                                                %       fprintf_header - converts header information into a character
                                                                                                                                                %	       string by writing it to a direct access file
                                                                                                                                                %	       and  reading it back as a character string
                                                                                                                                                c
                                                                                                                                                %       Calls no other routines.
                                                                                                                                                c
                                                                                                                                                %					M. Sambridge, June 1999
                                                                                                                                                c
                                                                                                                                                %----------------------------------------------------------------------------
                                                                                                                                                c
                                                                                                                                                %                                               open direct access
                                                                                                                                                %					temporary file
                                                                                                                                                Subroutine fprintf_header
                                                                                                                                                &  	   (lu,fnme,len,nd,nh,
                                                                                                                                                &              range,scales,n1,n2,n3,header)
                                                                                                                                                
                                                                                                                                                real*4		  range(2,nd)
                                                                                                                                                real*4		  scales(nd+1)
                                                                                                                                                character         header(nh)
                                                                                                                                                character*256 	  fnme
                                                                                                                                                
                                                                                                                                                open(lu,file=fnme,status='unknown',
                                                                                                                                                &       form='unformatted',access='direct',recl=len)
                                                                                                                                                
                                                                                                                                                fprintf(lu,rec=1)n1,n2,n3,range,scales,header
                                                                                                                                                
                                                                                                                                                close(lu)
                                                                                                                                                
                                                                                                                                                return
                                                                                                                                            end
                                                                                                                                            c
                                                                                                                                            %----------------------------------------------------------------------------
                                                                                                                                            c
                                                                                                                                            %       NA_display - fprintfs out current best fit model
                                                                                                                                            %	     to LU lu_det with misfit information
                                                                                                                                            c
                                                                                                                                            %       Calls no other routines.
                                                                                                                                            c
                                                                                                                                            %					M. Sambridge, Sept. 1999
                                                                                                                                            c
                                                                                                                                            %----------------------------------------------------------------------------
                                                                                                                                            c
                                                                                                                                            Subroutine NA_display
                                                                                                                                            &             (lu_det,model_opt,it,nd,ntot,
                                                                                                                                            &              mfitmin,mfitminc,mfitmean,mopt)
                                                                                                                                            
                                                                                                                                            real*4		model_opt(*)
                                                                                                                                            real*4		mfitmin,mfitmean,mfitminc
                                                                                                                                            
                                                                                                                                            fprintf(lu_det,*)' Iteration              ',
                                                                                                                                            &                 '           : ',it
                                                                                                                                            fprintf(lu_det,*)' Total number of samples ',
                                                                                                                                            &                 '          : ',ntot
                                                                                                                                            fprintf(lu_det,*)
                                                                                                                                            &  ' Minimum misfit                    : ',mfitmin
                                                                                                                                            fprintf(lu_det,*)
                                                                                                                                            &  ' Minimum misfit for this iteration : ',mfitminc
                                                                                                                                            fprintf(lu_det,*)
                                                                                                                                            &  ' Mean misfit in this iteration     : ',mfitmean
                                                                                                                                            fprintf(lu_det,*)
                                                                                                                                            &  ' Index of best fitting model       : ',mopt
                                                                                                                                            fprintf(lu_det,*)
                                                                                                                                            
                                                                                                                                            %					fprintf out optimum model
                                                                                                                                            
                                                                                                                                            fprintf(lu_det,*)' Parameter   best fitting model'
                                                                                                                                            for i=1:nd
                                                                                                                                                fprintf(lu_det,*)' ',i,' : ',model_opt(i)
                                                                                                                                            end
                                                                                                                                            fprintf(lu_det,100)
                                                                                                                                            
                                                                                                                                            100    format(/,72("-")/)
                                                                                                                                            
                                                                                                                                            
                                                                                                                                            return
                                                                                                                                        end
                                                                                                                                        c
                                                                                                                                        %----------------------------------------------------------------------------
                                                                                                                                        c
                                                                                                                                        %       jumble - randomly re-arranges input array
                                                                                                                                        c
                                                                                                                                        %       Calls ran3 and assumes that it has been initialized.
                                                                                                                                        c
                                                                                                                                        %                                               M. Sambridge, Oct. 1999.
                                                                                                                                        c
                                                                                                                                        %----------------------------------------------------------------------------
                                                                                                                                        c
                                                                                                                                        Subroutine jumble(iarr,arr,n)
                                                                                                                                        
                                                                                                                                        integer		iarr(n)
                                                                                                                                        real*4		arr(n)
                                                                                                                                        
                                                                                                                                        %disp(' input to jumble'
                                                                                                                                        %disp((iarr(k),k=1,n)
                                                                                                                                        
                                                                                                                                        rn = n
                                                                                                                                        for j=1:n
                                                                                                                                            val = ran3(iseed)
                                                                                                                                            k = 1 + int(val*rn)
                                                                                                                                            if (k==n+1)
                                                                                                                                                k = n
                                                                                                                                            else(k>n)
                                                                                                                                                disp(' error in jumble k'
                                                                                                                                                &        ,k,' val',val,' rn',rn
                                                                                                                                                call na_abort
                                                                                                                                            end
                                                                                                                                            ival = iarr(j)
                                                                                                                                            iarr(j) = iarr(k)
                                                                                                                                            iarr(k) = ival
                                                                                                                                            val = arr(j)
                                                                                                                                            arr(j) = arr(k)
                                                                                                                                            arr(k) = val
                                                                                                                                        end
                                                                                                                                        
                                                                                                                                        %disp(' output of jumble'
                                                                                                                                        %disp((iarr(k),k=1,n)
                                                                                                                                        
                                                                                                                                        return
                                                                                                                                    end
                                                                                                                                    c
                                                                                                                                    %----------------------------------------------------------------------------
                                                                                                                                    c
                                                                                                                                    %       na_abort - randomly re-arranges input array
                                                                                                                                    %					Used in place of a stop
                                                                                                                                    %					to ensure that all
                                                                                                                                    %					processes are stopped
                                                                                                                                    %					in MPI mode
                                                                                                                                    c
                                                                                                                                    %----------------------------------------------------------------------------
                                                                                                                                    c
                                                                                                                                    Subroutine na_abort
                                                                                                                                    
                                                                                                                                    cifdef MPI
                                                                                                                                    include "mpif.h"
                                                                                                                                    integer	ierr
                                                                                                                                    call mpi_abort(MPI_COMM_WORLD, ierr)
                                                                                                                                    cendif
                                                                                                                                    stop
                                                                                                                                end
