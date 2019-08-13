c###############################################################################
	subroutine delphi(acent_xdim, energy_xdim, site_xdim, 
     & py_igrid, py_scale, py_repsin, py_repsout, py_radprb, py_conc, 
     & py_ibctyp, py_res2, py_nlit, py_acent, py_energy, py_site, 
     & py_nonit, py_relfac, py_relpar, py_pbx, py_pby,
     & py_in_frc, natom, nmedia, nobject,	
     & py_i_atpos, py_i_rad3, py_i_chrgv4, py_i_atinf, py_i_medeps,
     & py_i_iatmmed, py_i_dataobject, py_i_phimap4, py_scale1,
     & py_out_phi, py_i_sitpot, esolvation, py_isurftype, py_parallel)	

c       PARSING ARGUMENTS FROM PYTHON
c	python stores floats as doubles
c       py_variables are not used directly but rather reassigned
	integer*8 acent_xdim, energy_xdim, site_xdim
	integer*8 py_igrid
	real*8  py_scale
	real*8  py_repsin
	real*8  py_repsout
	real*8  py_radprb
	real*8  py_conc
	integer*8 py_ibctyp
	real*8  py_res2
	integer*8 py_nlit
	real*8  py_acent(acent_xdim)
	character*(*) py_energy(energy_xdim), py_site(site_xdim)

c	Membrane arguments
	integer*8 py_nonit
	real*8    py_relfac, py_relpar
	logical*8 py_pbx, py_pby
		
        character*(*) py_in_frc
	integer*8 natom
	integer*8 nmedia
	integer*8 nobject

c	pointers of arrays
	integer*8 py_i_atpos
	integer*8 py_i_rad3
	integer*8 py_i_chrgv4
	integer*8 py_i_atinf
	integer*8 py_i_medeps
	integer*8 py_i_iatmmed
	integer*8 py_i_dataobject

c	output
	integer*8 py_i_phimap4
	real*8    py_scale1
	logical*8 py_out_phi
	integer*8 py_i_sitpot
	real*8    esolvation

c       nanoshaper
	integer*8 py_isurftype

c	solver
	logical*8 py_parallel
	
c       END PARSING ARGUMENTS FROM PYTHON
	
Cf2py   threadsafe
Cf2py   intent(in) py_igrid
Cf2py   intent(in) py_scale
Cf2py   intent(in) py_repsin
Cf2py   intent(in) py_repsout
Cf2py   intent(in) py_radprb
Cf2py   intent(in) py_conc
Cf2py   intent(in) py_ibctyp
Cf2py   intent(in) py_res2
Cf2py   intent(in) py_nlit
Cf2py   intent(in) py_acent
Cf2py   intent(in) py_energy
Cf2py   intent(in) py_site
Cf2py   intent(in) py_nonit
Cf2py   intent(in) py_relfac
Cf2py   intent(in) py_relpar
Cf2py   intent(in) py_pbx
Cf2py   intent(in) py_pby
Cf2py   intent(in) py_in_frc	
Cf2py   intent(in) natom
Cf2py   intent(in) nmedia
Cf2py   intent(in) nobject
Cf2py   intent(in) py_i_atpos
Cf2py   intent(in) py_i_rad3
Cf2py   intent(in) py_i_chrgv4
Cf2py   intent(in) py_i_iatmmed
Cf2py   intent(in) py_i_dataobject
Cf2py   intent(in) py_i_phimap4
Cf2py   intent(in) py_scale1
Cf2py   intent(in) py_out_phi	
Cf2py   intent(in) py_i_sitpot
Cf2py   intent(in,out) esolvation
Cf2py   intent(in) py_isurftype
Cf2py   intent(in) py_parallel
c       b++++++++++++++++++++++++++++++++++++++++++++
cpgi$ IF DEFINED (PC)
#ifdef PC
cDEC$ IF DEFINED (PC)
        use dflib
        use dfport
cDEC$ END IF
#endif
c       e++++++++++++++++++++++++++++++++++++++++++++

	pointer (py_i_sitpot,sitpot)

	include 'qdiffpar4.h'
	include 'qlog.h'

c------------------------------------------------------------------
	real sitpot(natom)
	
c       MODIFIED EXISTING ARRAYS
	real         atpos(natom*3)
	real         rad3(natom)
	real         chrgv4(natom)
	character*15 atinf(natom)
	real         medeps(0:nmedia+1)
	integer      iatmmed(natom+nobject)
	character*96 dataobject(nobject,2)

c 	ORIGINAL VARIABLES
     	character*50 fnma(30)
	character*80 argnam
	real cmin(3),cmax(3),cran(3),xo(3)
	dimension xn(3),xn2(1),ibgrd(3,1),dbval(0:1,0:6,0:1)
	dimension chrgv2(1),idpos(1),db(6,1)
	dimension sf1(1),sf2(1),qval(ngcrg)
	dimension sfd(5,0:1),scspos(3,1)
	dimension cqplus(3),cqmin(3)
	integer ifnma(30),iqpos(ngcrg),arglen
        character*24 day
c b+++++++++++++++++++
        integer crgatn(1),nqgrdtonqass(1),extracrg
        integer ndistr,idimtmp
	logical uniformdiel
        integer idirectalg,numbmol,icount1a,icount1b
        real atmeps(1),qfact
	dimension ibndx(1),ibndy(1),ibndz(1),bndx1(1),bndx2(1),bndx3(1)
        dimension bndx4(1),qmap1(1),qmap2(1),debmap1(1),debmap(2)
c       character*96 dataobject(nobjectmax,2)
c       character*80 datadistr(ndistrmax)
c       character*80 strtmp
c
c pointer based arrays
        real limobject(1),cgbp(5,1),chgpos(3,1)
c e+++++++++++++++++++

	integer iepsmp(1)
        logical*1 idebmap(1)
        real phimap(1),atmforce(3,1)
        real phimap1(1),phimap2(1),phimap3(1)
	integer neps(1),keps(1)
	
!!!!!!!!!!!!!!!!! Nano Shaper variables !!!!!!!!!!!!!!!!!
! NanoShaper flag
	integer useNanoShaper  
! Name of the loaded module 
!(shared object name without file extension)
	character*80 surfModuleFile      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

c       jbc !!!!!!!!!!!!!!!!!! New GPU-enabled solver !!!!!!!!!!!!!!!  jbc
	logical*1 imanual2
	character*80 solverModuleFile

c	profiling purposes
	character*25 filename_test
	start = cputime(start)

c	finish = finish - start
	write(6,*)'  '
	write(6,*)'loading variables time was (sec) ',finish
	write(6,*)'  '

	useNanoShaper = 0

c jbc !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
c------------------------------------------------------------------
	icount1a=0
	icount1b=0
	start = cputime(0.0)
c	call wrt(1)


c b+++++++++++++++++++++++++++++++++++++++++
c       tolerance in small numbers
	tol=1.e-7
c       debug flag
	debug=.false.
	if (debug) write(6,*)"WARNING: working in DEBUGging mode"
c e+++++++++++++++++++++++++++++++++++++++++
c
c	call rdlog(fnma,ifnma); find logical links..
c

c	i=iargc()
c	argnam=" "
c	if(i.gt.0)then
c	   call getarg(1,argnam)
c	   call namlen(argnam,arglen)
c	end if
c
c	i=1
c	argnam=prm
c	call namlen(argnam,arglen)
c	call qqint(i,argnam,arglen)
        
	i_atpos=py_i_atpos
	i_rad3=py_i_rad3
	i_chrgv4=py_i_chrgv4
	i_atinf=py_i_atinf
	i_medeps=py_i_medeps
	i_iatmmed=py_i_iatmmed
	i_dataobject=py_i_dataobject
	i_phimap4=py_i_phimap4	
	call defprm()

c
c       print parsed arguments
c

	PRINT *, 'igrid     -> ' , igrid     , py_igrid   
	PRINT *, 'scale     -> ' , scale     , py_scale   
	PRINT *, 'repsin    -> ' , repsin    , py_repsin  
	PRINT *, 'repsout   -> ' , repsout   , py_repsout 
	PRINT *, 'radprb(1) -> ' , radprb(1) , py_radprb  
	PRINT *, 'conc(1)   -> ' , conc(1)   , py_conc    
	PRINT *, 'ibctyp    -> ' , ibctyp    , py_ibctyp  
	PRINT *, 'res2      -> ' , res2      , py_res2    
	PRINT *, 'nlit      -> ' , nlit      , py_nlit
	PRINT *, 'nnit      -> ' , nnit      , py_nonit
	PRINT *, 'py_relfac -> ' , uspec     , py_relfac
        PRINT *, 'py_relpar -> ' , relpar    , py_relpar
        PRINT *, 'py_pbx    -> ' , iper(1)   , py_pbx
        PRINT *, 'py_pby    -> ' , iper(2)   , py_pby
	PRINT *, 'py_acent  > (', acent(1), acent(2), acent(3), ') (', py_acent(1), py_acent(2), py_acent(3), ')'
	PRINT *, 'energy      -> ', py_energy
	PRINT *, 'sites -> ', py_site
	PRINT *, 'natom -> ', natom
        PRINT *, 'nmedia -> ', nmedia
        PRINT *, 'nobject -> ', nobject



        


	

c	statement parsing
c

	igrid     = py_igrid   
	scale     = py_scale   
	repsin    = py_repsin  
	repsout   = py_repsout 
	radprb(1) = py_radprb
	radprb(2) = py_radprb
c	salt
	rionst    = py_conc    
	conc(1)   = py_conc
c	dfact     = 0.01990076478*sqrt(temperature*repsout)	
c	deblen    = dfact/sqrt(rionst)
c	boundary condition type
	ibctyp    = py_ibctyp
c       maxc
	res2      = py_res2
c	linit
	nlit      = py_nlit    
	iautocon  = .false.
c	
c	function parsing
c
c       grid offset
	acent(1) = py_acent(1)
	acent(2) = py_acent(2)
	acent(3) = py_acent(3)	
	iacent=.true.

c	coulombic energy	
	if ( ANY( py_energy=="c" ) ) logc=.true.
c	solvation energy	
	if ( ANY( py_energy=="s" ) ) logs=.true.

c	atom
	if ( ANY( py_site=="a" ) ) isita=.true.
c	charge
	if ( ANY( py_site=="q" ) ) isitq=.true.
c	potential
	if ( ANY( py_site=="p" ) ) isitp=.true.	   

	if ((isita .eqv. .true.) .or.
     &  (isitq .eqv. .true.) .or.
     &  (isitp .eqv. .true.)) then
	   isite=.true.
	endif
c
c
c
	
c
c       input file names parsing
c
	ias=1


	if (py_in_frc == 'self') iself=.true.


	nnit = py_nonit	
	if ((nnit .ne. 0) .and. (nnit.le.20)) then
	   write(6,*)'At least 30 nonlinear iterations'
	   nnit=30
	end if

	
	uspec = py_relfac
	if (py_relfac .ne. 0.0) then
	   iuspec = .true.
	endif
	
	relpar = py_relpar
	if (py_relpar .ne. 0.0) then
	   imanual=.true.
	endif
	
	if (py_pbx .eqv. .true.) then
	   iper(1)=.true.
	endif
	
	if (py_pby .eqv. .true.) then
	   iper(2)=.true.
	endif


	isurftype=py_isurftype

	
	if (py_parallel .eqv. .true.) then	                       
	   solverModuleFile = 'libDelphiGpuSolver'//CHAR(0)
	   isolvarch=1
	else
	   solverModuleFile = 'surfacemodule'//CHAR(0)
	   isolvarch=0
	endif

	
	call rdprm

c       b++++++++++++++++++++++++++++++++++++++
c       create suitable pdb file, if necessary
	if (icreapdb) call creapdb(repsout,numbmol)
c       e++++++++++++++++++++++++++++++++++++++
c       
c       set hash tables from size and charge files
c
c	call rdhrad(*999)
c	if(isolv) call rdhcrg(*999)

c       read in pdbfile and, if necesssary, assign charges/radii
c
c       b++++++++++++++++++++++++++++++
        uniformdiel=.true.
c       i_iatmmed=memalloc(i_iatmmed,4,natmax)
c       e+++++++++++++++++++++++++++++

	ndistr=0
	i_datadistr=memalloc(i_datadistr,4,80*ndistr)
	
c	call setrc(natom,*999,nmedia,nobject,ndistr)

        idirectalg=1
        if (nmedia.gt.1) then
	   write(6,*)'Attention, many dielectrics! not all the surface
     & charge is facing the solution!!'
	   idirectalg=1
        end if
        write(6,*)'Direct mapping of epsilon: (0/1)(n/y)',idirectalg
c       b++++++++++++++++++++++++++++++
c       i_iatmmed=memalloc(i_iatmmed,4,natom+nobject)	
c       epsin = repsin /epkt	!!!!!
	epsin=medeps(1)
c       e+++++++++++++++++++++++++++++
c

	finish=cputime(start)
	write(6,*) 'time to read in and/or assign rad/chrg=',finish
c       b+++++++++++++++++++++
        i_limobject=memalloc(i_limobject,realsiz,nobject*3*2) 
c       find box enclosing each object at this moment expressed in Angstrom
        call extrmobjects(nobject,scale,natom,numbmol,verbose)
c       e+++++++++++++++++++++
c       find extrema
c
	call extrm(natom,igrid,cmin,cmax,nobject)

	PRINT *, 'cmin', cmin
	PRINT *, 'cmax', cmax

c       b++++++++++++++++++++
        pmid(1)=(cmax(1)+cmin(1))/2.
        pmid(2)=(cmax(2)+cmin(2))/2.
        pmid(3)=(cmax(3)+cmin(3))/2.
c
c       calculate offsets
c
        call off(oldmid,pmid,*999)

        cran(1)=cmax(1)-cmin(1)
        cran(2)=cmax(2)-cmin(2)
        cran(3)=cmax(3)-cmin(3)
c       rmaxdim=max(cran(1),cran(2),cran(3))
c       ++++modificato per convenienza di dimensionamento
        rmaxdim=2.*max(abs(cmax(1)-oldmid(1)),abs(cmin(1)-oldmid(1)),
     &          abs(cmax(2)-oldmid(2)),abs(cmin(2)-oldmid(2)),
     &          abs(cmax(3)-oldmid(3)),abs(cmin(3)-oldmid(3)))
c       e+++++++++++++++++++
c
c       atpos=atom positions, oldmid=current midpoints, rmaxdim=largest dimension
c       rad3= radii in angstroms. calculate scale according to them and
c       to the percent box fill
c
	if(igrid.eq.0)then

	   if(scale.eq.10000.)scale=2.0
	   if(perfil.eq.10000.)perfil=80.
	   igrid=scale*100./perfil*rmaxdim
	   
	elseif(scale.eq.10000.)then
	   
	   if(perfil.eq.10000.)then
	      scale=2.0
	      perfil=100.*rmaxdim*scale/float(igrid-1)
	   else
	      scale=float(igrid-1)*perfil/(100.*rmaxdim)
	   endif
	else
	   perfil=100.*rmaxdim*scale/float(igrid-1)
	endif
c       calcolo quante deblen di soluzione sono contenute nel box
	if (deblen.lt.1000000.) then
	   debnum=(100./perfil-1.)*rmaxdim/deblen
	   write(6,*)'Debye Lengths contained in the finite diff. box',debnum
	end if
c       if(rionst.gt.0.0.and.exrad.lt.1.e-6)exrad=2.0
c       
	irm=mod(igrid,2)
	if(irm.eq.0)igrid=igrid+1
c       
	if(igrid.gt.ngrid)then
	   igrid=ngrid
	   write(6,*)'igrid = ',igrid,' exceeds ','ngrid = ',ngrid,'
     & so reset'
	   scale=float(igrid-1)*perfil/(100.*rmaxdim)
	endif
c       b+++++++++++++ added dependence by nmedia                              
c       medeps(0)=epsout
	call wrtprm(nmedia,cmin,cmax,cran)
	
c       e+++++++++++++
	ngp=igrid*igrid*igrid+1
	nhgp=ngp/2
	nbgp=igrid*igrid+1
	i_iepsmp=memalloc(i_iepsmp,4,3*igrid*igrid*igrid)
	i_idebmap=memalloc(i_idebmap,1,igrid*igrid*igrid)
c
c       calculate offsets
c
c       call off(oldmid,*999) vecchio posto
c
        xl1=oldmid(1)-(1.0/scale)*float(igrid+1)*0.5
        xl2=oldmid(2)-(1.0/scale)*float(igrid+1)*0.5
        xl3=oldmid(3)-(1.0/scale)*float(igrid+1)*0.5
        xr1=oldmid(1)+(1.0/scale)*float(igrid+1)*0.5
        xr2=oldmid(2)+(1.0/scale)*float(igrid+1)*0.5
        xr3=oldmid(3)+(1.0/scale)*float(igrid+1)*0.5

	if(cmin(1).lt.xl1.or.cmin(2).lt.xl2.or.cmin(3).lt.xl3.or.
     &	cmax(1).gt.xr1.or.cmax(2).gt.xr2.or.cmax(3).gt.xr3)
     &	write(6,*)
     &	'!!! WARNING: part of system outside the box!'
c       convert atom coordinates from angstroms to grid units
c
	i_xn2=memalloc(i_xn2,realsiz,3*natom)
	call grdatm(natom,igrid,scale,oldmid)
	
!!!!!!!!!!!!!!!!!! NanoShaper !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	if (isurftype.ne.-1) then
	   useNanoShaper = 1
	   write(*,*), "Surface id -> ",isurftype
!       write custom settings file for NanoShaper
c	   open (unit = 010, file = "custom.prm")
c	   if(isurftype.eq.0)  write(010,*), "Surface = ses"
c	   if(isurftype.eq.1)  write(010,*), "Surface = skin"
c	   if(isurftype.eq.2)  write(010,*), "Surface = blobby"
c	   if(isurftype.eq.3)  write(010,*), "Surface = mesh"
c	   if(isurftype.eq.4)  write(010,*), "Surface = msms"
c	   close(010)
	endif

	
	if (useNanoShaper.eq.1) then
	   ibmx = 2000000
! preallocate a big enough vector of bgps
	   i_ibgrd = memalloc(i_ibgrd,realsiz,ibmx*3)
! preallocate a big enough vector of bgps
	   i_atsurf = memalloc(i_atsurf,realsiz,ibmx)
!allocate(ns_ibgrd(ibmx*3))
	   i_scspos = memalloc(i_scspos,realsiz,ibmx*3)
!allocate(ns_scspos(ibmx*3))
	   i_scsnor = memalloc(i_scsnor,realsiz,ibmx*3)
!allocate(ns_scsnor(ibmx*3))
	endif
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	
c       verify if dielectric is uniform 
	
        do i = 0,nmedia-1
	   if (medeps(i).ne.medeps(i+1)) uniformdiel=.false.
        end do

	
c       b++++++++write visual.txt to dialog with the GUI
c       call wrtvisual(natom,igrid,nobject,scale,oldmid,nmedia,epkt)
	
c       
c       make some charge arrays for boundary conditions etc.
c       
        if(isolv)then
c       b++++++++++increased the starting dimension of crgatn and..+++++
	   extracrg=0
	   if (ndistr.gt.0)  extracrg=igrid**3
	   i_crgatn=memalloc(i_crgatn,4,natom+extracrg)
	   i_chrgv2=memalloc(i_chrgv2,realsiz,4*(natom+extracrg))
	   i_nqgrdtonqass=memalloc(i_nqgrdtonqass,4,natom+extracrg)
	   i_atmeps=memalloc(i_atmeps,realsiz,natom+extracrg)
	   i_chgpos=memalloc(i_chgpos,realsiz,3*ncrgmx)
c       e+++++++++++++++++++++++++++++++++++++++++++++++++++++++
c       
	   call crgarr(ncrgmx,cqplus,cqmin,atpos,igrid,natom,nqass
     &  ,nqgrd,qmin,qnet,qplus,nmedia,ndistr,scale,oldmid,nobject,
     &  radpolext,extracrg,realsiz,verbose)
c
c       b+++++++++++++++++++++++++++++
	   i_crgatn=memalloc(i_crgatn,4,nqass)
	   i_nqgrdtonqass=memalloc(i_nqgrdtonqass,4,nqgrd)
	   i_atmeps=memalloc(i_atmeps,realsiz,nqass)
	   i_chrgv2=memalloc(i_chrgv2,realsiz,4*nqgrd)
c       wr and sd bug fix: too early deallocation
c       if(.not.isite) i_iatmmed=memalloc(i_iatmmed,0,0)
	   i_chgpos=memalloc(i_chgpos,realsiz,3*nqass)
	   
	   if(logs.or.lognl)then
c       e++++++++++++++++++++++++++++
	      ico=0
	      do ic=1,nqass
		 cx1=chgpos(1,ic)
		 cx2=chgpos(2,ic)
		 cx3=chgpos(3,ic)
		 if(cx1.lt.xl1.or.cx1.gt.xr1.or.cx2.lt.xl2.or.
     &     cx2.gt.xr2.or.cx3.lt.xl3.or.cx3.gt.xr3)then
c
		    if (crgatn(ic).lt.0) then
c       b+++++++++++++++++
		       write(6,
     &  '(''!WARNING: distribution'',I4,'' outside the box'')')
     &  (-crgatn(ic))
		    else
		       if (crgatn(ic).gt.natom) then
                  write(6,
     &'(''WARNING:crg'',I4,'',object'',I4,'' outside the box'',3f8.3)'
     &  )ic,(crgatn(ic)-natom),cx1,cx2,cx3  
	       else
c       e+++++++++++++++++
                  write(6,
     &  '(''!!! WARNING : charge '',a15,'' outside the box'')')
     &  atinf(crgatn(ic))
	       endif
	    endif
c
	    ico=1
	 endif
	end do
	if(ico.gt.0.and.ibctyp.ne.3)then
	   write(6,*)'CHARGES OUTSIDE THE BOX AND NOT DOING FOCUSSING,
     &THEREFORE STOP'
	   stop
	end if
        endif
        endif
!!!!!!!!!!!!!!!!!! NanoShaper !!!!!!!!!!!!!!!!!!!!!!
	if (useNanoShaper.eq.1) then
	   surfModuleFile = 'libDelphiSurface'//CHAR(0)

	   call surfacemodule(surfModuleFile,
     &  oldmid(1)-(1.0/scale)*float(igrid-1)*0.5,
     &  oldmid(2)-(1.0/scale)*float(igrid-1)*0.5,
     &  oldmid(3)-(1.0/scale)*float(igrid-1)*0.5,
     &  oldmid(1)+(1.0/scale)*float(igrid-1)*0.5,
     &  oldmid(2)+(1.0/scale)*float(igrid-1)*0.5,
     &  oldmid(3)+(1.0/scale)*float(igrid-1)*0.5,
     &  oldmid(1),oldmid(2),oldmid(3),rmaxdim,
     &  perfil,i_iepsmp,igrid,scale,i_scspos,i_scsnor,
     &  i_idebmap,i_ibgrd,5+natom,ibnum,
     &  natom,i_atpos,i_rad3,0,ibmx,radprb(1),exrad,i_atinf,i_atsurf)
      
	else        
  
	   call epsmak(ibnum,natom,oldmid,uniformdiel,
     &  nobject,nmedia,numbmol)

	endif
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

c       sd e++++++++++++
c       write details

	
        call wrtadt(perfil,cmin,cmax,cran,oldmid,scale,natom
     &  ,nqass,qnet,qplus,cqplus,qmin,cqmin,isolv)

c       ++++++++++++++++++++++++++++++++++++
	if(isolv) then
	   write(6,*) 'number of dielectric boundary points',ibnum
	   if(iexun.and.(ibnum.eq.0)) then
	      write(6,*) "exiting as no boundary elements and"
	      write(6,*) "uniform dielectric exit flag has been set"
	      goto 998
	   end if
c###################
	   call dbsfd(dbval,sfd)
c
c       nsp=ibnum+1000 this on average saves memory but might encounter problems
	   nsp=2*(ibnum+1)
	   i_db= memalloc(i_db,realsiz,6*nsp)
	   i_idpos= memalloc(i_idpos,4,nsp)
	   i_sf1= memalloc(i_sf1,realsiz,nhgp)
	   i_sf2= memalloc(i_sf2,realsiz,nhgp)
	   
	   call mkdbsf(ibnum,nsp,dbval,icount2a,icount2b,sfd,natom,
     &  nmedia,idirectalg,nobject)

c
c write dielectric map
c
	   if(epswrt)then
	      imaxwrd = igrid/16 + 1
	      i_neps= memalloc (i_neps,4,3*imaxwrd*igrid*igrid)
	      i_keps= memalloc (i_keps,4,imaxwrd*igrid*igrid)
c       b++++++++++++++++++++++++++++++++Aug 2011++++++++++++++++
	      idimtmp=0
	      if (epsfrm.eq.1) then
		 i_tmpmap=memalloc(i_tmpmap,4,3*igrid*igrid*igrid)
		 idimtmp=igrid
	      end if
	      call wrteps(imaxwrd,natom+nobject+2,nmedia,idimtmp)
	      i_neps= memalloc (i_neps,0,0)
	      i_keps= memalloc (i_keps,0,0)
	      if (epsfrm.eq.1) 
     &       i_tmpmap=memalloc(i_tmpmap,0,0)
c e++++++++++++++++++++++++++++++++++++++++++++++++++++++++
	   endif
c	call chkeps(natom,atpos,rad3,ibnum,ibgrd,igrid,scale,radprb)



c
c       make qval and other linear charge arrays for the solver
c
	   i_phimap=memalloc(i_phimap,realsiz,igrid*igrid*igrid)
	   if(isph) then
	      call setfcrg(nqgrd,nqass,icount1a,icount1b,nmedia,natom,
     &idirectalg,nobject)
	   else
	      call setcrg(nqgrd,nqass,icount1a,icount1b,nmedia,natom,
     &idirectalg,nobject)
	   endif
c       b++++++++++++++++++++
	   i_nqgrdtonqass=memalloc(i_nqgrdtonqass,0,0)
c       e++++++++++++++++++++
	   finish=cputime(start)
	   write(6,*) 'iepsmp to db, and charging done at', finish
	   write(6,*) 'number of grid points assigned charge', icount1b

c       I can''t get rid of iepsmp if I have to calculate the nonlin energy
	   i_iepsmp= memalloc(i_iepsmp,0,0)
c
c       calculate boundary conditions
c
	   call setbc(qplus,qmin,cqplus,cqmin,nqass,natom,ibnum,py_scale1)
c
	   i_phimap1= memalloc(i_phimap1,realsiz,nhgp)
	   i_phimap2= memalloc(i_phimap2,realsiz,nhgp)
	   i_phimap3= memalloc(i_phimap3,realsiz,ngp)

	   i_ibndx= memalloc(i_ibndx,4,nbgp)
c       trasformati in interi
	   i_ibndy= memalloc(i_ibndy,4,nbgp)
	   i_ibndz= memalloc(i_ibndz,4,nbgp)
	   if(iuspec) then
	      spec=uspec
	      write(6,*) "using entered value for relaxation of: ",spec
	   else
	      call relfac(idpos,db,sf1,sf2,icount2a,icount2b,spec,nsp,
     &	phimap1,phimap2,phimap3,ibndx,ibndy,ibndz,idirectalg)
	   end if
c
	   noit=int(7.8/log(1.0 + sqrt(1-spec)))
	   write(6,*) 'estimated iterations to convergence',noit
	   if(iautocon) nlit = noit
c       
	   i_bndx1= memalloc(i_bndx1,realsiz,nbgp)
	   i_bndx2= memalloc(i_bndx2,realsiz,nbgp)
	   i_bndx3= memalloc(i_bndx3,realsiz,nbgp)
	   i_bndx4= memalloc(i_bndx4,realsiz,nbgp)
c   

	   finish=cputime(start)
	   write(6,*)'  '
	   write(6,*)'setup time was (sec) ',finish
	   write(6,*)'  '

c
c       iterate
c
c       b++++++++++++++++++++++++++++++++++++++++++++++++++++++++
	   if (float(ibc)/ibnum.gt..3.and.iautocon) then
	      nlit=nlit*ibc/(.3*ibnum)
	      write(6,*)'Re-estimated iterations now :',nlit
	   end if
c       e++++++++++++++++++++++++++++++++++++++++++++++++++++++++
	   print *, 'isolvarch', isolvarch
	   if(nnit.eq.0.or.rionst.lt.1.e-6) then
	      if (isolvarch==0) then

        call itit(idpos,db,sf1,sf2,iqpos,qval,icount2a,icount2b,
     &	icount1a,icount1b,spec,nsp,phimap,phimap1,phimap2,phimap3,
     &	ibndx,ibndy,ibndz,bndx1,bndx2,bndx3,bndx4,gval,idirectalg,cgbp)

	      else

        call solvermodule(solverModuleFile,idpos,i_db,i_sf1,i_sf2,i_iqpos,i_qval,icount2a,icount2b,
     &	icount1a,icount1b,spec,i_phimap,i_phimap1,i_phimap2,i_phimap3,
     &	i_ibndx,i_ibndy,i_ibndz,i_bndx1,i_bndx2,i_bndx3,i_bndx4,i_gval,
     &  ibctype,isolvarch,rionst,igrid,scale,nlit,res2)
		 
	      end if
	   else

	      i_qmap1=memalloc(i_qmap1,realsiz,nhgp)
	      i_qmap2=memalloc(i_qmap2,realsiz,nhgp)
	      i_debmap1=memalloc(i_debmap1,realsiz,nhgp)
	      i_debmap2=memalloc(i_debmap2,realsiz,nhgp)
	
	      if (noit.gt.50) then
		 nlit=noit/2
		 write(6,*)'Re-estimated iterations now :',nlit
	      end if
	      qfact=abs(qnet)*float(ibc)/ibnum

	      if (isolvarch==0) then


		 call nitit(idpos,db,sf1,sf2,iqpos,qval,icount2a,icount2b,
     &	icount1a,icount1b,spec,nsp,phimap,phimap1,phimap2,phimap3,
     &  idebmap,ibndx,ibndy,ibndz,bndx1,bndx2,bndx3,bndx4,
     &	qmap1,qmap2,debmap1,debmap2,idirectalg,qfact)

		 i_qmap1=memalloc(i_qmap1,0,0)
		 i_qmap2=memalloc(i_qmap2,0,0)
		 i_debmap1=memalloc(i_debmap1,0,0)
		 i_debmap2=memalloc(i_debmap2,0,0)
        
	      else

      	call nlsolvermodule(solverModuleFile,idpos,i_db,i_sf1,i_sf2,i_iqpos,i_qval,icount2a,icount2b,
     &	icount1a,icount1b,spec,i_phimap,i_phimap1,i_phimap2,i_phimap3,
     &	i_ibndx,i_ibndy,i_ibndz,i_bndx1,i_bndx2,i_bndx3,i_bndx4,i_gval,
     &  ibctype,isolvarch,rionst,igrid,scale,nlit,res2,i_idebmap,imanual,nnit,ival(1),ival(2),ival2(1),
     &  ival2(2),conc(1),conc(2),i_qmap1,i_qmap2,qfact,epsout,deblen,relpar)

	endif

        i_qmap1=memalloc(i_qmap1,0,0)
        i_qmap2=memalloc(i_qmap2,0,0)
        i_debmap1=memalloc(i_debmap1,0,0)
        i_debmap2=memalloc(i_debmap2,0,0)

	end if
c       
	   i_bndx1= memalloc(i_bndx1,0,0)
	   i_bndx2= memalloc(i_bndx2,0,0)
	   i_bndx3= memalloc(i_bndx3,0,0)
	   i_bndx4= memalloc(i_bndx4,0,0)
	   
	   i_ibndx= memalloc(i_ibndx,0,0)
	   i_ibndy= memalloc(i_ibndy,0,0)
	   i_ibndz= memalloc(i_ibndz,0,0)
c       
	   i_phimap1= memalloc(i_phimap1,0,0)
	   i_phimap2= memalloc(i_phimap2,0,0)
	   i_phimap3= memalloc(i_phimap3,0,0)
c
	   i_db= memalloc(i_db,0,0)
	   i_idpos= memalloc(i_idpos,0,0)
	   i_sf1= memalloc(i_sf1,0,0)
	   i_sf2= memalloc(i_sf2,0,0)
c       b++++++++++++++++++++++++++++++++
c       ++++++now encalc takes nmedia++++
	   call encalc(icount1b,nqass,natom,icount2b,nmedia,nqgrd,nobject,esolvation)
	   i_chrgv2=memalloc(i_chrgv2,0,0)
c       e++++++++++++++++++++++++++++++++++++++++++++
	   if(isite) then
	      i=0
	      if(isitsf) i=1
c	      call wrtsit(nqass,icount2b,atpos,chrgv4,rad3,natom,
c     &  ibnum,nmedia,nobject,i)
c       i_iatmmed=memalloc(i_iatmmed,0,0)
	      do i=1, natom
		 xo(1)=atpos(i*3-2)
		 xo(2)=atpos(i*3-1)
		 xo(3)=atpos(i*3)
		 call ctog(xo,xn)
		 call phintp(xn,vphi)
c		 PRINT *, 'phi ', vphi
		 sitpot(i)=vphi
	      enddo
	   endif

c	finish=cputime(start)
c b++++++++++++++++++write potential map for the GUI
c       call wrtphiForGUI
c e+++++++++++++++++++++++++++++++++++++++++++++++++++++++
c
c	if phiwrt set true then write potential map
c       
	   phiwrt=py_out_phi
	   if(phiwrt) call wrtphi
c b++++++++++++++++++++++++++++++++++++++++++++++++
c	if(idebwrt) call wrtdebyemap
c e+++++++++++++++++++++++++++++++++++++++++++++++
	i_idebmap= memalloc(i_idebmap,0,0)
	i_phimap= memalloc(i_phimap,0,0)
c
	   

	endif
	goto 998
c
c
 999	continue
	stop
c
 998	continue
c	PRINT *, 'esolvation', esolvation
c	esolvation=esolvation
c	PRINT *, 'esolvation', esolvation
	finish = cputime(start)
	finish = finish - start
	write(6,*)'  '
	write(6,*)'total cpu time was (sec) ',finish
	write(6,*)'  '
	call datime(day)
	write(6,*)'DelPhi exited at ',day(12:19)
	end
	
	
	
	
	
	

	
        
