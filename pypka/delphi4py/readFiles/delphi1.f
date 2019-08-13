c###############################################################################
	subroutine delphi(acent_xdim, py_igrid, py_scale, py_repsin, 
     & py_repsout, py_acent, py_in_pdb, py_in_crg, py_in_siz, natom, 
     & nobject, py_i_atpos, py_i_rad3, py_i_chrgv4, py_i_atinf, 
     & py_i_medeps, py_i_iatmmed, py_i_dataobject, py_rmaxdim)

c       PARSING ARGUMENTS FROM PYTHON
c	python stores floats as doubles
c       py_variables are not used directly but rather reassigned
	integer*8 acent_xdim
	integer*8 py_igrid
	real*8  py_scale
	real*8  py_repsin
	real*8  py_repsout
	real*8  py_acent(acent_xdim)
        character*(*) py_in_pdb, py_in_crg, py_in_siz
	integer*8 natom
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
	real*8 py_rmaxdim

c       END PARSING ARGUMENTS FROM PYTHON

Cf2py   intent(in) py_igrid
Cf2py   intent(in) py_scale
Cf2py   intent(in) py_repsin
Cf2py   intent(in) py_repsout
Cf2py   intent(in) py_acent
Cf2py   intent(in) py_in_pdb
Cf2py   intent(in) py_in_crg
Cf2py   intent(in) py_in_siz
Cf2py   intent(in) natom
Cf2py   intent(in) nobject
Cf2py   intent(in) py_i_atpos
Cf2py   intent(in) py_i_rad3
Cf2py   intent(in) py_i_chrgv4
Cf2py   intent(in) py_i_iatmmed
Cf2py   intent(in) py_i_dataobject
Cf2py   intent(in,out) py_rmaxdim
c       b++++++++++++++++++++++++++++++++++++++++++++
cpgi$ IF DEFINED (PC)
#ifdef PC
cDEC$ IF DEFINED (PC)
        use dflib
        use dfport
cDEC$ END IF
#endif
c       e++++++++++++++++++++++++++++++++++++++++++++

	include 'qdiffpar4.h'
	include 'qlog.h'

c------------------------------------------------------------------
c       MODIFIED EXISTING ARRAYS
	real         atpos(natom*3)
	real         rad3(natom)
	real         chrgv4(natom)
	character*15 atinf(natom)
	real         medeps(0:nmediamax)
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
c       b+++++++++++++++++++
        integer crgatn(1),nqgrdtonqass(1),extracrg
        integer nmedia,ndistr,idimtmp
	logical uniformdiel
        integer idirectalg,numbmol,icount1a,icount1b
        real atmeps(1),qfact
	dimension ibndx(1),ibndy(1),ibndz(1),bndx1(1),bndx2(1),bndx3(1)
        dimension bndx4(1),qmap1(1),qmap2(1),debmap1(1),debmap(2)
c       character*96 dataobject(nobjectmax,2)
c       character*80 datadistr(ndistrmax)
c       character*80 strtmp
c
c       pointer based arrays
        real limobject(1),cgbp(5,1),chgpos(3,1)
c       e+++++++++++++++++++

	integer iepsmp(1)
        logical*1 idebmap(1)
        real phimap(1),atmforce(3,1)
        real phimap1(1),phimap2(1),phimap3(1)
	integer neps(1),keps(1)


	icount1a=0
	icount1b=0

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
c       b_p++++++++++++++++++++++++++++++++++++++
c       
c	i=iargc()
c	argnam=" "
c	if(i.gt.0)then
c	   call getarg(1,argnam)
c	   call namlen(argnam,arglen)
c	end if
c

c	i=1
	
c	reseting variables

	i_atpos=py_i_atpos
	i_rad3=py_i_rad3
	i_chrgv4=py_i_chrgv4
	i_atinf=py_i_atinf
	i_medeps=py_i_medeps
	i_iatmmed=py_i_iatmmed
	i_dataobject=py_i_dataobject
	
	call defprm()

c	icount1a=0
c	icount1b=0
c	i_atpos=atom
c	i_rad3=radius
c	i_chrgv4=charge
c	i_atinf=atominfo

c 	lookup table number of the medium -> epsilon of the medium
c	i_medeps=mediumeps
c 	lookup table atom -> number of the medium
c	i_iatmmed=atommedia
	
c	i_polariz=polar

c	epsmap is the grid, it has the atom numbers, and the media numbers


c	statement parsing
c
	igrid     = py_igrid   
	scale     = py_scale   
	repsin    = py_repsin  
	repsout   = py_repsout 
	radprb(1) = py_radprb

c	
c	function parsing
c
c       grid offset
	acent(1) = py_acent(1)
	acent(2) = py_acent(2)
	acent(3) = py_acent(3)	
	iacent=.true.
	itest2=.true.
	
c
c       file names parsing
c
	ias=1
	pdbnam=py_in_pdb
	call namlen(pdbnam,pdblen)

	crgnam=py_in_crg
	call namlen(crgnam,crglen)

	siznam=py_in_siz
	call namlen(siznam,sizlen)
c
c       print parsed arguments
c
	call rdprm
		
c       e_p++++++++++++++++++++++++++++++++++++++
	
c       b++++++++++++++++++++++++++++++++++++++
c       create suitable pdb file, if necessary
c       if (icreapdb) call creapdb(repsout,numbmol)
c       e++++++++++++++++++++++++++++++++++++++
c       
c       set hash tables from size and charge files
c
	call rdhrad(*999)
	if(isolv) call rdhcrg(*999)
c
c       read in pdbfile and, if necesssary, assign charges/radii
c
c       b++++++++++++++++++++++++++++++
        uniformdiel=.true.
c       i_iatmmed=memalloc(i_iatmmed,4,natmax)
c       e+++++++++++++++++++++++++++++

	call setrc(natom,*999,nmedia,nobject,ndistr)

c	PRINT *, 'natom -> ', natom
c	do i=1, natom
c	   PRINT *, atinf(i)
c	end do

	
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
        medeps(0)=epsout	
c       e+++++++++++++++++++++++++++++
c

c	PRINT *, medeps(0), medeps(1), medeps(2)
	
c	finish=cputime(start)
c	write(6,*) 'time to read in and/or assign rad/chrg=',finish

	
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


	py_rmaxdim = rmaxdim

	goto 998
c       
 999	continue
	stop
c
 998	continue

c	finish = cputime(start)
c	finish = finish - start
c	write(6,*)'  '
c	write(6,*)'total cpu time was (sec) ',finish
c	write(6,*)'  '
c       call datime(day)
c       write(6,*)'DelPhi exited at ',day(12:19)
	
	end
