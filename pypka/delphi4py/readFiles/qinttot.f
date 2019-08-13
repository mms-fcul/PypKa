	subroutine defprm
c
c logical parameters used in qdiff, mainly from parameter file
c
	include "qlog.h"
c
c     data conversion to kT/charge at 25 celsius
      data fpi /12.566370614359/
c
	toplbl="qdiffxas: qdiffxs4 with an improved surfacing routine"
	isch=.false.
	isen=.false.
	iacs=.false.
	irea=.false.
	iautocon=.true.
	iper(1)=.false.
	iper(2)=.false.
	iper(3)=.false.
        iper(4)=.false.
        iper(5)=.false.
        iper(6)=.false.
	iconc=.false.
	ibios=.false.
	isite=.false.
	iatout=.false.
	diff=.false.
	isph=.false.
	phiwrt=.false.
	logs=.false.
	logc=.false.
	loga=.false.
	logg=.false.
c b++++++++++++++++++++++++++++++w Oct 2000
        logions=.false.
c flag for energy calculation of contribution from the solvent
        lognl=.false.
c flag for non linear energy calculation
        icreapdb=.false.
c flag for automatic insertion of objects
        imanual=.false.
c flag for manual assignment of relaxation parameterr
        verbose=.true.
c flag for removing part of the standard output
        idebwrt=.false.
c flag for saving idebmap
c e++++++++++++++++++++++++++++++++++++++++
	logas=.false.
	ipdbrd=.false.
	epswrt=.false.
	iacent=.false.
	ipdbwrt=.false.
	ifrcwrt=.false.
	ifrcrd=.false.
	ipoten=.false.
	igraph=.false.
	imem=.false.
	ihome=.false.
	ibufz=.false.
	inrgwrt=.false.
	iwgcrg=.false.
	iuspec=.false.
	icheb=.false.
	idbwrt=.false.
	ixphird=.false.
	iyphird=.false.
	izphird=.false.
	ixphiwt=.false.
	iyphiwt=.false.
	izphiwt=.false.
	isita=.false.
	isitq=.false.
	isitp=.false.
c b+++++++++++++++++++++++
        isitap=.false.
	  isitmd=.false.
	  isittf=.false.
	  isitsf=.false.
c e+++++++++++++++++++++++
	isitf=.false.
	isitr=.false.
	isitc=.false.
	isitx=.false.
	isiti=.false.
	iself=.false.
	isitrf=.false.
	isitcf=.false.
	isitt=.false.
	isolv=.true.
	isrf=.false.
c
	icon1=10
	icon2=1
	epsnam="fort.17"
	debnam="debmap.dat"
	phinam="fort.14"
	srfnam="grasp.srf"
	frcnam="fort.16"
	mpdbnam="fort.19"
	updbnam="fort.20"
	ufrcnam="fort.21"
	centnam="fort.15"
	pdbnam="fort.13"
	crgnam="fort.12"
	siznam="fort.11"
	phiinam="fort.18"
	frcinam="fort.15"
	prmnam="fort.10"
	scrgnam="scrg.dat"
	nrgnam="energy.dat"
	gcrgnam="crg.dat"
	dbnam="db.dat"
	xphinam="fort.31"
	yphinam="fort.32"
	zphinam="fort.33"
	xphonam="fort.34"
	yphonam="fort.35"
	zphonam="fort.36"
        meshnam="triangulatedSurf.off"
c
	epslen=7
	debnamlen=10
	philen=7
	srflen=9
	frclen=7
	mpdblen=7
	updblen=7
	ufrclen=7
	centlen=7
	pdblen=7
	crglen=7
	sizlen=7
	phiilen=7
	frcilen=7
	scrglen=8
	nrglen=10
	prmlen=7
	gcrglen=7
	dblen=6
	xphilen=7
	yphilen=7
	zphilen=7
	xpholen=7
	ypholen=7
	zpholen=7
        meshlen=20
c
	phifrm=0
	epsfrm=0
	frcfrm=0
	mpdbfrm=0
	updbfrm=0
	ufrcfrm=0
	pdbfrm=0
	crgfrm=0
	sizfrm=0
	prmfrm=0
	phiifrm=0
	frcifrm=0
	scrgfrm=0
	nrgfrm=0
	gcrgfrm=0
        radprb(1)=1.4
	scale=10000.
	exrad=2.0
	perfil=10000.
c b++++++++++++++++++
        radpolext=1.0
        radprb(2)=-1.0
	conc(1)=0.0
        conc(2)=0.0
        rionst=0.0
        relpar=1.0 
c ival are the valencies of salts
        ival(1)=1
        ival(2)=1
        ival2(1)=0
        ival2(2)=0
        res1=0.0
        res2=0.0
        res5=0.0
        vdropx=0.0
        vdropy=0.0
        vdropz=0.0
        atompotdist=0.5
        temperature=297.3342119
        realsiz=4
#ifdef PC
cDEC$ IF DEFINED (PC)
cDEC$ IF DEFINED (DP)
cDEC$ ELSE
      realsiz=0
cDEC$ END IF
cDEC$ END IF
#else	
#ifdef DP
        realsiz=realsiz+4
#endif
#endif
          print*,'realsiz == ', realsiz
	  resnummax=0
c e++++++++++++++++++
	repsout=80
	epsout=80
	repsin=2
	epsin=2
	gten=0.0
	uspec=0.9975
	offset(1)=0.
	offset(2)=0.
	offset(3)=0.
	acent(1)=0.0
	acent(2)=0.0
	acent(3)=0.0
	do j=1,3
	bufz(1,j)=0
	bufz(2,j)=0
	end do
	igrid=0
	nlit=0
	nnit=0
	ibctyp=2

c jbc +++++++++++
        isolvarch=0
c jbc +++++++++++

	return
	end
c ***********************************************************
	subroutine rdflnm(j1,j2,direc,fnam,flen)
c
	character*100 direc,calph
	integer un(10),unum,flen,six,five
	character*10 cnumb
	character*80 fnam
c
	cnumb="1234567890"
	calph=" "
	calph(1:40)="ABCDEFGHIJKLMNOPQRSTUVWXYZ.:_-+=!@#$^123"
	calph(41:80)="4567890abcdefghijklmnopqrstuvwxyz|\/?><;"
      fnam=" "
	flen=0
c
	if(j1.ne.0) then
	j=j1
	k=j+4
	n=1
150	k=k+1
	l=index(cnumb,direc(k:k))
	if(l.ne.0) then
	un(n)=l
	n=n+1
	goto 150
	end if
	unum=0
	fnam(1:5)="fort."
	do i=1,n-1
	fnam(5+i:5+i)=cnumb(un(i):un(i))
	end do
	flen=4+n
c	write(6,*) "(filename)unit number,length= ",fnam,flen
	end if
c
c j2=position of first letter in file descriptor, i.e. the letter
c f of "file"
	if(j2.ne.0) then
	j=j2
	k=j+5
	if(index(calph,direc(k:k)).ne.0) k=k-1
	n=1
	m=k+1
160	k=k+1
	l=index(calph,direc(k:k))
	if(l.ne.0) then
	n=n+1
	goto 160
	end if
	fnam(1:n-1)=direc(m:j+4+n)
	flen=n-1
c	write(6,*) "filename,length= ",fnam(1:n-1),flen
	end if
c
	return
	end
c *********************************************************************
	subroutine rdprm
c
	include "qlog.h"
c
	character*80 line,asci,enc,filnam
	character*10 cnumb
	character*400 mline
	dimension lin1(20),itype(20)
	logical iext
c b+++++++++++++++++++++++++++++++
        real cz1,cz2,z1p,z1m,z2p,z2m
c e++++++++++++++++++++++++++++++
c
c read in run parameters
c
204	format(a80)
	asci="1234567890 .-+#,$asdfghjklzxcvbnmqwertyuio
     &pASDFGHJKLZXCVBNMQWERTYUIOP)(}{][/"
c     &pASDFGHJKLZXCVBNMQWERTYUIOP)(}{][/\"
	cnumb="0123456789"
c

	if((repsin.lt.0.).or.(repsout.lt.0.)) then
	repsin=abs(repsin)
	repsout=abs(repsout)
	diff=.true.
	end if

c        dfact = 3.047*sqrt(repsout/80.)
      dfact=0.01990076478*sqrt(temperature*repsout)
c b++++++++++++++++++++++
c define a number that indicates whether or not there is some salt 
        z1p=ival(1)
        z1m=ival(2)
        z2p=ival2(1)
        z2m=ival2(2)
c now cz1 and cz2 are concentration of positive ion !!!
        cz1=conc(1)*z1m
        cz2=conc(2)*z2m
        rionst=(cz1*z1p*(z1p+z1m)+cz2*z2p*(z2p+z2m))/2.
c coefficients in Taylor series of the charge concentration
c  apart from n! (order >=1)
c        chi1=-2.*rionst
c        chi2=cz1*z1p*(z1p**2-z1m**2)+cz2*z2p*(z2p**2-z2m**2)
c        chi2=chi2/2.
c        chi3=cz1*z1p*(z1p**3+z1m**3)+cz2*z2p*(z2p**3+z2m**3)
c        chi3=-chi3/6.
c        chi4=cz1*z1p*(z1p**4-z1m**4)+cz2*z2p*(z2p**4-z2m**4)
c        chi4=chi4/24.
c        chi5=cz1*z1p*(z1p**5+z1m**5)+cz2*z2p*(z2p**5+z2m**5)
c        chi5=-chi5/120.
c 2012-04-24 chuan Correct coefficients in Taylor series. NOT in 
c compact form in order to compare with eqn. 12 in the paper of
c
c Rocchia, W.; Alexov, E.; Honig, B. "Extending the applicability of 
c the nonlinear Poisson-Boltzmann equation: Multiple dielectric 
c constants and multivalent ions" J Phys.Chem. B 105, 6507-6514 (2001)
        chi1=-2.*rionst
        chi2=cz1*z1p*z1p**2-cz1*z1m*z1m**2+cz2*z2p*z2p**2-cz2*z2m*z2m**2
        chi2=chi2/2.
        chi3=cz1*z1p*z1p**3+cz1*z1m*z1m**3+cz2*z2p*z2p**3+cz2*z2m*z2m**3
        chi3=-chi3/6.
        chi4=cz1*z1p*z1p**4-cz1*z1m*z1m**4+cz2*z2p*z2p**4-cz2*z2m*z2m**4
        chi4=chi4/24.
        chi5=cz1*z1p*z1p**5+cz1*z1m*z1m**5+cz2*z2p*z2p**5+cz2*z2m*z2m**5
        chi5=-chi5/120.

c convert ionic strength to debye length
c
	if(rionst.gt.1.e-6) then
          deblen = dfact/sqrt(rionst) 
          if (nnit.gt.0) lognl=.true.
	else
          logions=.false.
	  deblen = 1.e6
	end if
        
c       write(6,*)'non linear energy:', lognl
c e++++++++++++++++++++++
c
c test for unformatted pdb and frc files
c
	if(.not.ipdbrd) then
 	open(13,file=pdbnam(:pdblen),form='formatted')
        read(13,204,iostat=n)line
	ias=0
	do 600 i=1,80
	if(index(asci,line(i:i)).eq.0) ias=ias+1
600	continue
	if(ias.gt.10) ipdbrd=.true.
 	close(13)
	end if
c
	if(ifrcwrt) then
 	open(15,form='formatted')
        read(13,204,iostat=n)line
	ias=0
	do 610 i=1,80
	if(index(asci,line(i:i)).eq.0) ias=ias+1
610	continue
	if(ias.gt.10) ifrcrd=.true.
 	close(15)
	end if
c
c epkt assignment as a function of temperature
c	epkt=166804.4928/temperature
        epkt=167100.9162872952/temperature
c set epsin and epsout (=epkt adjusted dielectrics such that
c all distances are in angstroms, charges in e)
c
	epsin = repsin/epkt
	epsout = repsout/epkt
c
c go back to main
c
	return
	end
c ***************************************************************
