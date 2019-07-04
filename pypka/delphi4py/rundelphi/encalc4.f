	subroutine encalc(icount1b,nqass,natom,ibnum,nmedia,nqgrd,
     &  nobject,esolvation)
c
	include "qdiffpar4.h"
	include "qlog.h"
c
c       integer igridout !lilin test
	real*8 esolvation
        real phimap(igrid,igrid,igrid)
	dimension gchrg(icount1b)
	dimension rad3(natom),xn2(3,natom),scspos(3,ibnum)
	integer gchrgp(3,icount1b),ibgrd(3,ibnum),nqgrd
	logical ido
c comments on energies
c	logs= solvation energy
c	logc= coulombic energy
c	loga= analytic energy
c	logg= grid energy
c     logions= ionic contribution
c     lognl=ionic contribution in non linear case
c b++++++++++++++
        integer igridout
        real*8 ergs,ergas,ergnl,ergc,ergg
c        real sout(4,100) !lilin test
c	real medeps(0:nmedia)

c        do i=1,10
c               print *,'lilin_encalc',i,sout(1,i)
c        enddo

c
	if(inrgwrt) open(42,file=nrgnam(:nrgfrm))

	if(loga)  then
c	  call anagrd(icount1b,epsin*epkt,erga,scale)
        erga=0.0
	  write(6,*) 'analytic grid energy is no longer available'
	  stop
	  if(inrgwrt) write(42,*) 'analytic grid energy is',erga,' kt'
	end if

	if(logg) then
	  ergg=0.0
          limx1=2+bufz(1,1)
	  limx2=igrid-1-bufz(2,1)
	  limy1=2+bufz(1,2)
	  limy2=igrid-1-bufz(2,2)
	  limz1=2+bufz(1,3)
	  limz2=igrid-1-bufz(2,3)
	  do 587 i=1,icount1b
	    ix=gchrgp(1,i)
	    iy=gchrgp(2,i)
	    iz=gchrgp(3,i)
	    ido=.true.
	    if((ix.lt.limx1).or.(ix.gt.limx2)) ido=.false.
	    if((iy.lt.limy1).or.(iy.gt.limy2)) ido=.false.
	    if((iz.lt.limz1).or.(iz.gt.limz2)) ido=.false.
	    if(ido) ergg=ergg + phimap(ix,iy,iz)*gchrg(i)
587	  continue
	  ergg=ergg/2.0
c	  write(6,*) ' '
	  write(6,*) 'total grid energy          :     ',ergg,' kt'
	  if(inrgwrt) write(42,*)'total grid energy: ',ergg,' kt'
	end if
c
	if(logg.and.loga) then
	  write(6,*) 'difference energy, in kt, is',(ergg-erga)
	  write(6,*) 'difference energy, in kcals, is',(ergg-erga)*0.6
	end if
c
c b+++++++++++++++++++++++++++++++++++++w Oct 2000
      if (ibctyp.eq.5) then
	  write(6,*)
     &"WARNING!!!Not completely tested routine for polarized membrane!!"

	  if(lognl.or.logas) then
          write(6,*)
     &"This option is not yet working with fixed potential difference!"
	    ergas=0.0
          ergnl=0.0
	  end if

c	calcolo contributo facce "condensatore" assunte facce zeta!!
c     uso formula 0.5*(|q| media sulle due facce)*deltaV

	  deltaphi=0.0
	zfieldown=0.
	zfieldup=0.
	        open(52,file='fields.txt',form='formatted')
        do ix=2,igrid-1
	  do iy=2,igrid-1
	   phiup    = phimap(ix,iy,igrid)
	   phidown  = phimap(ix,iy,1)
	   deltaphi = deltaphi+ (phiup-phimap(ix,iy,igrid-1))
	zfieldup=zfieldup-(phiup-phimap(ix,iy,igrid-1))*scale
c	write(52,*)ix,iy,phimap(ix,iy,igrid-2)-(38.5-0.0025547)
c     &,phimap(ix,iy,2)-0.0025547
	write(52,*)ix,iy,phimap(ix,iy,igrid-1),phimap(ix,iy,igrid-2)

c	write(52,*)ix,iy,phiup,phidown
c         write(52,*)ix,iy,(phimap(ix,iy,igrid-1)-phiup)*scale,
c     &(phidown-phimap(ix,iy,2))*scale

c      averaging over upper and lower cube faces in order to cancel numerical errors
	   deltaphi = deltaphi- (phidown-phimap(ix,iy,2))
c        write(6,*)"phiup= ",phiup,"phidown= ",phidown
	zfieldown=zfieldown+(phidown-phimap(ix,iy,2))*scale
c	write(6,*)ix,iy,zfieldup,zfieldown,0.5*(zfieldup+zfieldown)
	  end do
	  end do
	  enface=deltaphi*vdropz*epsout*((igrid-1.)/(igrid-2.))**2/
     &(4.*fpi*scale)
c      ci sarebbe *((igrid-1)/(igrid-4))**2 ma ininfluente
        write(6,*)"Energy contribution from voltage drop=",enface,"kt"
c		write(6,*)deltaphi,vdropz,epsout

        close (52)
        open(52,file='potcen.txt',form='formatted')
c        write(6,*)"medeps",medeps(0)*epkt,medeps(1)*epkt
	  write(6,*)"fieldup medio: ",zfieldup/(igrid-2.)**2
	  write(6,*)"fieldown medio: ",zfieldown/(igrid-2.)**2
         do iz=1,igrid
           write(52,*)iz,phimap((igrid+1.)/2,(igrid+1.)/2,iz)
c           write(52,*)iz,phimap(10,10,iz)
          end do
        close (52)
 	end if   

c	else
	  if(irea.or.logs.or.lognl.or.logas.or.isen.or.isch) then
	
          ergs=0.0
          ergas=0.0
          ergnl=0.0
          ergest=0.0
c ergest=interaction energy of the solvent and the fixed charges
	    if(diff) ibc=0
	    iisitpot=0
	    if(isitpot) iisitpot=1
	    call react(nqass,icount1b,ibnum,ergs,ergas,natom,nmedia,
     &               nobject,iisitpot)
	    esolvation = ergs
	 end if
c
	  if(logc.and.(.not.logions.or..not.lognl)) then
          ergc=0.0
          if (logions) then
c         linear case
            ergest=0.0
c            print *,'lilin1',ergest
            call clbtot(nqass,ergest,ergc)
c            print *,'lilin2',ergest
            write(6,*)'solvent contribution to fixed charges'
            write(6,*)'respectively inside and outside the cube:',
     &ergest,'kt',ergestout,'kt'
            write(6,*)'total ionic direct contribution :',ergest+
     &ergestout,'kt'
          else
            if (nmedia.eq.1) then
 	        call clb(nqass,ergc)
              ergc=ergc/epsin
            else
              call clbmedia(nqass,ergc)
            end if
          end if
	    write(6,*) 'coulombic energy :              ',ergc,' kt'
	    if(inrgwrt) write(42,*) 'total coulombic energy:',ergc,' kt'
        end if

        if (lognl) then
         call nlener(ergnl,igridout)
         ergc=0.0
         ergest=0.0
c         print *,'lilin1',ergest
         if (logions) then
c                print *,'lilin2',ergest
          call clbnonl(nqass,ergc,ergest,igridout)
c          do i=1,10
c               print *,'lilin_encalc',i,sout(1,i)
c        enddo

c               print *,'lilin3',ergest
          write(6,*)'direct ionic contrib. inside the box:',ergest,' kt'
          write(6,*) 'coulombic energy:                     ',ergc,' kt'
          if(inrgwrt) write(42,*) 'total coulombic energy:',ergc,' kt'
         end if
        end if

        if (logs.and.logions) then
          write(6,*)'Energy arising from solvent and boundary pol.',
     &ergnl+ergs+ergest+ergestout,' kt'
        end if

        if (lognl.and.logg) then
          write(6,*)'Total non linear grid energy:',ergg+ergnl,' kt'
        end if

        ergtot=ergnl+ergc+ergs+ergest+ergestout
	  if (logs.or.logc) then
         write(6,*)'All required energy terms but grid and self_react.:'
     &,ergtot,'kt'
         if(inrgwrt) write(42,*)'total required energy (everything
     & calculated but grid and self_reaction energies: ',ergtot,'kt'
	  end if
c	end if
c e+++++++++++++++++++++++++++++++++++++
c
	if(logas.and.loga.and.logg) then
	write(6,*) "excess grid energy= ",ergg-ergas-erga
	end if
c
	finish=cputime(start)
	if(verbose)write(6,*) 'energy calculations done at',finish
c
	if(inrgwrt) close(42)
c
	return
	end
