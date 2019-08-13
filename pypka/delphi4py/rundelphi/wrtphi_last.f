	subroutine wrtphi
c
cNB uses phimap3 as a temp storage for phimap in case
c want to write not potentials bu salt concentraions
c
	include "qdiffpar4.h"
	include "qlog.h"
        parameter (mgrid=65)
c
	real phimap(igrid,igrid,igrid)
	character*80 filnam
	character*10 nxtlbl
	character*20 uplbl
	character*16 botlbl
c b+++++++++++++++++++++++++++for back compatibility with single precision phi readers
	real*4 phimap4(1),phimap3(ngp)
	real*4 scalesingle,scalesingle1,oldmidsingle(3),oldmidsingle1(3)
	real*4 minim,maxim,somma,average
	real*4 origin(3),coeff,stepsize
	integer igrid1      

	scalesingle=scale
	oldmidsingle(1)=oldmid(1)
	oldmidsingle(2)=oldmid(2)
	oldmidsingle(3)=oldmid(3)


      if (realsiz.ne.4.and.phifrm.ne.2) then
c	   i_phimap4=memalloc(i_phimap4,4,igrid*igrid*igrid)
	   i=1
         do iz=1,igrid
         do iy=1,igrid
         do ix=1,igrid
            phimap4(i)=phimap(ix,iy,iz)
	      i=i+1
         end do
         end do
         end do
      end if
c	  write(6,*) phimap(5,5,5)
c	  write(6,*) phimap4(5+(5-1)*igrid+(5-1)*igrid*igrid)

	
c b++++++++++++++walter++++write potential map for the GUI
c       open(52,FILE="phimap",form='unformatted')
c       write(52)phimap
c       close(52)
c e+++++++++++++++++++++++++++++++++++++++++++++++++++++++
        return
        end


