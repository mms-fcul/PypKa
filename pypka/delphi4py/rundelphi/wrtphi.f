	subroutine wrtphi

	include "qdiffpar4.h"
	include "qlog.h"

	real phimap(igrid,igrid,igrid)
	real*4 phimap4(igrid*igrid*igrid)

	if (realsiz.ne.4.and.phifrm.ne.2) then
	   i=1
	   do iz=1,igrid
	      do iy=1,igrid
		 do ix=1,igrid
		    phimap4(i)=phimap(ix,iy,iz)
		    i=i+1
		 end do
	      end do
	   end do
	else if (realsiz.eq.4) then
	   i=1
	   do iz=1,igrid
	      do iy=1,igrid
		 do ix=1,igrid
		    phimap4(i)=phimap(ix,iy,iz)
		    if ((ix .eq. 46) .and. (iy .eq. 46) .and. (iz .eq. 91)) then
		       PRINT *, phimap4(i), phimap(ix,iy,iz), ix,iy,iz
		    else if ((ix .eq. 46) .and. (iy .eq. 46) .and. (iz .eq. 1)) then
		       PRINT *, phimap4(i), phimap(ix,iy,iz), ix,iy,iz
		    endif

		    i=i+1
		 end do
	      end do
	   end do	   	   
	end if
	

        return
        end
