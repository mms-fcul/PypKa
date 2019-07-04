	subroutine wrt(ipgn)
	character*24 day

c subroutine to handle log file writes
c
      if(ipgn.eq.1) then
         write(6,*)'   '
         write(6,*)' _____________DelPhi V. 5.1_Patched_____________   '
         write(6,*)'|                                                | '
         write(6,*)'| A program to solve the PB equation             | '
         write(6,*)'| in 3D, using non-linear form, incorporating    | '
         write(6,*)'| many dielectric regions, multisalt ionic       | '
         write(6,*)'| strength, different probe radii, periodic      | '
         write(6,*)'| and focussing boundary conditions, utilizing   | '
         write(6,*)'| stripped optimum successive over-relaxation    | '
         write(6,*)'| and an improved algorithm for mapping the      | '
         write(6,*)'| Mol. Surface to the finite-Difference grid     | '
         write(6,*)'|                                                | '
         write(6,*)'|    If there is any question, please go to:     | '
         write(6,*)'|       http://compbio.clemson.edu/forum/        | '
         write(6,*)'|     June 2012,by DelPhi Development Team       | '
         write(6,*)'|                                                | '
         write(6,*)'| This patched version exposes an interface to   | '
         write(6,*)'| CudaSolver and NanoShaper modules.             | '
         write(6,*)'| For more info and the proper references goto   | '
         write(6,*)'|       http://www.electrostaticszone.eu         | '
         write(6,*)'|_________________             __________________| '
         write(6,*)'              DelPhi V. 5.1 Patched                '
	     write(6,*)'  '
	call datime(day)
        write(6,*)' program started on ',day(1:10)//day(20:24)
        write(6,*)'             at ',day(12:19)
	end if
c
	return
	end
