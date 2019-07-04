	subroutine setrc(natom,*,nmedia,nobject,ndistr)
c
	include "qdiffpar4.h"
	include "qlog.h"
c
	character*80 line,filnam
	character*6 head,str6
	character*5 atnum 
	character*24 crdstr
	character*15 radstr
	character*16 crgstr
	character*5 atmstr
	dimension atpos(3*1),rad3(1),chrgv4(1),xo(3)
	logical ifrm,iatinf,iatrad,iatcrg
	character*15 atinf(1)


c b++++++++++++++++++++
        integer nmedia,nobject,ndistr
c e++++++++++++++++++++
 204	format(a80)
 205	format(3f8.3)
 206	format(F6.2,F7.3)
c
	write(6,*) "assigning charges and radii..."
	write(6,*) " "
	natom=0
c
c begin read
c       b++++++++++++++++++++

c	call getatm(pdbnam,pdblen,ifrm,idfrm,iatinf,iatrad,iatcrg,natom,
c     &  nmedia,nobject,ndistr,repsin,epkt,realsiz,ionlymol,pdbfrm)
	if (natom.gt.natmax) then 
	   write(6,*)'number of atom coordinates read',natom,'is larger'
        write(6,*)'than natmax parameter, it must be increased!!'
	  stop
	end if
c e++++++++++++++++++++
c
        print*,'in setrc, iatrad=', iatrad
	if(.not.iatrad) then
c
	do i=1,natom
c
	  atm = atinf(i)(1:5)
	  res = atinf(i)(7:9)
	  rnum = atinf(i)(12:15)
	  chn = atinf(i)(11:11)
	  
	  PRINT *, '1', atinf(i)(1:5), atinf(i)(7:9), atinf(i)(12:15), atinf(i)(11:11)
c
	  call up(atm,6)
	  call elb(atm,6)
	  call up(res,3)
	  call elb(res,3)
	  call up(rnum,4)
	  call elb(rnum,4)
	  call up(chn,1)
	  call elb(chn,1)
c
c assign radius, searching for decreasingly specific specification
c ending with generic atom type
c note all atoms must have an assignment
c
	call radass(atm,res,rnum,chn,rad,norad)
c
	if(norad.eq.1) then
	rad=0.0
c	if(atm(1:1).ne.'H')
c     &	write(6,'(''!!! WARNING: no radius record for '',a15)')atinf(i)

      write(6,'(''!!! WARNING: no radius record for '',a15)')atinf(i)
      print*,'atminfo=',atinf(i)
      print*,'rad3[',i,']=', rad3(i);
c
c need stop here
c
      elseif(rad.lt.1.e-6.and.(atm(1:1).ne.'H'.and.atm(2:2).ne.'H'))
     &  then 
        write(6,'(''!!! WARNING: radius of heavy atom'',a15,''
     &	is set to zero'')')atinf(i)
	end if
c
c store rad,xn in rad3,  for later use
c
	rad3(i)=rad
c
c scale and assign charge to grid
c
	  call crgass(atm,res,rnum,chn,chrgv)
	  chrgv4(i)=chrgv
c
	end do
c
c write record to new coordinate file if required, with
c occupancy and temperature factor fields replaced by radius and charge
c
	end if
c
c check charge assignments (Sri April 2, 1996)
	call chkcrg(natom,i_atinf,chrgv4,resnummax,isitsf)
	do i=1,natom
	   PRINT *, '2', atinf(i)(1:5), atinf(i)(7:9), atinf(i)(12:15), atinf(i)(11:11)
	end do

	PRINT *, ipdbwrt, iatout
c unformatted write...
	if(ipdbwrt) then
	   open(20,file=updbnam(:updblen),form='unformatted')
	   do i=1,natom
	      xo(1)=atpos(3*i-2)
	      xo(2)=atpos(3*i-1)
	      xo(3)=atpos(3*i)
	      rad=rad3(i)
	      chrgv=chrgv4(i)
	      write(20) xo,rad,chrgv
	   end do
	   close(20)
	end if
c       
	if(iatout) then
c       
	   open(19,file=mpdbnam(:mpdblen))
c       
	   filnam = ' '
	   inquire(19,name = filnam)
	   write(6,*)
     &	  'atomic coordinates, charges and radii written to file'
	   write(6,*)filnam
	   write(6,*)'   '
	   if (mpdbfrm.eq.0) then
	      write(19,*)'DELPHI PDB FILE'
	      write(19,*)'FORMAT = 1', mpdbfrm
	      write(19,*)'HEADER output from qdiff'
	      write(19,*)'HEADER atom radii in columns 55-60'
	      write(19,*)'HEADER atom charges in columns 61-67'
	   end if
	   if (mpdbfrm.eq.1) then
	      write(19,*)'DELPHI PDB FILE'
	      write(19,*)'FORMAT = PQR'
	      write(19,*)'HEADER output from qdiff'
	      write(19,*)'HEADER atom charges in columns 56-61'
	      write(19,*)'HEADER atom radii   in columns 63-68'
	   end if
	   if (mpdbfrm.eq.40) then !LinLi: 4 digit precision
	      write(19,*)'DELPHI PDB FILE'
	      write(19,*)'FORMAT = 1', mpdbfrm
	      write(19,*)'4 digits precison'
	      write(19,*)'HEADER output from qdiff'
	      write(19,*)'HEADER atom radii in columns 55-60'
	      write(19,*)'HEADER atom charges in columns 61-67'
	   end if
	   if (mpdbfrm.eq.41) then !LinLi: 4 digit precision
	      write(19,*)'DELPHI PDB FILE'
	      write(19,*)'FORMAT = PQR'
	      write(19,*)'4 digits precison'
	      write(19,*)'HEADER output from qdiff'
	      write(19,*)'HEADER atom charges in columns 56-61'
	      write(19,*)'HEADER atom radii   in columns 63-68'
	   end if

c
	   line=' '
	   if (mpdbfrm.eq.0) then
	      line(1:11)="ATOM       "
	      do i=1,natom
		 rad=rad3(i)
		 chrgv=chrgv4(i)
		 xo(1)=atpos(3*i-2)
		 xo(2)=atpos(3*i-1)
		 xo(3)=atpos(3*i)
		 line(12:26)=atinf(i)
		 write(radstr,206)rad,chrgv
		 write(crdstr,205)xo
		 line(55:67) = radstr
		 line(31:54) = crdstr
		 line(27:30)='    '
		 write(19,204)line
	      end do
	   end if
c+++++++++++++LinLi++++++++++
	   if (mpdbfrm.eq.40) then
	      line(1:11)="ATOM       "
	      do i=1,natom
		 rad=rad3(i)
		 chrgv=chrgv4(i)
		 xo(1)=atpos(3*i-2)
		 xo(2)=atpos(3*i-1)
		 xo(3)=atpos(3*i)
		 line(12:26)=atinf(i)
		 write(radstr,"(F7.4,F8.4)")rad,chrgv
		 write(crdstr,205)xo
		 line(55:70) = radstr
		 line(31:54) = crdstr
		 line(27:30)='    '
		 write(19,204)line
	      end do
	   end if

c ++++++++++++++WWW++++++++++
	   if (mpdbfrm.eq.1) then
	      line(1:6)="ATOM  "
	      do i=1,natom
		 write(atmstr,208)i
		 line(7:11)=atmstr
		 rad=rad3(i)
		 chrgv=chrgv4(i)
		 xo(1)=atpos(3*i-2)
		 xo(2)=atpos(3*i-1)
		 xo(3)=atpos(3*i)
		 line(12:26)=atinf(i)
		 line(22:22)=' '
		 write(crdstr,205)xo
		 line(31:54) = crdstr
c       blank after coordinates
		 line(55:55)=' '
c       write(crgstr,207)chrgv,rad
		 write(str6,209)chrgv
		 line(56:61) = str6
		 line(62:62)=' '
		 write(str6,209)rad
		 line(63:68) = str6
c       line(55:70) = crgstr
		 line(27:30)='    '
c       clear the chain identifier
		 write(19,204)line
	      end do
	   end if

	

c ++++++++++++++LinLi++++++++++
	   if (mpdbfrm.eq.41) then
	      line(1:6)="ATOM  "
	      do i=1,natom
		 write(atmstr,208)i
		 line(7:11)=atmstr
		 rad=rad3(i)
		 chrgv=chrgv4(i)
		 xo(1)=atpos(3*i-2)
		 xo(2)=atpos(3*i-1)
		 xo(3)=atpos(3*i)
		 line(12:26)=atinf(i)
		 line(22:22)=' '
		 write(crdstr,205)xo
		 line(31:54) = crdstr
c       blank after coordinates
c       line(55:55)=' ' !Lin
		 write(crgstr,'(F8.4,F7.4)')chrgv,rad !Lin
c       write(str6,209)chrgv
		 line(55:70) = crgstr
		 line(27:30)='    '
c       clear the chain identifier
		 write(19,204)line
	      end do
	   end if
	   
 207	   format(F8.4,F8.3)
 209	   format(F6.3)
 208	   format(I5)
c       ++++++++++++++WWW++++++++++
	
	   close(19)
	end if
c
 	do i=1,natom
	   PRINT *, '4', atinf(i)(1:5), atinf(i)(7:9), atinf(i)(12:15), atinf(i)(11:11)
	end do
 	do i=1,natom
	   PRINT *, '4', atinf(i)
	end do
	if((natom.eq.0).and.(nobject.eq.0)) goto 903
	return
903	write(6,*) "exiting due to non-existence of atom file 
     &  nor object data"
	return 1
	end



