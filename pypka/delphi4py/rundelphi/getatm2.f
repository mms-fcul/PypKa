	subroutine getatm(fname,nam1,ifrm,idfrm,iatinf,iatrad,iatcrg,atot,
     &  nmedia,nobject,ndistr,repsin,epkt,realsiz,ionlymol,pdbfrm)

	include "pointer.h"
        include "qdiffpar5.h"

	character*80  fname,line,fname2,asci
	character*6 head,cval1
	character*7 cval2
	character*8 cval3,cval4
	character*24 crdstr
	character*8 strtmp
	character*26 xfield
	character*12 xnum
	character*15 atminf
	character*10 cnum
	character*26 cap1,cap2
	real xo(3),radx,crgx
	logical ifrm,iatinf,iatrad,iatcrg
	

        real repsintmp,epkt,repsin
        integer iatmmed(natmax),nmedia,objecttype,realsiz,pdbfrm           
        integer tmpiatmmed(nobjectmax),nobject
        character*80 datadistr(ndistrmax)
	character*50 datastr
        logical ifirst,ionlymol

	integer atot

        ndistr=0
        ionlymol=.true.

	idfrm=0
	cnum="0123456789"

	ifrm=.true.

	i=0

 204	format(a80)
 205	format(3F8.3)
 206	format(F7.3)
 2062	format(F7.3)
 2061	format(F6.2)
 2063	format(F8.4)
 207	format(F8.3)
 208	format(I3)

	call form(fname,nam1,ifrm)

	if(ifrm) then
	   open(9,file=fname(:nam1),form='formatted')
	   read(9,204,end=4000,err=4000)line
	   if(index(line,"DELPHI").ne.0) then
	      write(6,*) "Reading a Delphi-style pdb file"
	      read(9,204,end=4000,err=4000)line
	      j=index(line,"=")
	      read(line(j+1:80),'(i5)') idfrm		 
	   end if
	   
	   nmedia=1
	   nobject=1
	   objecttype=0
	   imedianumb=1
	   repsintmp=repsin
	   tmpiatmmed(nobject)=imedianumb
	   
	   atot=0
 502	   continue
	   head = line(1:6)
	   call up(head,6)
	   if((head.ne.'ATOM  ').and.(head.ne.'HETATM').and.(head.ne.
     &     'OBJECT').and.(head.ne.'MEDIA ').and.(head.ne.'CRGDST')) then
	      read(9,204,end=2000,err=4000)line
	      goto 502
	   end if
	   atot=atot+1
	   crdstr=line(31:54)
	   read(crdstr,205)xo
	   read(9,204,end=2000,err=4000)line
	   goto 502
	   
4000	continue
	write(6,*) "an error occurred in reading this formatted file"

2000	continue 
	write(6,*) "number of atoms read in = ",atot
	end if

c        i_datadistr=memalloc(i_datadistr,4,80*ndistr)

999	continue
c
	return
	end
