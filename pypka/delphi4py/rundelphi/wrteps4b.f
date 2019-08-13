      subroutine wrteps(imaxwrd,epsdim,nmedia,idimtmp)
c	kim sharp/9 feb 88
c	reformats epsilon array epsmap
c	to compact array (MKG format)
c	old format: epsmap(65,65,65,3) real*4
c	new format: neps(5,65,65,3) integer*2
c	where first index 1-65 now compressed
c	into 1-5 plus offset into 16 bit words
c	compact format also contains oldmid, the
c	center of the protein in real coordinates
c	compaction is effected by storing
c	real eps (which take values of 0. and 1.)
c	as bits in a 16 bit word
c	access is via pointers idimx and ioffset
c	thus x arrary indices of reps 0-15 -> word 1
c	16-31 -> word 2 etc
c	
	include 'qdiffpar4.h'
	include 'qlog.h'

c	parameter  (imaxwrd = igrid/16 + 1)

	integer iepsmp(igrid,igrid,igrid,3),idimtmp,epsdim,nmedia
	integer*4 tmpmap(idimtmp,idimtmp,idimtmp,3)
	real medeps(0:nmedia)
	logical*1 idebmap(igrid,igrid,igrid)



c attenzione!!! le dimensioni non sono quelle che sembrano! int*2!!!
	dimension neps(imaxwrd,igrid,igrid,3)	
	dimension keps(imaxwrd,igrid,igrid)	




c	compact fine epsilon map
	dimension idimx(ngrid),ip2(ngrid)
c	array of pointers to words
	dimension ioffset(ngrid)		
c	array of pointers to bit offsets
	integer*2  i,j,ioffset
c	integer*2  neps,i,j ,ioffset,keps
c 	integer  neps,ioffset
	character*80 filnam
c
c	integer epsmap(65,65,65)
c---------------------------------------------------------------------
        print *,"imaxwrd,igrid",imaxwrd,igrid
      if (epsfrm.eq.0) then
      write(6,*)' setting up pointers...'
	do 9000 ix = 1, igrid
	  idimx(ix) = ix/16 + 1
	  ioffset(ix) = mod(ix,16)
	  ip2(ix)=2**(mod(ix,16))
9000	continue
c	  ip2(15)=-2**15
c	  ip2(31)=-2**15
c	  ip2(47)=-2**15
c	  ip2(63)=-2**15
c	  ip2(79)=-2**15
c	  ip2(95)=-2**15
c	  ip2(111)=-2**15
c	  ip2(127)=-2**15
	do id=15,igrid,16
	ip2(id)=-2**15
	end do
      write(6,*)' clearing bits...'
      do 9001 idir = 1, 3
        do 9002 iz=1,igrid
          do 9003 iy=1,igrid
		do 9004 ix=1,imaxwrd
              neps(ix,iy,iz,idir) = 0
              keps(ix,iy,iz) = 0
9004        continue
9003	    continue
9002	  continue
9001	continue

      write(6,*)' generating compact fine epsilon array...'
        do 9006 iz=1,igrid
          do 9008 ix=1,igrid
	      i=idimx(ix)
            do 9007 iy=1,igrid
	        j1=0
        	  j2=0
        	  j3=0
           	  k1=0
c divide solvente da non solvente
	        if((iepsmp(ix,iy,iz,1)/epsdim).ne.0) j1=ip2(ix)
	        if((iepsmp(ix,iy,iz,2)/epsdim).ne.0) j2=ip2(ix)
	        if((iepsmp(ix,iy,iz,3)/epsdim).ne.0) j3=ip2(ix)
	        if(idebmap(ix,iy,iz)) k1=ip2(ix)
	        neps(i,iy,iz,1)=neps(i,iy,iz,1)+j1
	        neps(i,iy,iz,2)=neps(i,iy,iz,2)+j2
	        neps(i,iy,iz,3)=neps(i,iy,iz,3)+j3
	       keps(i,iy,iz)=keps(i,iy,iz)+k1
9007	      continue
9008	    continue
9006	  continue

	kmap = 1

      write(6,*)' writing to compact epsilon file'
	open(17,file=epsnam(:epslen),form='unformatted')
	filnam = ' '
	inquire(17,name = filnam)
	imap = 0
	write (17) kmap, scale, oldmid
	write (17) neps
	write (17) keps
	close (17)

c     end if of epsfrm 0, that is compact epswrite
      end if

c b++++++Aug 2011++++++++++++++++++++++++++++++++++++
      if (epsfrm.eq.1) then
c converting epsmap in a map containing only media numbers
        do iz=1,igrid
          do ix=1,igrid
            do iy=1,igrid
	        tmpmap(ix,iy,iz,1)=iepsmp(ix,iy,iz,1)/epsdim
              tmpmap(ix,iy,iz,2)=iepsmp(ix,iy,iz,2)/epsdim
	        tmpmap(ix,iy,iz,3)=iepsmp(ix,iy,iz,3)/epsdim
            end do
	    end do
	  end do

        write(6,*)' writing extended epsilon file'
	  open(17,file=epsnam(:epslen),form='unformatted')
	  filnam = ' '
	  inquire(17,name = filnam)
	  write (17) igrid, scale, oldmid
        write (17) nmedia, medeps
	  write (17) tmpmap
	  close (17)

c     end if of epsfrm 1, that is non compact epswrite
      end if
c e++++++++++++++++++++++++++++++++++++++++++


	write(6,*)' '
	write(6,*)'dielectric map written to file'
	write(6,*)filnam
	write(6,*)' '

c
999	continue
c
	return
	end

c b++++++Aug 2011++++++++++++++++++++++++++++++++++++
      subroutine wrtdebyemap
c	
	include 'qdiffpar4.h'
	include 'qlog.h'

	logical*1 idebmap(igrid,igrid,igrid)
      character*80 filnam

      write(6,*)' writing Debye boolean map'
	open(52,file=debnam(:debnamlen),form='unformatted')
	filnam = ' '
	inquire(52,name = filnam)
	write (52) igrid, scale, oldmid
	write (52) idebmap
	close (52)


	write(6,*)' '
	write(6,*)'Debye boolean map written to file'
	write(6,*)filnam
	write(6,*)' '
c
	return
	end

c e++++++++++++++++++++++++++++++++++++++++++
