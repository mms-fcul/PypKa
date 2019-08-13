	subroutine chkcrg(natom,i_atinf,chrgv4,resnummax,isitsf)
	
      pointer (i_atinf,atinf)
      real chrgv4(natom)
	character*15 atinf(natom)
      character*80 line
      character*4 err,strng,resnme,resnum,rsnm
	integer resnummax,tmp
	logical isitsf

      crgsm=0.0
      rsnm='    '
      resnme='    '
c	write(6,*)' '

c b+++++++++++++++++++++++++++++++++++++++++++++++
      if (isitsf) then
	  do i=1,natom
	    crg=chrgv4(i)
          resnum=atinf(i)(12:15)
          strng=atinf(i)(7:10)
          if(resnum.ne.rsnm.or.strng.ne.resnme)then
            if(abs(crgsm).gt.1.0e-4)then
              eror=abs(crgsm)-1.0
              if(abs(eror).gt.1.0e-4)
     &	      write(6,10)resnme,rsnm,crgsm
            endif
            resnme=strng
            rsnm=resnum
            crgsm=crg
          else
            crgsm=crgsm+crg
          endif 
c specific section to find maximum residue number
          read(resnum,'(i4)')itmp
	    if(itmp.gt.resnummax) resnummax=itmp
	  end do
	else
c e++++++++++++++++++++++++++++++++++++++++++++
	  do i=1,natom
	    crg=chrgv4(i)
          resnum=atinf(i)(12:15)
          strng=atinf(i)(7:10)
          if(resnum.ne.rsnm.or.strng.ne.resnme)then
            if(abs(crgsm).gt.1.0e-4)then
              eror=abs(crgsm)-1.0
              if(abs(eror).gt.1.0e-4)
     &	      write(6,10)resnme,rsnm,crgsm
            endif
            resnme=strng
            rsnm=resnum
            crgsm=crg
          else
            crgsm=crgsm+crg
          endif 
	  end do
	endif

      if(abs(crgsm).gt.1.0e-4)then
        eror=abs(crgsm)-1.0
        if(abs(eror).gt.1.0e-4)
     &    write(6,10)resnme,rsnm,crgsm
      endif
c	write(6,*)' '
10	format('!!! WARNING: ',2a4,' has a net charge of ',f8.4)
      end
