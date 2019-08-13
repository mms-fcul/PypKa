	subroutine crgarr(ncrgmx,cqplus,cqmin,xn1,igrid,natom,
     &	nqass,nqgrd,qmin,qnet,qplus,nmedia,ndistr,scale,cmid,nobject,
     &  radpolext,extracrg,realsiz,verbose)

	include "pointer.h"
c
	dimension xn2(3,natom),xn1(3,natom)
        integer extracrg
	dimension chrgv2(natom+extracrg)
c b++++++++++++
      integer nqgrdtonqass(natom+extracrg),ncrgmx,realsiz
      real atmeps(natom+extracrg),atmcrg(4,1),chgpos(3,ncrgmx)
      dimension iatmmed(natom+nobject),iepsmp(igrid,igrid,igrid,3)
      real medeps(0:nmedia),rx,ry,rz,cmid(3)
      real rad3(natom)
      integer imed,jx,jy,jz,epsdim
      logical ipassed,verbose
	character*15 atinf(natom)
c e++++++++++++
	dimension xo(3),chrgv4(natom),cqplus(3),cqmin(3)
	integer crgatn(natom+extracrg)
c
	ic1=0
c
c atmcrg contains grid positions of all charges AND the charge in the 4th field
c atmeps6 contains 6*epsilon/epkt as a function of ic2-th charge internal
c atmeps6  to the grid NO LONGER USED
c nqgrdtonqass maps ic2 to ic1
c atmeps contains epsilon/epkt as a funcion of ic1-th general charge
c
c b+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        if(ncrgmx.le.natom) then
          write(*,*)'Increase ncrgmx, at least as the number of atoms'
          stop
        end if 
        i_atmcrg=memalloc(i_atmcrg,realsiz,4*ncrgmx)
        epsdim=nobject+natom+2
c find charge moments for dipole approximation
c
	qnet=0.0
	qplus=0.0
	qmin=0.0
	do 2106 ix=1,3
	  cqplus(ix)=0.0
	  cqmin(ix)=0.0
2106	continue


	do 3104 ix=1,natom
	  if(abs(chrgv4(ix)).gt.1.e-6) then
	    ic1=ic1+1
	    chrg=chrgv4(ix)
	    atmcrg(1,ic1)=xn2(1,ix)
	    atmcrg(2,ic1)=xn2(2,ix)
	    atmcrg(3,ic1)=xn2(3,ix)
	    chgpos(1,ic1)=xn1(1,ix)
	    chgpos(2,ic1)=xn1(2,ix)
	    chgpos(3,ic1)=xn1(3,ix)
	    atmcrg(4,ic1)=chrg
	    crgatn(ic1)=ix

          qnet=qnet + chrg
          if(chrg.gt.0.) then
            qplus=qplus + chrg
	      cqplus(1)=cqplus(1) + chrg*atmcrg(1,ic1)
	      cqplus(2)=cqplus(2) + chrg*atmcrg(2,ic1)
	      cqplus(3)=cqplus(3) + chrg*atmcrg(3,ic1)
	    else
	      qmin=qmin + chrg
	      cqmin(1) =cqmin(1)  + chrg*atmcrg(1,ic1)
	      cqmin(2) =cqmin(2)  + chrg*atmcrg(2,ic1)
	      cqmin(3) =cqmin(3)  + chrg*atmcrg(3,ic1)
	    end if

	  end if
3104	continue
       
      if(verbose)
     &   write(6,*)"number of charges coming from molecules ",ic1

c insert charges from charge distributions
      if (ndistr.gt.0) call distrTOpoint(igrid,scale,ncrgmx,ndistr,ic1
     &  ,cmid,natom,qplus,qmin,cqplus,cqmin,qnet)

c debug++++++++++++++++++++++++WWW
      if (.false.) then
        open(52,file='Charges.txt',form='formatted')
        do iiii=1,ic1
          write(52,*) iiii,chgpos(1,iiii),chgpos(2,iiii),chgpos(3,iiii)
        end do
        close (52)
      end if
c end debug+++++++++++++++++++
c e+++++++++++++++
c
c assign charges for boundary conditions
c
c ic1 = number of charges

c divide by charge totals
c
	if(qplus.gt.1.e-6) then
	  do 2110 k = 1,3
	    cqplus(k) = cqplus(k)/qplus
2110	  continue
	end if
	if(abs(qmin).gt.1.e-6) then
	  do 2111 k = 1,3
	    cqmin(k) = cqmin(k)/qmin
2111	  continue
	end if
c
c select those charges which will be charging the grid
c
	nqass=ic1
c b+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        if(ncrgmx.le.nqass) then
          write(*,*)'Increase ncrgmx, to be at least equal to:',nqass
          stop
        end if 
        i_atmcrg=memalloc(i_atmcrg,realsiz,4*nqass)
c e+++++++++++++++++++++++++++++++++++++++++++++++++++++++
	rgrid=igrid
	ic2=0

	do 2102 ix=1,nqass
c b++++++++++
c crgatn(crg number)=atom number or natom+objectnumber or -distr.number
          ii=crgatn(ix)

          if(ii.lt.0) then
c now we have to consider charge distributions
c in this case the distribution is not linked to any object
c (jx,jy,jz)=coordinates of closest grid point to charge
c (rx,ry,rz)=coordinates of the charge relatives to the current grid point
              jx=int(atmcrg(1,ix)+0.5)
              jy=int(atmcrg(2,ix)+0.5)
              jz=int(atmcrg(3,ix)+0.5)
              rx=atmcrg(1,ix)-jx
              ry=atmcrg(2,ix)-jy
              rz=atmcrg(3,ix)-jz
              if (rz.gt.rx) then
                if (rz.gt.-rx) then
                  if (rz.gt.ry) then
                    if (rz.gt.-ry) then
                      imed=iepsmp(jx,jy,jz,3)
                    else
                      imed=iepsmp(jx,jy-1,jz,2)
                    end if
                  else
                    imed=iepsmp(jx,jy,jz,2)
                  end if
                else
                  if (ry.gt.rx) then
                    if (ry.gt.-rx) then
                      imed=iepsmp(jx,jy,jz,2)
                    else
                      imed=iepsmp(jx-1,jy,jz,1)
                    end if
                  else
                    imed=iepsmp(jx,jy-1,jz,2)
                  end if
                end if
              else
                if (rz.gt.-rx) then
                  if(ry.gt.rx) then
                    imed=iepsmp(jx,jy,jz,2)
                  else
                    if (ry.gt.-rx) then
                      imed=iepsmp(jx,jy,jz,1)
                    else
                      imed=iepsmp(jx,jy-1,jz,2)
                    end if
                  end if
                else
                  if (rz.gt.ry) then
                    imed=iepsmp(jx,jy-1,jz,2)
                  else
                    if (rz.gt.-ry) then
                      imed=iepsmp(jx,jy,jz,2)
                    else
                      imed=iepsmp(jx,jy,jz-1,3)
                    end if
                  end if
                end if
              end if
            imed=imed/epsdim
          else
            imed=iatmmed(ii)
          end if
          atmeps(ix)=medeps(imed)
c e+++++++++++++++++++++++++++++++++

	  if((atmcrg(1,ix).gt.1.).and.(atmcrg(1,ix).lt.rgrid)) then
	  if((atmcrg(2,ix).gt.1.).and.(atmcrg(2,ix).lt.rgrid)) then
	  if((atmcrg(3,ix).gt.1.).and.(atmcrg(3,ix).lt.rgrid)) then
	    ic2=ic2+1
            chrgv2(4*ic2-3)=atmcrg(1,ix)
	    chrgv2(4*ic2-2)=atmcrg(2,ix)
	    chrgv2(4*ic2-1)=atmcrg(3,ix)
	    chrgv2(4*ic2)=atmcrg(4,ix)
c b+++++++++++++++++++++++++++++++++
            nqgrdtonqass(ic2)=ix
c           atmeps6(ic2)=6.0*atmeps(ix)
c e+++++++++++++++++++++++++++++++++
          end if
          end if
	  end if
2102	continue
	nqgrd=ic2

        ipassed=.false.
        do i=1,nqass
          ii=crgatn(i)
          if(ii.gt.0.and.ii.le.natom) then
            if (rad3(ii).le.0.) then
             ipassed=.true.
             write(6,'(I4,A16,'' is charged! Radius moved from zero to''
     &            ,f5.1)')ii,atinf(ii),radpolext
             rad3(ii)=radpolext
            end if
          end if
        end do
	  if(ipassed) write(6,*)'BE CAREFUL!! A WRONG ASSIGNMENT FOR THE 
     &RADIUS MIGHT LEAD TO INACCURATE REACTION FIELD ENERGY !!!'

	return
	end
