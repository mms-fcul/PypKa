	subroutine rforce(afield,natom,scspos,scrg,
     &  atsurf,ibnum,atmcrg,xn1,chgpos,nqass,crgatn,nobject)
c
	dimension afield(3,natom+nobject),atmcrg(4,nqass),
     &	chgpos(3,nqass),xn1(3,natom)
	dimension scspos(3,ibnum),scrg(ibnum),sfield(3)
	integer atsurf(ibnum),scale
	integer crgatn(nqass)
c
	do i=1,natom+nobject
	  afield(1,i)=0.0
	  afield(2,i)=0.0
	  afield(3,i)=0.0
	end do

	do i=1,ibnum
	  sfield(1)=0.0
	  sfield(2)=0.0
	  sfield(3)=0.0
c
	  x=scspos(1,i)
	  y=scspos(2,i)
	  z=scspos(3,i)
	  sc=scrg(i)
c erased some comments (Walter)
c calculate total field on this surface element due to ALL charges
	  do j=1,nqass
	    dist1=(x-chgpos(1,j))
	    dist2=(y-chgpos(2,j))
	    dist3=(z-chgpos(3,j))
c
	    dist=dist1**2+dist2**2+dist3**2
	    sdist=sqrt(dist)
	    temp=atmcrg(4,j)/(dist*sdist)
c
	    sfield(1)=sfield(1)+temp*dist1
	    sfield(2)=sfield(2)+temp*dist2
	    sfield(3)=sfield(3)+temp*dist3
c add negative of this to the charged atom
	    iat=crgatn(j)
c b+++++++++++++++++++++++++++++
          if (iat.lt.0) go to 10
c e++++++++++++++++++++++++++++++
	    if(iat.eq.0)then
	      write(6,*)'problems with crgatn '
	      stop
	    endif

	    afield(1,iat)=afield(1,iat)-sc*temp*dist1
	    afield(2,iat)=afield(2,iat)-sc*temp*dist2
	    afield(3,iat)=afield(3,iat)-sc*temp*dist3

10      continue
	  end do
c
	  sfield(1)=sfield(1)*scrg(i)
	  sfield(2)=sfield(2)*scrg(i)
	  sfield(3)=sfield(3)*scrg(i)
c
	  j=atsurf(i)
c
	  afield(1,j)=afield(1,j)+sfield(1)
	  afield(2,j)=afield(2,j)+sfield(2)
	  afield(3,j)=afield(3,j)+sfield(3)
	end do
	return
	end


c b+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
	subroutine rforceeps1(afield,natom,scspos,scrg,scale,
     &atsurf,ibnum,atmcrg,xn1,chgpos,nqass,crgatn,nobject,scsnor,oldmid)
c
	dimension afield(3,natom+nobject),atmcrg(4,nqass),
     &	chgpos(3,nqass),xn1(3,natom),oldmid(3)
	dimension scspos(3,ibnum),scsnor(3,ibnum),scrg(ibnum),sfield(3)
	integer atsurf(ibnum),p
	integer crgatn(nqass)
	real fact
c       era un 4 bytes
c
	do i=1,natom+nobject
	  afield(1,i)=0.0
	  afield(2,i)=0.0
	  afield(3,i)=0.0
	end do
c      fact=-2.*3.1416*80/(79.)
      fact=-2.*3.14159265359*80/(79.) !LinLi test
	sum=0.0
	cont=0
	zax=0.0
	ipmax=0
	realds=0.0
	rmyds=0.0
	do p=1,ibnum
	  sfield(1)=0.0
	  sfield(2)=0.0
	  sfield(3)=0.0
c
	  x=scspos(1,p)
	  y=scspos(2,p)
	  z=scspos(3,p)
	ix=nint( (x - oldmid(1))*scale )
	iy=nint( (y - oldmid(2))*scale )
	iz=nint( (z - oldmid(3))*scale )
	xx=float(ix)/scale + oldmid(1)
	yy=float(iy)/scale + oldmid(2)
	zz=float(iz)/scale + oldmid(3)

	  sc=0.5*scrg(p)
      
	fact1=0.8/(scale**2)-
     &     ((x-xx)*scsnor(1,p)+(y-yy)*scsnor(2,p)+(z-zz)*scsnor(3,p))**2
c      deltas=3.1416*fact1
      deltas=3.14159265359*fact1 !LinLi test
      sigmap=scrg(p)/deltas
	realds=realds+scrg(p)/0.0393
	rmyds=rmyds+deltas
c	if (abs(sigmap-0.0393).gt.zax) then
c	  ipmax=p
c	  zax=abs(sigmap-0.0393)
c	end if
	trullo=fact*sigmap*sigmap*deltas
      if(p.eq.400) then
	write(6,*)"normale",scsnor(1,p),scsnor(2,p),scsnor(3,p)
	write(*,*)"P",x-xx,y-yy,z-zz
	write(*,*)x,y,z
      write(*,*)"area",scrg(p)/0.0393,scrg(p)*scale*scale/0.0393
	end if
	  j=atsurf(p)
c erased some comments (Walter)
	  do iat=1,natom
          if(iat.eq.j) then
	cont=cont+1
  	      afield(1,iat)=afield(1,iat)+trullo
	      afield(2,iat)=afield(2,iat)
	      afield(3,iat)=afield(3,iat)
c	 sum=sum+scrg(p)*scrg(p)*scale*scale*scale*scale/(0.0393*0.0393)
c	 write(6,*)"dp",scrg(p),0.0393,scrg(p)*0.8/scale/scale/0.0393,
c     &0.0393*0.8/scale/scale*ibnum
c  	      afield(1,iat)=afield(1,iat)+fact1*scsnor(1,p)
c	      afield(2,iat)=afield(2,iat)+fact1*scsnor(2,p)
c	      afield(3,iat)=afield(3,iat)+fact1*scsnor(3,p)
	    end if
	  end do
c calculate total field on this surface element due to ALL charges
	  do i=1,-nqass
	    dist1=(x-chgpos(1,i))
	    dist2=(y-chgpos(2,i))
	    dist3=(z-chgpos(3,i))
c
	    dist=dist1**2+dist2**2+dist3**2
	    sdist=sqrt(dist)
	    temp=atmcrg(4,i)/(dist*sdist)
c
	    sfield(1)=sfield(1)+temp*dist1
	    sfield(2)=sfield(2)+temp*dist2
	    sfield(3)=sfield(3)+temp*dist3
c add negative of this to the charged atom
	    iat=crgatn(i)
c b+++++++++++++++++++++++++++++
          if (iat.lt.0) go to 10
c e++++++++++++++++++++++++++++++
	    if(iat.eq.0)then
	      write(6,*)'problems with crgatn '
	      stop
	    endif
          if(iat.eq.j) then
  	      afield(1,iat)=afield(1,iat)-sc*temp*dist1+fact1*scsnor(1,p)
	      afield(2,iat)=afield(2,iat)-sc*temp*dist2+fact1*scsnor(2,p)
	      afield(3,iat)=afield(3,iat)-sc*temp*dist3+fact1*scsnor(3,p)
	    end if
	  end do
10    continue
c
	  sfield(1)=sfield(1)*sc
	  sfield(2)=sfield(2)*sc
	  sfield(3)=sfield(3)*sc
c
c	  afield(1,j)=sfield(1)+afield(1,j)
c	  afield(2,j)=sfield(2)+afield(2,j)
c	  afield(3,j)=sfield(3)+afield(3,j)
	end do
c      write(*,*)"supcalc",realds/(16*3.1416),"mia",rmyds/(16*3.1416)
      write(*,*)"supcalc",realds/(16*3.14159265359),"mia",
     &rmyds/(16*3.14159265359) !LinLi test
	return
	end
c ***********************************************************
	subroutine rforcenew(natom,nobject,nmedia,ibnum,nqass)

	include "qdiffpar4.h"
	include "qlog.h"
c
	integer natom,ibnum,nqass,nobject,nmedia
c     atmcrg (1..3,nqass) atomic charge pos. in grid units, like xn2(1..3,natom)
	dimension rfield(3,natom+nobject),atmcrg(4,nqass),xn1(3,natom)
c     chgpos real charge pos in Amstrong
c     scspos, surface  charge position in Amstrong
	dimension scspos(3,ibnum),schrg(ibnum),chgpos(3,nqass)
	integer iatmmed(natom+nobject),crgatn(nqass),imed,q,p
	real polariz(natom),eps,medeps(0:nmedia),atmeps(nqass)
	real cost1, cost2,fact1,fact2,fact3,chrgv4(natom)
      real fx,fy,fz,ffx,ffy,ffz,qq,qi,xi,yi,zi,xq,yq,zq,xp,yp,zp
      real riq,riq2,rip,rip2
	real distiqx,distiqy,distiqz,distipx,distipy,distipz
c       erano un 4 bytes

      if (.not.ionlymol) then
	  write(6,*)"Not yet ready to give forcefield in case of objects"
	  stop
	end if

 	do 100 i=1,natom
	  fx=0.
	  fy=0.
	  fz=0.

        qi=chrgv4(i)
	  imed=iatmmed(i)
	  eps=medeps(imed)*epkt
	  cost1=3./(eps+2.)
	  cost2=polariz(i)*epkt*cost1
        xi=xn1(1,i)
        yi=xn1(2,i)
        zi=xn1(3,i)

c ciclo sulle cariche vere
	  do 200 q=1,nqass
	    if (crgatn(q).eq.i) goto 200
          qq=atmcrg(4,q)/atmeps(q)
          xq=chgpos(1,q)
          yq=chgpos(2,q)
          zq=chgpos(3,q)

	    distiqx=xi-xq
	    distiqy=yi-yq
	    distiqz=zi-zq

	    riq2=distiqx**2+distiqy**2+distiqz**2
		riq=sqrt(riq2)
          fact=qq/(riq*riq2)
          fact1=-2.*fact*fact

c     termine elettroforetico
          fx=fx+qi*fact*distiqx
          fy=fy+qi*fact*distiqy
          fz=fz+qi*fact*distiqz
c     termine dielettroforetico q=p
          ffx=fact1*distiqx
          ffy=fact1*distiqy
          ffz=fact1*distiqz
c     termine dielettroforetico q<>p col simmetrico
	    do 300 p=1,q-1
	      if (crgatn(p).eq.i) goto 300
            qp=atmcrg(4,p)/atmeps(p)
            xp=chgpos(1,p)
            yp=chgpos(2,p)
            zp=chgpos(3,p)

	      distipx=xi-xp
	      distipy=yi-yp
	      distipz=zi-zp

	      rip2=distipx**2+distipy**2+distipz**2
		  rip=sqrt(rip2)
	      fact2=qp*fact/(rip*rip2)
            fact3=(distiqx*distipx+
     &            distiqy*distipy+distiqz*distipz)/(riq2*rip2)

            ffx=ffx+fact2*(xp-xq-3.*(riq2*distipx+rip2*distiqx)*fact3)
            ffy=ffy+fact2*(yp-yq-3.*(riq2*distipy+rip2*distiqy)*fact3)
            ffz=ffz+fact2*(zp-zq-3.*(riq2*distipz+rip2*distiqz)*fact3)
300       continue
          fx=fx+ffx*cost2
          fy=fy+ffy*cost2
          fz=fz+ffz*cost2
200     continue

c ciclo sulle cariche di polarizzazione
	  do 210 q=1,ibnum
          qq=schrg(q)*epkt
          xq=scspos(1,q)
          yq=scspos(2,q)
          zq=scspos(3,q)

	    distiqx=xi-xq
	    distiqy=yi-yq
	    distiqz=zi-zq

	    riq2=distiqx**2+distiqy**2+distiqz**2
		riq=sqrt(riq2)
          fact=qq/(riq*riq2)
          fact1=-2.*fact*fact

c     termine elettroforetico
          fx=fx+qi*fact*distiqx
          fy=fy+qi*fact*distiqy
          fz=fz+qi*fact*distiqz
c     termine dielettroforetico p=q
          ffx=fact1*distiqx
          ffy=fact1*distiqy
          ffz=fact1*distiqz
c     termine dielettroforetico q<>p col simmetrico
c     prima parte: p scorre le cariche vere
	    do 305 p=1,nqass
	      if (crgatn(p).eq.i) goto 305
            qp=atmcrg(4,p)/atmeps(p)
            xp=chgpos(1,p)
            yp=chgpos(2,p)
            zp=chgpos(3,p)

	      distipx=xi-xp
	      distipy=yi-yp
	      distipz=zi-zp

	      rip2=distipx**2+distipy**2+distipz**2
		  rip=sqrt(rip2)
	      fact2=qp*fact/(rip*rip2)
            fact3=(distiqx*distipx+
     &            distiqy*distipy+distiqz*distipz)/(riq2*rip2)

            ffx=ffx+fact2*(xp-xq-3.*(riq2*distipx+rip2*distiqx)*fact3)
            ffy=ffy+fact2*(yp-yq-3.*(riq2*distipy+rip2*distiqy)*fact3)
            ffz=ffz+fact2*(zp-zq-3.*(riq2*distipz+rip2*distiqz)*fact3)
305       continue

c     seconda parte, p scorre le cariche di polarizzazione
	    do 310 p=1,q-1
            qp=schrg(p)*epkt
            xp=scspos(1,p)
            yp=scspos(2,p)
            zp=scspos(3,p)

	      distipx=xi-xp
	      distipy=yi-yp
	      distipz=zi-zp

	      rip2=distipx**2+distipy**2+distipz**2
		  rip=sqrt(rip2)
	      fact2=qp*fact/(rip*rip2)
            fact3=(distiqx*distipx+
     &            distiqy*distipy+distiqz*distipz)/(riq2*rip2)
            ffx=ffx+fact2*(xp-xq-3.*(riq2*distipx+rip2*distiqx)*fact3)
            ffy=ffy+fact2*(yp-yq-3.*(riq2*distipy+rip2*distiqy)*fact3)
            ffz=ffz+fact2*(zp-zq-3.*(riq2*distipz+rip2*distiqz)*fact3)
310       continue

          fx=fx+ffx*cost2
          fy=fy+ffy*cost2
          fz=fz+ffz*cost2
210     continue
	  
	  rfield(1,i)=fx*cost1
	  rfield(2,i)=fy*cost1
	  rfield(3,i)=fz*cost1
100   continue

	return
	end
c e+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
