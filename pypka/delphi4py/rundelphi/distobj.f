#ifdef PC
cDEC$ IF DEFINED (PC)
      recursive subroutine distobj(xp,dx,dy,dz,nobject,ii,prbrad,dist,
     &inout,zeta,axdist)
cDEC$ ELSE
#elif (IFC) || (AIX) || (GNU)
      recursive subroutine distobj(xp,dx,dy,dz,nobject,ii,prbrad,dist,
     &inout,zeta,axdist)
#else
      subroutine distobj(xp,dx,dy,dz,nobject,ii,prbrad,dist,
     &inout,zeta,axdist)
#endif
cDEC$ END IF

c this routine calculate the outgoing normal vector onto the surface
c  (possibly extended by prbrad) 
c be careful, dist is an oriented distance!!supposed to be positive 
c if the point is external to the object
c in this context, the vector d is (G-H)/dist
c where G is the point and H its projection onto the surface
c   (Walter, last update August 5, 2001) 

      include 'pointer.h'

      integer nobject,count,objecttype,ii,sign
      logical inout,imovepoint
      character*96 dataobject(nobject,2)
      character*96 strtmp,strtmp1
      real xb(3),xa(3),tmp,radius,vectz(3),vectx(3),vecty(3),
     &  alpha,xc(3),xd(3),modul2,modx,mody,prx,pry,prz,tmp1,tmp2
      real xp(3),prbrad,dist,modul,zeta,axdist,dot,talpha,pntsh,tet
	real extent
	real Pi,fact,tang,czeta,xxx
    
c      Pi=3.1415926536
      Pi=3.14159265359 !LinLi test
c       inout is option to calculate only dist and not dx,dy,dz
      strtmp=dataobject(ii,1)
      strtmp1=dataobject(ii,2)
      read(strtmp(8:10),'(I3)')objecttype
c     write(6,*)'objectype',strtmp

	extent=0.0

      if (objecttype.eq.1) then
c sphere
         read(strtmp(20:80),*)xb,radius
         dx=xp(1)-xb(1)
         dy=xp(2)-xb(2)
         dz=xp(3)-xb(3)
         tmp=dx*dx+dy*dy+dz*dz
         if (tmp.eq.0.) then 
c      point exactly in the center
           dist=-radius
	     dx=-1.0
           goto 999
         end if
         tmp1=sqrt(tmp)
         dist=tmp1-radius-extent
         if (.not.inout) then
            dx=dx/tmp1
            dy=dy/tmp1
            dz=dz/tmp1
         end if
      else
      if (objecttype.eq.2) then
c cylinder
         read(strtmp(20:80),*)xa,xb,radius
         read(strtmp1,'(5f8.3)')vectz,modul,modul2
          
         dx=xp(1)-xb(1)
         dy=xp(2)-xb(2)
         dz=xp(3)-xb(3)
c dot=(P-B)(A-B)
         dot=dx*vectz(1)+dy*vectz(2)+dz*vectz(3)
         zeta=dot/modul
c tmp=|P-B|**2
         tmp=(dx*dx+dy*dy+dz*dz)
c axdist is |P-B|sin(teta)
         tmp1=tmp-zeta*zeta
         if (tmp1.lt.0.0.and.tmp1.gt.-1.e-3)then
	     axdist=0.0
	   else
           axdist=sqrt(tmp1)
	   end if

         dist=-extent-zeta
         if (.not.inout) then
            if (dist.ge.0.) then
c            inferiore
               if (axdist.gt.radius+extent) then
c               laterale inferiore
                  tmp1=(radius+extent)/axdist
                  tmp2=(extent+zeta*tmp1)/modul
                  dx=dx*(1-tmp1)+vectz(1)*tmp2
                  dy=dy*(1-tmp1)+vectz(2)*tmp2
                  dz=dz*(1-tmp1)+vectz(3)*tmp2
                  dist=sqrt(dx*dx+dy*dy+dz*dz)
                  dx=dx/dist
                  dy=dy/dist
                  dz=dz/dist
               else
c               assiale inferiore
                  dx=-vectz(1)/modul
                  dy=-vectz(2)/modul
                  dz=-vectz(3)/modul
               end if
            else
               dist=zeta-(modul+extent)
               if (dist.ge.0.) then
c            superiore
                  if (axdist.gt.radius+extent) then
c            laterale superiore
                     tmp1=(radius+extent)/axdist
                     tmp2=(-extent+zeta*tmp1)/modul-1.
                     dx=dx*(1-tmp1)+vectz(1)*tmp2
                     dy=dy*(1-tmp1)+vectz(2)*tmp2
                     dz=dz*(1-tmp1)+vectz(3)*tmp2
                     dist=sqrt(dx*dx+dy*dy+dz*dz)
                     dx=dx/dist
                     dy=dy/dist
                     dz=dz/dist
                  else
c            assiale superiore
                     dx=vectz(1)/modul
                     dy=vectz(2)/modul
                     dz=vectz(3)/modul
                  end if
               else
c            z-centrale
                  tmp1=axdist-radius-extent
                  if(tmp1.lt.dist.or.tmp1.lt.-extent-zeta) then
c               vicini alle basi da dentro
                     if (-dist.le.zeta+extent) then
c                  base superiore
                        dx=vectz(1)/modul
                        dy=vectz(2)/modul
                        dz=vectz(3)/modul
                     else
c                  base inferiore
                        dist=-zeta-extent
                        dx=-vectz(1)/modul
                        dy=-vectz(2)/modul
                        dz=-vectz(3)/modul
                     end if
                  else
c               vicino ai lati (laterale), da dentro o da fuori
                    if (axdist.eq.0.) then
c                       point exactly onto the axis
                       dist=-radius
                       goto 999
                    end if
                    dist=tmp1
                    tmp2=zeta/modul                   
                    dx=(dx-vectz(1)*tmp2)/axdist
                    dy=(dy-vectz(2)*tmp2)/axdist
                    dz=(dz-vectz(3)*tmp2)/axdist
                  end if
               end if
            end if
         else
c -------------------------------------------------------------------
            if (dist.ge.0.) then
c            inferiore, assiale inf. gia' a posto
               if (axdist.gt.radius+extent) then
c               laterale inferiore
                  tmp1=(radius+extent)/axdist
                  tmp2=(extent+zeta*tmp1)/modul
                  dx=dx*(1-tmp1)+vectz(1)*tmp2
                  dy=dy*(1-tmp1)+vectz(2)*tmp2
                  dz=dz*(1-tmp1)+vectz(3)*tmp2
                  dist=sqrt(dx*dx+dy*dy+dz*dz)
               end if
            else
               dist=zeta-(modul+extent)
               if (dist.ge.0.) then
c            superiore
                  if (axdist.gt.radius+extent) then
c               laterale superiore, assiale gia' a posto
                     tmp1=(radius+extent)/axdist
                     tmp2=(-extent+zeta*tmp1)/modul-1.
                     dx=dx*(1-tmp1)+vectz(1)*tmp2
                     dy=dy*(1-tmp1)+vectz(2)*tmp2
                     dz=dz*(1-tmp1)+vectz(3)*tmp2
                     dist=sqrt(dx*dx+dy*dy+dz*dz)
                  end if
               else
c            z-centrale
                  tmp1=axdist-radius-extent
                  if(tmp1.lt.dist.or.tmp1.lt.-extent-zeta) then
c               vicini alle basi da dentro, base sup. gia' a posto
                     if (-dist.gt.zeta+extent) dist=-zeta-extent
c               base inferiore
                  else
c               vicino ai lati, da dentro o da fuori
                     dist=tmp1
                  end if
               end if
            end if
         end if
c -------------------------------------------------------------------
      else
c else di object type 2
      if (objecttype.eq.3) then
c cone, xb is the tip
         read(strtmp(20:80),*)xa,xb,alpha
c conversion degrees --> radiants
         alpha=alpha*Pi/180
         talpha=tan(alpha*.5)
         read(strtmp1,'(5f8.3)')vectz,modul,modul2
         radius=modul*talpha
	   calpha=cos(alpha*.5)
	   tang=tan((Pi-alpha)/4.)
	
          
         dx=xp(1)-xb(1)
         dy=xp(2)-xb(2)
         dz=xp(3)-xb(3)
c dot=(P-B)(A-B)
         dot=dx*vectz(1)+dy*vectz(2)+dz*vectz(3)
         zeta=dot/modul
c czeta= complement to zeta
	   czeta=modul-zeta
c tmp=|P-B|**2
         tmp=(dx*dx+dy*dy+dz*dz)
c axdist is |P-B|sin(teta)
         axdist=sqrt(tmp-zeta*zeta)
	   xxx=axdist-radius
         if (.not.inout) then
           if (xxx.le.0.) then
	       if (czeta.lt.0.) then
c exactly below the basis (zona 1)
	         dist=-czeta
               dx=vectz(1)/modul
               dy=vectz(2)/modul
               dz=vectz(3)/modul
	         goto 999
	       elseif(czeta.le.tang*(radius-axdist)) then
c internal, close to the basis (zona 5)
	         dist=-czeta
               dx=-vectz(1)/modul
               dy=-vectz(2)/modul
               dz=-vectz(3)/modul
	         goto 999
	       elseif(axdist.le.zeta*talpha) then
	         if (axdist.gt.0.)then
c internal, close to lateral (zona 6)
                 tmp1=-calpha/axdist
                 tmp2=calpha*(talpha+zeta/axdist)/modul
                 dx=dx*tmp1+vectz(1)*tmp2
                 dy=dy*tmp1+vectz(2)*tmp2
                 dz=dz*tmp1+vectz(3)*tmp2
                 dist=calpha*(axdist-zeta*talpha)
		       goto 999
	         else
                 if (zeta.gt.0.) then
c exactly on the axis
                   fact=zeta*talpha*calpha*calpha
	             tmp1=modul2-vectz(1)*vectz(1)
			     if(tmp1.gt.0.)then
c vectz not aligned with x axis 
                     tmp1=1./sqrt(tmp1)
	               tmp2=fact*(talpha-vectz(1)*tmp1)/modul
                     dx=tmp2*vectz(1)+fact*modul*tmp1
                     dy=tmp2*vectz(2)
                     dz=tmp2*vectz(3)
                     dist=-zeta*talpha*calpha
		           goto 999
	             else
c vectz aligned with x axis 
	               tmp2=fact*talpha/modul
                     dx=fact*talpha*vectz(1)/modul
                     dy=fact*talpha*vectz(2)/modul+fact
                     dz=fact*talpha*vectz(3)/modul
                     dist=-zeta*talpha*calpha
		           goto 999
	             end if
	           else
c exactly on tip (apex)
                   dx=vectz(1)/modul
                   dy=vectz(2)/modul
                   dz=vectz(3)/modul
                   dist=0.
		         goto 999
	           endif
	         end if
             end if
           elseif(czeta.lt.0.) then
c below basis, laterally (part of zona 2)
             tmp1=radius/axdist
c here axdist > 0 for sure, dist probably too
             tmp2=(zeta*tmp1)/modul-1.
             dx=dx*(1.-tmp1)+vectz(1)*tmp2
             dy=dy*(1.-tmp1)+vectz(2)*tmp2
             dz=dz*(1.-tmp1)+vectz(3)*tmp2
             dist=sqrt(dx*dx+dy*dy+dz*dz)
             dx=dx/dist
             dy=dy/dist
             dz=dz/dist	    
		   goto 999
	     end if
	     tmp1=calpha*(zeta+axdist*talpha)
c          tmp1 = |T-B|
	     if(tmp1.le.0.) then
c below the tip (zona 4)
             dist=sqrt(tmp)
             dx=dx/dist
             dy=dy/dist
             dz=dz/dist	    
		   goto 999
           elseif(tmp1.ge.modul/calpha) then
c beyond basis,  external     (remaining of zone 2)
             tmp1=radius/axdist
c here axdist > 0 for sure, dist probably too
             tmp2=(zeta*tmp1)/modul-1.
             dx=dx*(1.-tmp1)+vectz(1)*tmp2
             dy=dy*(1.-tmp1)+vectz(2)*tmp2
             dz=dz*(1.-tmp1)+vectz(3)*tmp2
             dist=sqrt(dx*dx+dy*dy+dz*dz)
             dx=dx/dist
             dy=dy/dist
             dz=dz/dist	    
		   goto 999
           elseif(axdist/zeta.gt.talpha) then
c lateral external (zone 3)
             tmp1=calpha/axdist
             tmp2=calpha*(-talpha-zeta/axdist)/modul
             dx=dx*tmp1+vectz(1)*tmp2
             dy=dy*tmp1+vectz(2)*tmp2
             dz=dz*tmp1+vectz(3)*tmp2
             dist=calpha*(axdist-zeta*talpha)
            goto 999
	     end if
c ----------------------------------------------------------------------
         else
           if (xxx.le.0.) then
             if(czeta.le.tang*(radius-axdist)) then
c exactly below the basis (zona 1)
c internal, close to the basis (zona 5)
	         dist=-czeta
	         goto 999
	       elseif(axdist.le.zeta*talpha) then
	         if (axdist.gt.0.)then
c internal, close to lateral (zona 6)
                 dist=calpha*(axdist-zeta*talpha)
		       goto 999
	         else
                 if (zeta.ge.0.) then
c exactly on the axis
                   dist=-zeta*talpha*calpha
		         goto 999
	           endif
	         end if
             end if
           elseif(czeta.lt.0.) then
c below basis, laterally (part of zona 2)
             tmp1=radius/axdist
c here axdist > 0 for sure, dist probably too
             tmp2=(zeta*tmp1)/modul-1.
             dx=dx*(1.-tmp1)+vectz(1)*tmp2
             dy=dy*(1.-tmp1)+vectz(2)*tmp2
             dz=dz*(1.-tmp1)+vectz(3)*tmp2
             dist=sqrt(dx*dx+dy*dy+dz*dz)	    
		   goto 999
	     end if
	     tmp1=calpha*(zeta+axdist*talpha)
c          tmp1 = |T-B|
	     if(tmp1.le.0.) then
c below the tip (zona 4)
             dist=sqrt(tmp)
		   goto 999
           elseif(tmp1.ge.modul/calpha) then
c beyond basis,  external     (remaining of zone 2)
             tmp1=radius/axdist
c here axdist > 0 for sure, dist probably too
             tmp2=(zeta*tmp1)/modul-1.
             dx=dx*(1.-tmp1)+vectz(1)*tmp2
             dy=dy*(1.-tmp1)+vectz(2)*tmp2
             dz=dz*(1.-tmp1)+vectz(3)*tmp2
             dist=sqrt(dx*dx+dy*dy+dz*dz)
             dx=dx/dist
             dy=dy/dist
             dz=dz/dist	    
		   goto 999
           elseif(axdist/zeta.gt.talpha) then
c lateral external (zone 3)
             dist=calpha*(axdist-zeta*talpha)
            goto 999
	     end if
         end if
c   -------------------------------------------------------------------
      else
c else di objecttype 3
         if(objecttype.eq.4) then
c box
           read(strtmp(20:80),*)xa,xb,xc,xd
           read(strtmp1,'(12f8.3)')vectz,modul,vecty,mody,vectx,
     &  modx
c normalize vectors
           do i=1,3
             vectx(i)=vectx(i)/modx
             vecty(i)=vecty(i)/mody
             vectz(i)=vectz(i)/modul
           end do
c calulate projections
           dx=xp(1)-xa(1)
           dy=xp(2)-xa(2)
           dz=xp(3)-xa(3)
c dot=(P-A)(D-A) new notation
           dot=dx*vectz(1)+dy*vectz(2)+dz*vectz(3)
           prz=dot
c dot=(P-A)(C-A) new notation
           dot=dx*vecty(1)+dy*vecty(2)+dz*vecty(3)
           pry=dot
c dot=(P-A)(B-A) new notation
           dot=dx*vectx(1)+dy*vectx(2)+dz*vectx(3)
           prx=dot
c calculate components of the projection of P onto the extended surface
c change temporarily modulus values only for convenience
           modx=modx+extent
           mody=mody+extent
           modul=modul+extent
           sign=1
           if(prx.lt.-extent) then 
             prx=-extent
             sign=0
           else
             if(prx.gt.modx) then
               prx=modx
               sign=0
             end if
           end if
           if(pry.lt.-extent) then
             pry=-extent
             sign=0
           else
             if(pry.gt.mody) then
               pry=mody
               sign=0
             end if
           end if
           if(prz.lt.-extent) then
             prz=-extent
             sign=0
           else
             if(prz.gt.modul) then
               prz=modul
               sign=0
             end if
           end if
c now if sign=1 we are inside the box so dist will eventually be negative
           if (sign.eq.1) then
             i=1
             dist=prx+extent
             if (prx+extent.gt.modx-prx) then
               i=-1
               dist=modx-prx
             end if
             if (pry+extent.lt.dist) then
               i=2
               dist=pry+extent
             end if
             if (mody-pry.lt.dist) then
               i=-2
               dist=mody-pry
             end if
             if (prz+extent.lt.dist) then
               i=3
               dist=prz+extent
             end if
             if (modul-prz.lt.dist) then
               i=-3
               dist=modul-prz
             end if
             dist=-dist
c            i moduli li lascio sbagliati perche' non li riuso
c            mody=mody-extent
c            modul=modul-extent
             if(inout) go to 90
             go to(20,30,40,50,60,70,80) i+4
50           continue
               write(6,*)'problems in distobject'
               stop
20           continue
               dx=vectz(1)
               dy=vectz(2)
               dz=vectz(3)
               go to 90
30           continue
               dx=vecty(1)
               dy=vecty(2)
               dz=vecty(3)
               go to 90
40           continue
               dx=vectx(1)
               dy=vectx(2)
               dz=vectx(3)
               go to 90
60           continue
               dx=-vectx(1)
               dy=-vectx(2)
               dz=-vectx(3)
               go to 90
70           continue
               dx=-vecty(1)
               dy=-vecty(2)
               dz=-vecty(3)
               go to 90
80           continue
               dx=-vectz(1)
               dy=-vectz(2)
               dz=-vectz(3)
90           continue
           else    
c            i moduli li lascio sbagliati perche' non li riuso
c            modx=modx-extent
c            mody=mody-extent
c            modul=modul-extent
             dx=prx*vectx(1)+pry*vecty(1)+prz*vectz(1)-dx
             dy=prx*vectx(2)+pry*vecty(2)+prz*vectz(2)-dy
             dz=prx*vectx(3)+pry*vecty(3)+prz*vectz(3)-dz
             dist=sqrt(dx*dx+dy*dy+dz*dz)
             if (.not.inout) then
               if (dist.eq.0.) write(6,*)"Wrong step in distobj!"
                 dx=-dx/dist
                 dy=-dy/dist
                 dz=-dz/dist
               end if
             end if
           end if
         end if 
       end if
      end if
888   continue
      goto 999



c   bgp lies on extended surface or center or axdist=0, correcting 
      pntsh=0.
      if (imovepoint) then
         pntsh=1.e-5
         if (vectz(2).eq.0.) tet=1.5707963268
         else tet=atan(vectz(3)/vectz(2))
         xp(2)=xp(2)+pntsh*sin(tet)
         xp(3)=xp(3)+pntsh*cos(tet)
      end if
      call distobj(xp,dx,dy,dz,nobject,ii,prbrad+1.e-5-pntsh,dist,inout
     &,zeta,axdist)
      if (.not.imovepoint) dist=0.0

999   continue

      dist=dist-prbrad

      return
      end

