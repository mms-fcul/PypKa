      subroutine phintp(gp,phi)
c
c interpolates the potential at any point inside
c a cubical volume using the potential values at the
c 8 vertices by means of a trilinear function:
C	
C W = A1.XYZ + A2.XY + A3.XZ + A4.YZ + A5.X +
C		A6.Y + A7.Z + A8
C	
C where Ai coefficients are linear combinations of Wi at
C the cube corners
c
	include 'qdiffpar4.h'
	include 'qlog.h'
c---------------------------------------------------

      dimension gp(3),phimap(igrid,igrid,igrid)
c
c return 0.0 if outside grid
c
	rgrid = float(igrid)
	do 9000 k = 1,3
	  if((gp(k).lt.1).or.(gp(k).gt.rgrid)) then
	    phi = 0.0
            write(6,*)"Pay attention, point out of the cube!!"
            write(6,*)"i=",k,"Value:",gp(k),"Igrid:",igrid
	    return
	  end if
c	end do
9000	continue
c
c find lower left bottom grid point
c
      nx = int(gp(1))
	nx1 = nx+1
	if(nx1.gt.igrid)nx1=nx
      ny = int(gp(2))
	ny1 = ny+1
	if(ny1.gt.igrid)ny1=ny
      nz = int(gp(3))
	nz1 = nz+1
	if(nz1.gt.igrid)nz1=nz
c
c calculate cube coordinates of point
c
	xg = nx
        yg = ny
	zg = nz
	xgr = gp(1) - xg
	ygr = gp(2) - yg
	zgr = gp(3) - zg
c
c calculate coefficients of trilinear function
c
	a8 = phimap(nx,ny,nz)
	a7 = phimap(nx,ny,nz1) - a8
	a6 = phimap(nx,ny1,nz) - a8
	a5 = phimap(nx1,ny,nz) - a8
	a4 = phimap(nx,ny1,nz1) - a8 - a7 - a6
	a3 = phimap(nx1,ny,nz1) - a8 - a7 - a5
	a2 = phimap(nx1,ny1,nz) - a8 - a6 - a5
	a1 = phimap(nx1,ny1,nz1) - a8 - a7 - a6
     &   - a5 - a4 - a3 - a2
c
c determine value of phi
c
	phi = a1*xgr*ygr*zgr + a2*xgr*ygr + a3*xgr*zgr
     &     + a4*ygr*zgr + a5*xgr + a6*ygr + a7*zgr + a8

      return
      end

      subroutine readgreen

	real green(-10:10,-10:10,-10:10)

	common /gree/green


      open(52,file='Green10.bin',form='unformatted')
	   read(52)green
      close(52)

      return
      end

      subroutine phintpgreen(gp,phi,phiclean)
c
c interpolates the potential at any point inside
c a cubical volume using the potential values at the
c 8 vertices by means of a trilinear function:
C	
C W = A1.XYZ + A2.XY + A3.XZ + A4.YZ + A5.X +
C		A6.Y + A7.Z + A8
C	
C where Ai coefficients are linear combinations of Wi at
C the cube corners
c
	include 'qdiffpar4.h'
	include 'qlog.h'
c---------------------------------------------------

      dimension gp(3),phimap(igrid,igrid,igrid)
	real phi,phigrid,phiclean,green(-10:10,-10:10,-10:10)

	common /gree/green
c
      q=1.0
	eps=2.0
c
	rgrid = float(igrid)
	do 9000 k = 1,3
	  if((gp(k).lt.1).or.(gp(k).gt.rgrid)) then
	    phi = 0.0
            write(6,*)"Pay attention, point out of the cube!!"
            write(6,*)"i=",k,"Value:",gp(k),"Igrid:",igrid
	    return
	  end if
c	end do
9000	continue
c
c find lower left bottom grid point
c
      nx = int(gp(1))
	nx1 = nx+1
	if(nx1.gt.igrid)nx1=nx
      ny = int(gp(2))
	ny1 = ny+1
	if(ny1.gt.igrid)ny1=ny
      nz = int(gp(3))
	nz1 = nz+1
	if(nz1.gt.igrid)nz1=nz
c
c calculate cube coordinates of point
c
	xg = nx
      yg = ny
	zg = nz
	xgr = gp(1) - xg
	ygr = gp(2) - yg
	zgr = gp(3) - zg
c
c calculate coefficients of trilinear function
c
	a8 = phimap(nx,ny,nz)
	a7 = phimap(nx,ny,nz1) - a8
	a6 = phimap(nx,ny1,nz) - a8
	a5 = phimap(nx1,ny,nz) - a8
	a4 = phimap(nx,ny1,nz1) - a8 - a7 - a6
	a3 = phimap(nx1,ny,nz1) - a8 - a7 - a5
	a2 = phimap(nx1,ny1,nz) - a8 - a6 - a5
	a1 = phimap(nx1,ny1,nz1) - a8 - a7 - a6
     &   - a5 - a4 - a3 - a2
c
c determine value of phi
c
	phi = a1*xgr*ygr*zgr + a2*xgr*ygr + a3*xgr*zgr
     &     + a4*ygr*zgr + a5*xgr + a6*ygr + a7*zgr + a8

c determine grid position with respect to grid center
c     pay attention! is valid only if acenter is null!!!
      icx=(igrid+1)/2
      icy=icx
      icz=icx
	inx=nx-icx
      iny=ny-icy
	inz=nz-icz

	inx1 = inx+1
	iny1 = iny+1
	inz1 = inz+1
	if(iabs(inx1).ge.10.or.iabs(iny1).ge.10.or.iabs(inz1).ge.10)
     &write(6,*)"we are far!!!!"

c calculate coefficients of trilinear function
c
	a8 = green(inx,iny,inz)
	a7 = green(inx,iny,inz1) - a8
	a6 = green(inx,iny1,inz) - a8
	a5 = green(inx1,iny,inz) - a8
	a4 = green(inx,iny1,inz1) - a8 - a7 - a6
	a3 = green(inx1,iny,inz1) - a8 - a7 - a5
	a2 = green(inx1,iny1,inz) - a8 - a6 - a5
	a1 = green(inx1,iny1,inz1) - a8 - a7 - a6
     &   - a5 - a4 - a3 - a2
c     14099...=561*4*pi*2
c     attention: somewhere 561 has been modified and elsewhere it has not
	phigreen =14099.46782931099*q*(a1*xgr*ygr*zgr + a2*xgr*ygr + 
     &a3*xgr*zgr+ a4*ygr*zgr + a5*xgr + a6*ygr + a7*zgr + a8)/eps

c     calculate phiclean
      r=sqrt((gp(1)-icx)**2+(gp(2)-icy)**2+(gp(3)-icz)**2)/scale
      phiclean= phi+561.9969367787*q/(eps*r)-phigreen
      write(6,*)"phi:",phi,561.9969367787*q/(eps*r),phigreen
      write(6,*)"G0:",green(0,0,0)
      return
      end

      subroutine debtp(gp,deb)
c
c interpolates the debye map at any point inside
c a cubical volume using idebmap values at the
c 8 vertices by means of a trilinear function:
C	
C W = A1.XYZ + A2.XY + A3.XZ + A4.YZ + A5.X +
C		A6.Y + A7.Z + A8
C	
C where Ai coefficients are linear combinations of Wi at
C the cube corners
c
	include 'qdiffpar4.h'
	include 'qlog.h'
c---------------------------------------------------

      dimension gp(3)
	logical*1 idebmap(igrid,igrid,igrid)
	real debm(8)

c
c return 0.0 if outside grid
c
	rgrid = float(igrid)
	do 9000 k = 1,3
	  if((gp(k).lt.1).or.(gp(k).gt.rgrid)) then
	    deb = 0.0
            write(6,*)"Pay attention, point out of the cube!!"
            write(6,*)"i=",k,"Value:",gp(k),"Igrid:",igrid
	    return
	  end if
c	end do
9000	continue
c
c find lower left bottom grid point
c
      nx = int(gp(1))
	nx1 = nx+1
	if(nx1.gt.igrid)nx1=nx
      ny = int(gp(2))
	ny1 = ny+1
	if(ny1.gt.igrid)ny1=ny
      nz = int(gp(3))
	nz1 = nz+1
	if(nz1.gt.igrid)nz1=nz
c
c calculate cube coordinates of point
c
	xg = nx
      yg = ny
	zg = nz
	xgr = gp(1) - xg
	ygr = gp(2) - yg
	zgr = gp(3) - zg

c convert logical idebmap to real
	do ii = 1,8
	   debm(ii)=0.0
	end do
      if (idebmap(nx,ny,nz)) debm(1)=1.0
	if (idebmap(nx,ny,nz1)) debm(2)=1.0
	if (idebmap(nx,ny1,nz)) debm(3)=1.0
	if (idebmap(nx1,ny,nz)) debm(4)=1.0
	if (idebmap(nx,ny1,nz1)) debm(5)=1.0
	if (idebmap(nx1,ny,nz1)) debm(6)=1.0
	if (idebmap(nx1,ny1,nz)) debm(7)=1.0
	if (idebmap(nx1,ny1,nz1)) debm(8)=1.0

c
c calculate coefficients of trilinear function
c
	a8 = debm(1)
	a7 = debm(2) - a8
	a6 = debm(3) - a8
	a5 = debm(4) - a8
	a4 = debm(5) - a8 - a7 - a6
	a3 = debm(6) - a8 - a7 - a5
	a2 = debm(7) - a8 - a6 - a5
	a1 = debm(8) - a8 - a7 - a6
     &   - a5 - a4 - a3 - a2
c
c determine value of phi
c
	deb = a1*xgr*ygr*zgr + a2*xgr*ygr + a3*xgr*zgr
     &     + a4*ygr*zgr + a5*xgr + a6*ygr + a7*zgr + a8

      return
      end

	
