      SUBROUTINE TGETGR(NAR,afac,cg,V1,V2,V3,X,AREA,GR)

c This subroutine fills in the GR matrix of integrals over the BE elements.
c The integrals over diagonal elements are perform exactly. Aragon 8/2/2001.
c The integrals over off-diagonal blocks are done in very precisely in the 
c subroutine triangleint.  See that subroutine for details.  Aragon 7/6/2001
c  Major Changes to routine:  Aragon 8/2/2001
c  Subroutine computes the InCenter of each triangle. Exac diagonal element integration.
c  Subroutine computes the normal vector of each triangle as it needs it.
c	This obviates the need for the array xn(n,3).
c  Double precision functions are called explicitly to make sure other users will
c  always compile with the -r8 compiler option. Code made more efficient in divisions.
c  Aragon 7/14/2003

      dimension g(3,3),rr(3,3)
      dimension gr(3*nar,3*nar),x(nar,3),V1(NAR,3),V2(NAR,3),V3(NAR,3)
      pi = 3.1415926535897932385
      nar3=3*nar

c  The vertices are in vv1(n,j),vv2(n,j),vv3(n,j) for the nth triangle, coords. j
c  The patch normal is in normx,normy,normz
c  The patch InCenter is in x(n,3)
c  The patch area is in Area(n)
c  The lengths of the sides of the triangle are a,b,c
c  The angles opposite to these sides are aa,bb,cc.

c	COMPUTE THE DIAGONAL BLOCKS 
	g(1,3) = 0.0
	g(2,3) = 0.0
	g(3,1) = 0.0
	g(3,2) = 0.0

      DO N1=1,NAR
      I0=3*(N1-1)
c Compute the InCenter of the triangle (this is used in triangleint).
	v21x = v2(n1,1)-v1(n1,1)
	v21y = v2(n1,2)-v1(n1,2)
	v21z = v2(n1,3)-v1(n1,3)
	v31x = v3(n1,1)-v1(n1,1)
	v31y = v3(n1,2)-v1(n1,2)
	v31z = v3(n1,3)-v1(n1,3)
	v32x = v3(n1,1)-v2(n1,1)
	v32y = v3(n1,2)-v2(n1,2)
	v32z = v3(n1,3)-v2(n1,3)
	a = dsqrt(v21x*v21x+v21y*v21y+v21z*v21z)
	b = dsqrt(v31x*v31x+v31y*v31y+v31z*v31z)
	c = dsqrt(v32x*v32x+v32y*v32y+v32z*v32z)
	cosC = (v21x*v31x+v21y*v31y+v21z*v31z)/a/b
	temp = (a+b-c)/(2.*(1.+cosC))
	x(n1,1) = v1(n1,1)+temp*(v21x/a+v31x/b)
	x(n1,2) = v1(n1,2)+temp*(v21y/a+v31y/b)
	x(n1,3) = v1(n1,3)+temp*(v21z/a+v31z/b)
c	print *, 'InCenter coordinates'
c	print *, x(n1,1),x(n1,2),x(n1,3)
c  Compute the normal and area of the triangle 
	xnorm = v31y*v21z-v31z*v21y
	ynorm = v31z*v21x-v31x*v21z
	znorm = v31x*v21y-v31y*v21x
 	area2 = Sqrt(xnorm**2 + ynorm**2 +znorm**2)/2.
c	Must compute norm with the corrected area
	xnorm = .5*xnorm/area2
	ynorm = .5*ynorm/area2
	znorm = .5*znorm/area2

c	COMPUTE THE DIAGONAL BLOCKS
c	Compute the angles aa,bb,cc, the InCircle radius
	cc = dacos(cosC)
        r = temp*dsin(cc)
c	print *, ' r = ',r
	cosB = -(v21x*v32x+v21y*v32y+v21z*v32z)/a/c
	bb = dacos(cosB)
	aa = pi -bb-cc
c	cosA = (v31x*v32x+v31y*v32y+v31z*v32z)/b/c
c	aa = acos(cosA)
c	Print *,'Pi = ',aa+bb+cc
c	Compute the incenter to vertex distances
	temp2 = 1./(2.*dcos(0.5*cc))
	dc = (a+b-c)*temp2
	db = dsqrt(a*a+(b-c)*(b-c)-2.*a*(b-c)*cosC)*temp2
	da = dsqrt(b*b+(a-c)*(a-c)-2.*b*(a-c)*cosC)*temp2
c	print *, 'Incenter to vertex distances'
c	print *, 'da = ',da,' db = ',db,' dc = ',dc

c	Values of the integrals in terms of ratio ha = r/da, etc.
	ha = r/da
	hb = r/db
	hc = r/dc
	dpa=(pi-aa)*(da-r)-2.0*(da*dacos(ha) -
     $  r*dlog((1.+dsqrt(1.-ha*ha))/ha))
	dpb=(pi-bb)*(db-r)-2.0*(db*dacos(hb) -
     $  r*dlog((1.+dsqrt(1.-hb*hb))/hb))
	dpc=(pi-cc)*(dc-r)-2.0*(dc*dacos(hc) -
     $  r*dlog((1.+dsqrt(1.-hc*hc))/hc))
c	print *,'dpa= ',dpa,' dpb= ',dpb,' dpc= ',dpc
	sca= 0.5*r*((dcos(2.0*cc)+dcos(2.0*bb))*si(ha) -
     $              (dsin(2.0*cc)+dsin(2.0*bb))*csi(ha))
	scb=0.5*r*((1.+dcos(2.0*bb))*si(hb)+dsin(2.0*bb)*csi(hb))
	scc=0.5*r*((1.+dcos(2.0*cc))*si(hc)+dsin(2.0*cc)*csi(hc))
c	print *,'sca= ',sca,' scb= ',scb,' scc= ',scc

	cca= -dsin(bb-cc)*r*(dcos(aa)*si(ha)+dsin(aa)*csi(ha))
	ccb= dsin(bb)*r*(dcos(bb)*si(hb)+dsin(bb)*csi(hb))
	ccc= -dsin(cc)*r*(dcos(cc)*si(hc)+dsin(cc)*csi(hc))
c	print *,'cca= ',cca,' ccb= ',ccb,' ccc= ',ccc

	cor = dsqrt(afac)
        g(1,1) = cor*((2.0+cg)*Pi*r+(1.0+0.5*cg)*cg*(dpa+dpb+dpc)
     $+cg*0.25*(sca+scb+scc))
        g(1,2) = 0.5*cor*(cca+ccb+ccc)*cg
        g(2,1) = g(1,2)
        g(2,2) = g(1,1) - 0.5*cor*cg*(sca+scb+scc)
        g(3,3) = cor*(2.0*Pi*r+(dpa+dpb+dpc))

c	print *,' g tensor '
c	print *, g(1,1),g(1,2),g(1,3)
c	print *, g(2,1),g(2,2),g(2,3)
c	print *, g(3,1),g(3,2),g(3,3)
c	pause
c	Rotate the G matrix to the lab frame. Compute the
c	direction cosine matrix = rotation matrix.
	rr(1,1) = v21x/a
	rr(1,2) = v21y/a
	rr(1,3) = v21z/a
	rr(3,1) = xnorm
	rr(3,2) = ynorm
	rr(3,3) = znorm
	rr(2,1) = (ynorm*v21z-znorm*v21y)/a
	rr(2,2) = (znorm*v21x-xnorm*v21z)/a
	rr(2,3) = (xnorm*v21y-ynorm*v21x)/a

c	print *,' rotation tensor '
c	print *, rr(1,1),rr(1,2),rr(1,3)
c	print *, rr(2,1),rr(2,2),rr(2,3)
c	print *, rr(3,1),rr(3,2),rr(3,3)
c	pause
	do i=1,3
	do j=1,3
	temp = 0.0
	do k=1,3
	temp = temp+rr(k,i)*(g(k,1)*rr(1,j)+g(k,2)*rr(2,j)+
     $  g(k,3)*rr(3,j))
	enddo
	gr(I0+i,I0+j) = temp
	enddo
	enddo
	
	enddo
c	COMPUTE THE OFF-DIAGONAL BLOCKS
c	The n1 index selects the external triangle, the n2 index
c	is the triangle we integrate over.  We need the InCenter
c	of the external triangle, and we need the vertices of the
c	triangle being integrated over. These are formed inside triangleint.
c	The pairs are independent.  When compiled with OpenMP, the rows N1 are
c	shared out among threads; each thread has its own copy of COMMON /dots/.

!$OMP PARALLEL DO PRIVATE(N2,I0,J0,J,G) SCHEDULE(DYNAMIC,8)
       DO N1=1,NAR
       DO N2=1,NAR
       IF (N1.NE.N2) THEN
       I0=3*(N1-1)
       J0=3*(N2-1)
       CALL triangleint(G,n2,n1,cg,V1,V2,V3,X,AREA,NAR)
c Note:  the area correction is done in the Main program by rescaling each triangle
c  area.  Thus, when triangleint multiplies the integrals by the area, it is done.
       DO J=1,3    
       GR(I0+1,J0+J)=G(1,J)
       GR(I0+2,J0+J)=G(2,J)
       GR(I0+3,J0+J)=G(3,J)
       ENDDO
       ENDIF
       ENDDO
       ENDDO
!$OMP END PARALLEL DO
       RETURN
       END

	double precision function si(x)
	real*8 x
	si = 2.0*(dlog((1.0+dsqrt(1.0-x*x))/x)-dsqrt(1.0-x*x))
	return
	end

	double precision function csi(x)
	real*8 x
	csi = 3.0-1.0/x-2.0*x
	return
	end

