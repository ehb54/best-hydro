      SUBROUTINE triangleint(G,n2,n1,cg,V1,V2,V3,X,AREA,NAR)
c This subroutine integrates the Oseen Tensor over a reference
c  triangle which is a distance away from another triangle.
c  The first of the double integrals is done analytically. The
c  second integral is done via 21st order Gaussian quadrature.
c  This is necessary to compute the integrals to double precision
c  accuracy.
c  Aragon.  7/10/2003

       common /dots/ d212c,d2132,d2121,d2c2c,d322c,d3232,d2c1c,
     $ d211c,d2131,d1c1c,d1c31,d3131,d3c31,d3c3c,d3231,d3c32
       dimension v21(3),v31(3),v32(3),v1c(3),v2c(3),v3c(3)
       DIMENSION G(3,3),w(21),y(21),x(nar,3),V1(NAR,3),V2(NAR,3),
     $ V3(NAR,3),AREA(NAR)

       DO I=1,3
       G(1,i)=0.0
       G(2,i)=0.0
       G(3,i)=0.0
       enddo

c  select the vertices of the triangle being integrated 
c  Form all the subtractions the integrations will need.

	v21(1)= v2(n2,1)-v1(n2,1)
	v21(2)= v2(n2,2)-v1(n2,2)
	v21(3)= v2(n2,3)-v1(n2,3)
	v31(1)= v3(n2,1)-v1(n2,1)
	v31(2)= v3(n2,2)-v1(n2,2)
	v31(3)= v3(n2,3)-v1(n2,3)
	v32(1)= v3(n2,1)-v2(n2,1)
	v32(2)= v3(n2,2)-v2(n2,2)
	v32(3)= v3(n2,3)-v2(n2,3)
	v1c(1)= v1(n2,1)-x(n1,1)
	v1c(2)= v1(n2,2)-x(n1,2)
	v1c(3)= v1(n2,3)-x(n1,3)
	v2c(1)= v2(n2,1)-x(n1,1)
	v2c(2)= v2(n2,2)-x(n1,2)
	v2c(3)= v2(n2,3)-x(n1,3)
	v3c(1)= v3(n2,1)-x(n1,1)
	v3c(2)= v3(n2,2)-x(n1,2)
	v3c(3)= v3(n2,3)-x(n1,3)

c  Form the dot products the integrations will need.

	d2121 = dot(v21,v21)	
	d2131 = dot(v21,v31)
	d2132 = dot(v21,v32)	
	d211c = dot(v21,v1c)
	d212c = dot(v21,v2c)
	d3131 = dot(v31,v31)
	d3232 = dot(v32,v32)
	d3231 = dot(v32,v31)
	d322c = dot(v32,v2c)
	d1c1c = dot(v1c,v1c)
	d1c31 = dot(v1c,v31)	
	d2c2c = dot(v2c,v2c)
	d2c1c = dot(v2c,v1c)		
	d3c31 = dot(v3c,v31)
	d3c32 = dot(v3c,v32)
	d3c3c = dot(v3c,v3c)	

c  These are the abscissas and weights for 21 st. Order Gaussian Q.

	y(1)= 3.1239146898052499D-3
	y(2)= 1.6386580716846853D-2
	y(3)= 3.9950332924799586D-2
	y(4)= 7.3318317708341358D-2
	y(5)= 0.1157800182621610
	y(6)= 0.1664305979012938
	y(7)= 0.2241905820563901
	y(8)= 0.2878289398962806
	y(9)= 0.3559893415987995
	y(10)=0.4272190729195525
	y(11)=0.5
	y(12)=0.5727809270804476
	y(13)=0.6440106584012005
	y(14)=0.7121710601037194
	y(15)=0.7758094179436099
	y(16)=0.8335694020987062
	y(17)=0.8842199817378390 
	y(18)=0.9266816822916586
	y(19)=0.9600496670752004
	y(20)=0.9836134192831531
	y(21)=0.9968760853101948

	w(1)= 1.601722825777433D-2
	w(2)= 3.695378977085249D-2
	w(3)= 5.713442542685720D-2
	w(4)= 7.610011362837930D-2
	w(5)= 9.344442345603390D-2
	w(6)= 0.1087972991671484
	w(7)= 0.1218314160537285
	w(8)= 0.1322689386333375
	w(9)= 0.1398873947910731
	w(10)=0.1445244039899701
	w(11)=0.1460811336496904
	w(12)=w(10)
	w(13)=w(9)
	w(14)=w(8)
	w(15)=w(7)
	w(16)=w(6)
	w(17)=w(5)
	w(18)=w(4)
	w(19)=w(3)
	w(20)=w(2)
	w(21)=w(1)


c The scalar part of the integral and the tensor part of the integrals.
c  First compute the seven integrals

	t0 = 0.0
	s0 = 0.0
	sx = 0.0
	sy = 0.0
	sxy = 0.0
	syy = 0.0
	sxx = 0.0

	do j=1,21
	t0 = t0+w(j)*h(y(j))
	s0 = s0+w(j)*hy(y(j))
	sx = sx+w(j)*hxy(y(j))
	sy = sy+w(j)*y(j)*hy(y(j))
	sxy = sxy+w(j)*y(j)*hxy(y(j))
	syy = syy+w(j)*y(j)*y(j)*hy(y(j))
	sxx = sxx+w(j)*y(j)*y(j)*hx(y(j))
	enddo

	t0 = 0.5*t0
	s0 = 0.5*s0
	sx = 0.5*sx
	sy = 0.5*sy
	sxy = 0.5*sxy
	syy = 0.5*syy
	sxx = 0.5*sxx

c	print *, 't0= ',t0
c	print *, 's0= ',s0
c	print *, 'sx= ',sx
c	print *, 'sy= ',sy
c	print *, 'sxy= ',sxy
c	print *, 'sxx= ',sxx
c	print *, 'syy= ',syy

c  Now compute the tensor components of Gij    	
	t11 = v1c(1)*v1c(1)*s0+(v1c(1)*v21(1)+v21(1)*v1c(1))*sx
     $ +(v31(1)*v1c(1)+v1c(1)*v31(1))*sy+
     $ (v21(1)*v31(1)+v31(1)*v21(1))*sxy+v21(1)*v21(1)*sxx+
     $ v31(1)*v31(1)*syy

	t12 = v1c(1)*v1c(2)*s0+(v1c(1)*v21(2)+v21(1)*v1c(2))*sx
     $ +(v31(1)*v1c(2)+v1c(1)*v31(2))*sy+
     $ (v21(1)*v31(2)+v31(1)*v21(2))*sxy+v21(1)*v21(2)*sxx+
     $ v31(1)*v31(2)*syy

	t13 = v1c(1)*v1c(3)*s0+(v1c(1)*v21(3)+v21(1)*v1c(3))*sx
     $ +(v31(1)*v1c(3)+v1c(1)*v31(3))*sy+
     $ (v21(1)*v31(3)+v31(1)*v21(3))*sxy+v21(1)*v21(3)*sxx+
     $ v31(1)*v31(3)*syy

	t22 = v1c(2)*v1c(2)*s0+(v1c(2)*v21(2)+v21(2)*v1c(2))*sx
     $ +(v31(2)*v1c(2)+v1c(2)*v31(2))*sy+
     $ (v21(2)*v31(2)+v31(2)*v21(2))*sxy+v21(2)*v21(2)*sxx+
     $ v31(2)*v31(2)*syy

	t23 = v1c(2)*v1c(3)*s0+(v1c(2)*v21(3)+v21(2)*v1c(3))*sx
     $ +(v31(2)*v1c(3)+v1c(2)*v31(3))*sy+
     $ (v21(2)*v31(3)+v31(2)*v21(3))*sxy+v21(2)*v21(3)*sxx+
     $ v31(2)*v31(3)*syy

	t33 = v1c(3)*v1c(3)*s0+(v1c(3)*v21(3)+v21(3)*v1c(3))*sx
     $ +(v31(3)*v1c(3)+v1c(3)*v31(3))*sy+
     $ (v21(3)*v31(3)+v31(3)*v21(3))*sxy+v21(3)*v21(3)*sxx+
     $ v31(3)*v31(3)*syy

      
      G(1,1)= 2.0*area(n2)*(t0+cg*t11)
      G(2,1)= 2.0*area(n2)*t12*cg
      G(2,2)= 2.0*area(n2)*(t0+cg*t22)
      G(3,1)= 2.0*area(n2)*t13*cg
      G(3,2)= 2.0*area(n2)*t23*cg
      G(3,3)= 2.0*area(n2)*(t0+cg*t33)
      G(1,2)=G(2,1)
      G(1,3)=G(3,1)
      G(2,3)=G(3,2)
      RETURN
      END

	double precision function dot(a,b)
	dimension a(3),b(3)
	dot = a(1)*b(1)+a(2)*b(2)+a(3)*b(3)
	return
	end

	double precision function h(y)
	common /dots/ d212c,d2132,d2121,d2c2c,d322c,d3232,d2c1c,
     $ d211c,d2131,d1c1c,d1c31,d3131,d3c31,d3c3c,d3231,d3c32
	h = DLog((D212c+D2132*y+DSqrt(D2121)*DSqrt(D2c2c+ 
     $ 2.0*y*D322c+y*y*D3232))/(D211c+D2131*y+DSqrt(D2121)*
     $ DSqrt(D1c1c+2*D1c31*y+D3131*y*y)))/DSqrt(D2121)
	return
	end

	double precision function hx(x)
	common /dots/ d212c,d2132,d2121,d2c2c,d322c,d3232,d2c1c,
     $ d211c,d2131,d1c1c,d1c31,d3131,d3c31,d3c3c,d3231,d3c32
	 a= (-D1c31 - x*D2131)/
     $ DSqrt(D1c1c+2.0*x*D211c+x*x*D2121)
	 b = (-x*D3231+D3c31)/
     $ DSqrt(x*x*D3232 - 2.0*x*D3c32 + D3c3c)
	 c = (-D1c31*D1c31 - 2.0*x*D1c31*D2131 - x*x*D2131*D2131 + 
     $ D1c1c*D3131 + 2.0*x*D211c*D3131 + x*x*D2121*D3131)
	hx = (a+b)/c 
	return
	end

	double precision function hy(y)
	common /dots/ d212c,d2132,d2121,d2c2c,d322c,d3232,d2c1c,
     $ d211c,d2131,d1c1c,d1c31,d3131,d3c31,d3c3c,d3231,d3c32
	a= (D211c + y*D2131)/
     $  DSqrt(D1c1c+2.0*y*D1c31+y*y*D3131)
	b= (y*D2132 + D212c)/
     $     DSqrt(D2c2c+2.0*y*D322c+y*y*D3232)
	c= (D211c*D211c - D1c1c*D2121 - 2.0*y*D1c31*D2121+ 
     $    2.0*y*D211c*D2131 + y*y*D2131*D2131 - y*y*D2121*D3131)
	hy = (a -b )/c
	return
	end

	double precision function hxy(y)
	common /dots/ d212c,d2132,d2121,d2c2c,d322c,d3232,d2c1c,
     $ d211c,d2131,d1c1c,d1c31,d3131,d3c31,d3c3c,d3231,d3c32
	a= DSqrt(D1c1c + 2.0*y*D1c31 + y*y*D3131)
	b= (y*D211c - D2c1c - 2.0*y*D1c31 - y*D2131 - 
     $   y*y*D3231)/DSqrt(D2c2c + 2.0*y*D322c + y*y*D3232)
	c= (-D211c*D211c + D1c1c*D2121 + 2.0*y*D1c31*D2121 - 
     $  2.0*y*D211c*D2131 - y*y*(D2131*D2131) + y*y*D2121*D3131)
	hxy = (a+b)/c  
	return
	end
