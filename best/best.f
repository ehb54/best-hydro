       PROGRAM BEST90

c      AUTHOR: Sergio Aragon
c			Professor of Chemistry
c			San Francisco State University
c			aragons@sfsu.edu

c      DESCRIPTION

c      BEST performs hydrodynamic Boundary Element (BE) calculations with triangular 
c      shaped boundary elements. Stick boundary conditions are assumed. The computation
c	 includes translational and rotational diffusion tensors and their coupling, the 
c	 viscosity factor, and the intrinsic viscosity.

c      METHOD
c      BEST implements a precise numerical approximation to the exact solution of the 
c      Stokes Equations as developed by Youngren and Acrivos J. Fluid Mech. 69,377 (1975).
c      The program performs a RESISTANCE calculation.  The friction tensors so obtained
c      are later inverted to obtain the DIFFUSION tensors.  The computation moves a body
c      in an otherwise quiescent unbounded fluid and computes the total Force and Torque
c      on the body for the assumed motion. Six BE calculations are done for three translations
c      and three rotations about the x,y,z axis respectively. From these 6 calculations,
c      the complete 6x6 translation/rotation friction tensor is computed.  The diagonal blocks
c      are the translational friction tensor Ktt and the rotational friction tensor Krr.
c      The off diagonal blocks are the rotation translation coupling friction tensor Krt
c      which is the transpose of Ktr. The asymmetry in Krt is used to translate Krr to the
c      Center of Resistance. Ktt is independent of the origin. At the Center of Resistance,
c      the friction tensors are inverted to compute the diffusion tensors.  The asymmetry
c      in the translation/rotation coupling diffusion tensor Drt is used to translate the
c      translational diffusion tensor Dtt to the Center of Diffusion.
c      The main reference to this work is S. Aragon, J. Comp. Chem. 25, 1191-1205 (2004).
c	 The viscosity calculation is done by default in the centroid. The reference for this
c	 implementation is D. Hahn & S. Aragon, J. Chem. Thoery Comput. 2, 1416-1428 (2006). 
c	A calculation in the center of viscosity can be requested by the -vc flag. 
c	It is more time consuming and there is little difference from the default in most cases.  

c      INPUT
c      The program requires a triangulation of the surface of the body as input.  The main
c      input file contains the number of triangles, followed by the vertex coordinates,
c      one line per triangle. See subroutine Input for the format, or the sample file provided 
c      in the distribution. In addition, the program implements the "area" correction by default.
c      Thus, a Connolly MSROLL output file *.c3v is assumed present.  The option -na is provided in
c      case the user wishes to perform an approximate calculation without the area correction. In
c      this case, the Connolly log file does not have to be present.  See the usage statement below.

c      OUTPUT
c      The program generates an output file called [filename].be, where "filename" is the name
c      of the input vertex file.  In addition the program writes a [filename].log file with
c      supplementary information produced during the run. In version 3.90, the flags are appended
c	to the output file name in order to keep a record of calculation type.  The appended names
c	are "n", for -na, "p" for -p, "v" for -v, "vc" for -vc, and "m" for -mv.  The flags are described
c	 below.

c      USAGE
c      Typing "best" without arguments produces the usage statement:
c      'Usage : best -f input_filename [-na] [-p] [-v] [-vc] [-mw]'
c	
c      The -na option causes the program to perform a calculation wihtout the area correction.
c		In this case, the Connolly .c3v file containing the precise surface area value is not 
c               needed.
c      The -p  option causes the program to calculate diffusion coefficients and intrinsic viscosity
c		with the pre-averaged hydrodynamic interaction. This option does not save time! 
c		It is provided for esoteric reasons.
c      The -v  option causes the program to calculate the intrinsic viscosity at the centroid
c      The -vc option causes the program to calculate the intrinsic viscosity with the origin at
c               the center of viscocity.  This option is more time consuming but is needed for
c               asymmetric shapes such as bent rods.
c      The -mw option causes the program to read the molecular weight in grams per mole
c               from the input line.

c      ACKNOWLEDGEMENTS
c      I thank Stuart Allison, Georgia State University, for the initial push in this work
c      and Tilman Rosales, M.S. SFSU (2002) for debugging and some nice features he added.
c	 The viscosity calculation and many enhancements have been programmed by David Hahn. 

C      HISTORY
C      1998,Sabbatical in Spain: Modification of be program to output full translational and rotational 
c      diffusion tensors.  Program requests and input filename  INFNM and generates the output filename 
c      OUTFNM by appending ".be" to INFNM.
c      6/26/01 S. Aragon.  Inserted centroid shift after coordinate read.
c      7/3/01 Added center of diffusion correction.
c      7/7/01  Added new Tgetgr and Triangleint subroutines to work with
c      triangles. Triangleint very precisely integrates the Oseen
c      tensor over two distinct triangles.  Tgetgr approximates
c      the integrals over the same triangle by an exact calculation
c      over a circle.  This will be improved upon later.
c      8/2/01 Aragon. Introduced an exact integration of the Oseen tensor over a
c      triangle, thus eliminating the circle approximation altogether.
c      The triangle centers are now the InCenter.  Modifications occur in
c      best.f (InCenter); tgetgr.f (exact diagonal element integration.)
c      11/27/01 Aragon. Version 3.0:Eliminated the Matrix inversion in favor of the direct
c      solution of the system of equation using the LAPACK computational routine
c      SGETRS and SGETRF. This requires the installation of the INTEL Math Kernel Library
c      or equivalent to access the BLAS routines. This reduces the storage
c      requirements by almost a factor of two.
c      6/6/02 Transformed the friction tensors to Center of Resistance. Computed Diffusion from
c      these tensors. Then transformed to Center of Diffusion. Output all tensors.
c      6/15/2003  D. Hahn and S. Aragon, Computation of the viscosity factor (-v flag), the intrinsic
c      viscosity (requires reading in the molecular weight -mw flag), and the computation of
c      the viscosity in the center of viscosity (-vc flag).
C      9/13/04 D. HAHN & S. Aragon: MODULAR, F90 VERSION CREATED (VERSION 3.90) (SUBROUTINES DIFCO,
C      VISCO, & OUTPUT); Flags appended to name of output file.
C	4/14/2025 Corrected difco.f so Dtr and Ktr matrices are symmetric.  Error noticed by Jerome
C	Fung in 2020.  This is new verson 4.0 distribution.
C	11/13/2026 E. H. Brookes: The -vc calculation no longer inverts GR.  DIFCO and VISCO
C	solve for all their flows at once with the LU factors of GR (DGETRS).  This is much
C	faster, and removes the rounding error of the explicit inverse from the -vc results.
C	11/13/2026 E. H. Brookes: TGETGR fills GR in parallel when compiled with OpenMP.
C	11/13/2026 E. H. Brookes: GR is deflated before it is factored.  On a closed surface
C	the single layer of the normal field vanishes (a uniform pressure makes no flow), so
C	the unit normals are a near-null vector of GR.  The area correction had been keeping
C	GR away from singular.  Without it (-na, or a surface whose own area is exact), the
C	solutions were dominated by that pressure mode, and the torques, the centers, D_r,
C	the viscosity and even D_t were wrong: CHI 1054 instead of 3.24 for a lysozyme
C	surface of 7948 triangles, and with -na [eta] 3.21 instead of 2.04 for 2eu7_02562.
C	GR + GDEF*XMU*XMU^T, with XMU the area times the unit normal of each triangle, has no
C	such mode.  With the area correction, D_t, CHI and [eta] change by 2.3e-6 or less and
C	D_r by 1.1e-5 or less.
C	11/13/2026 E. H. Brookes: DIFCO sums all three terms of Dtr = -Kttinv Krt^T Drr.  Only
C	the last term was kept, so Dtr, the center of diffusion and Dtt were wrong.
C	11/13/2026 E. H. Brookes: DIFCO's transform of Krt and Krr to the center of resistance
C	used Krt in place of its transpose, so Krt there was not symmetric.  Drr, Dtr, the
C	center of diffusion and Dtt now agree with the direct inverse of the 6x6 friction tensor.
	

c      ARRAY DIMENSIONS
c      The dimension of the arrays related to the number of BE elements
c      is allocated at run time from the input file. 

c      CAUTION: Results are inaccurate in single precision on 32 bit processors for 1000 or
c      more triangles.  With double precision in a 32 bit processor, one can safely use up to
c      5000 triangles.  For larger numbers of triangles, double precision in a 64 bit processor
c      is required. This should not be an issue except for very large macromolecules, exceeding
c      200 kDaltons.  BEST should be used in conjunction with COALESCE to produce a series of
c      calculations as a function of the number of triangles, and extrapolated vs 1/N for precise
c      results.  The area correction should be used for such precise work (that is, provide
c      the Connolly precise surface area.)

       INTEGER,DIMENSION(:),ALLOCATABLE :: IPIV
       REAL,DIMENSION(:),ALLOCATABLE :: AREA, XMU
       REAL,DIMENSION(:,:),ALLOCATABLE :: TQ, F, W, V1, V2, V3, X, GR,
     $ B, V0
       REAL,DIMENSION(:,:,:),ALLOCATABLE :: G
       REAL,DIMENSION(:,:,:,:),ALLOCATABLE :: H
       DIMENSION FTOTH(3),TTOTH(3),OMGA(3),D(3),V(3,3),U(3),DTR(3,3),
     $ DRR(3,3),DTT(3,3),UU(3,3),AUX(3,3),BB(3,3),BBI(3,3),TEMP01(5),
     $ TEMP02(5),TEMP03(5),TEMP11(5),TEMP12(5),TEMP13(5),TEMP21(5),
     $ TEMP22(5),TEMP23(5),TEMP31(5),TEMP32(5),TEMP33(5),GG(3,3),
     $ RR(3,3),V21(3),V31(3),V32(3),V1C(3),V2C(3),V3C(3),W1(21),Y(21),
     $ URC(3),UDR(3)
       CHARACTER*80 INLINE,outfnm,outfnm2,INFNM,INFNM3,ARG
       character*16 datetime,datetime2,datetime3,ver
       real Ktt(3,3),Ktr(3,3),Krt(3,3),Krr(3,3),Kttinv(3,3),tarray2,
     $ tarray3,tarray4

c      Change the version number here.

       ver = '4.0 SFSU/ARAGON'

C      I. INPUT SURFACE TRIANGULATION/CENTROID SHIFT
	
       call input(infnm,iarea,truearea,TRUEVOL,cg,
     $ IVISCO,ICOV,MW,WM)

c      read number of surface patches

       OPEN(UNIT=10,FILE=INFNM)      

       IF(IAREA.EQ.1) THEN                                                                         
       READ(10,*)NAR,TRUEAREA,TRUEVOL
       ELSE
       READ(10,*)NAR
       ENDIF
                
       NAR3 = 3*NAR
                                                                
       ALLOCATE(V1(NAR,3),V2(NAR,3),V3(NAR,3),X(NAR,3),
     $ GR(NAR3,NAR3),AREA(NAR),V0(NAR3,8),F(NAR,3),B(NAR3,6),
     $ TQ(NAR,3),W(5,NAR3),G(NAR,5,3),H(NAR,5,3,3),IPIV(NAR3),
     $ XMU(NAR3))
                                                                                
c      read input file:
c      We read-in only the triangle vertices and the number of triangles, NAR.
c      read vertex coordinates of all triangles
                                                                                
       DO N=1,NAR
       READ(10,*)v1(N,1),v1(N,2),v1(N,3),v2(N,1),v2(N,2),v2(N,3)
     $ ,v3(N,1),v3(N,2),v3(N,3)
       enddo
       close(unit=10)

C      OPEN OUTPUT FILE 'infnm.be'
                                                                                
       do i=1,80
       l=i
       if (infnm(l:l).eq.' ') goto 150
       enddo
150    IF(IVISCO.EQ.1) THEN
       INFNM(1:L) = INFNM(1:L-1)//'v'
       L = L + 1
       ENDIF 
       IF(ICOV.EQ.1) THEN
       INFNM(1:L) = INFNM(1:L-1)//'c'
       L = L + 1
       ENDIF
       IF(MW.EQ.1) THEN
       INFNM(1:L) = INFNM(1:L-1)//'m' 
       L = L + 1
       ENDIF 
       IF(IAREA.EQ.0) THEN
       INFNM(1:L) = INFNM(1:L-1)//'n'
       L = L + 1
       ENDIF
       IF(CG.EQ.0.0) THEN
       INFNM = INFNM(1:L-1)//'p' 
       L = L + 1
       ENDIF
       outfnm = infnm(1:l-1)//'.be'
       outfnm2=infnm(1:l-1)//'.log'
       open(9,file=outfnm,status='replace')
       open(12,file=outfnm2,status='replace')

       call fdate(datetime)
       call cpu_time(tarray2)
       write(12,*)'Start date: ',datetime

c      Compute the centroid of the object and move the coordinates
c      to origin at the centroid.

       cx = 0.0
       cy = 0.0
       cz = 0.0
       do j =1,nar
       cx = cx + (v1(j,1)+v2(j,1)+v3(j,1))/3
       cy = cy + (v1(j,2)+v2(j,2)+v3(j,2))/3
       cz = cz + (v1(j,3)+v2(j,3)+v3(j,3))/3
       enddo

       cx = cx/nar
       cy = cy/nar
       cz = cz/nar
       do n =1, nar
       v1(n,1) = v1(n,1)-cx
       v1(n,2) = v1(n,2)-cy
       v1(n,3) = v1(n,3)-cz
       v2(n,1) = v2(n,1)-cx
       v2(n,2) = v2(n,2)-cy
       v2(n,3) = v2(n,3)-cz
       v3(n,1) = v3(n,1)-cx
       v3(n,2) = v3(n,2)-cy
       v3(n,3) = v3(n,3)-cz
       enddo

c      Calculate the areas of the triangles

       totarea=0.0
       do n =1,nar
       v21x = v2(n,1)-v1(n,1)
       v21y = v2(n,2)-v1(n,2)
       v21z = v2(n,3)-v1(n,3)
       v31x = v3(n,1)-v1(n,1)
       v31y = v3(n,2)-v1(n,2)
       v31z = v3(n,3)-v1(n,3)
       area(n) = DSQRT((v31y*v21z-v31z*v21y)**2 + 
     $ (v31z*v21x-v31x*v21z)**2+(v31x*v21y-v31y*v21x)**2)*0.5
       totarea = totarea+area(n)
       enddo

c      AREA CORRECTION	

       close(12)
       open(12,file=outfnm2,status='old')
       if (iarea.eq.1) then
       afac = truearea/totarea
       do j=1,nar
       area(j) = afac*area(j)
       enddo
       write(12,*)'BEST ver.',ver
       write(12,*)' '
       write(12,*)'control gtensor value=',cg
       write(12,*) 'true area (A^2)= ',truearea
       else
       afac= 1.0
       endif
       write(12,*) 'total area (A^2) = ',totarea
       write(12,*) 'centroid coords = ',cx,cy,cz 
       write(12,*) ' '
       close(12)

c      II. CONSTRUCTION OF THE SUPER "G" MATRIX, GR
c      The next call constructs the coefficients of the linear equations to be solved.

       CALL TGETGR(NAR,afac,cg,V1,V2,V3,X,AREA,GR)

       open(12,file=outfnm2,status='old',position='append')
       write(12,*) 'first few rows of gr matrix'
       do j=1,9
       write(12,*)gr(j,1),gr(j,2),gr(j,3),gr(j,4),gr(j,5)
       enddo
       write(12,*) 'first few columns of gr matrix'
       do j=1,9
       write(12,*)gr(1,j),gr(2,j),gr(3,j),gr(4,j),gr(5,j)
       enddo
       write(12,*) ' '

c      Deflation of the normal field.  On a closed surface the Stokes single layer of the
c      normal field vanishes: a uniform pressure makes no flow.  The unit normals are then a
c      near-null vector of GR (right), and XMU, the area times the unit normal of each
c      triangle, a left one.  The finer and smoother the surface, the closer GR comes to
c      singular, and DGETRS returns tractions dominated by that pressure mode.  The area
c      correction hides this in part: TGETGR scales the diagonal blocks by SQRT(AFAC) and the
c      others by AFAC, which moves GR's Rayleigh quotient along the normals to a few percent of
c      its mean diagonal entry, negative, for AFAC of 1.02 to 1.06.  With -na, or a surface
c      whose own area is exact (AFAC = 1), the quotient is 1e-4 or less of the diagonal.
c      The rigid-body and strain flows solved below are compatible with the near-null space:
c      XMU.U vanishes for their surface velocities U (up to the discretization error).  Adding
c      GDEF*XMU*XMU^T to GR then selects the solution with no net normal traction (XMU.F = 0)
c      and leaves it otherwise unchanged.  GDEF*|XMU|^2 is the mean diagonal entry of GR, so
c      the added term is of the size of GR's own entries, far above the area correction's
c      shift.  With GDEF ten times larger, D_t, CHI and [eta] change by 3e-8 or less and D_r
c      by 2e-7 or less.

       xmu2 = 0.0
       diagm = 0.0
       do n=1,nar
       v21x = v2(n,1)-v1(n,1)
       v21y = v2(n,2)-v1(n,2)
       v21z = v2(n,3)-v1(n,3)
       v31x = v3(n,1)-v1(n,1)
       v31y = v3(n,2)-v1(n,2)
       v31z = v3(n,3)-v1(n,3)
       xnx = v31y*v21z-v31z*v21y
       xny = v31z*v21x-v31x*v21z
       xnz = v31x*v21y-v31y*v21x
       xnl = DSQRT(xnx**2+xny**2+xnz**2)
       xmu(3*n-2) = area(n)*xnx/xnl
       xmu(3*n-1) = area(n)*xny/xnl
       xmu(3*n) = area(n)*xnz/xnl
       xmu2 = xmu2+xmu(3*n-2)**2+xmu(3*n-1)**2+xmu(3*n)**2
       diagm = diagm+gr(3*n-2,3*n-2)+gr(3*n-1,3*n-1)+gr(3*n,3*n)
       enddo
       diagm = diagm/nar3
       gdef = diagm/xmu2
!$OMP PARALLEL DO PRIVATE(I,GDJ) SCHEDULE(STATIC)
       DO J=1,NAR3
       GDJ = GDEF*XMU(J)
       DO I=1,NAR3
       GR(I,J) = GR(I,J)+GDJ*XMU(I)
       ENDDO
       ENDDO
!$OMP END PARALLEL DO
       write(12,*) 'normal-field deflation: mean diagonal of GR = ',diagm
       write(12,*) ' added GDEF*XMU*XMU^T to GR, GDEF = ',gdef
       write(12,*) ' '
       write(12,*) 'about to enter dgetrf'
       close(12)       

       ierr=0

c      III. SOLUTION OF SYSTEM OF EQUATIONS WITH LAPACK
c      Call the LAPACK Computational routine to factor GR into LU form.
C      The LU factored form is placed back into GR.

       call dgetrf(nar3,nar3,GR,nar3,ipiv,info)

       call fdate(datetime2)

       open(12,file=outfnm2,status='old',position='append')
       write(12,*)'Just out of dgetrf: ',datetime2
       write(12,*)' '
       close(12)

       call cpu_time(tarray3)
       write(12,*)'Elapsedtime upto Lapack solution:',tarray3-tarray2

c      Test the variable info and abort if there is a problem.

       if(info.ne.0) then
       if(info.lt.0) then
       print *,"Parameter ",Abs(info)," had an illegal value"
       print *, "Computation Aborted."
       STOP
       else
       print *,"Factor U(",info,info,")=0. U is singular."
       print *, "Computation Aborted."
       STOP
       endif
       endif

C      IV. CALCULATE THE DIFFUSION TENSORS

       CALL DIFCO(OUTFNM2,NAR,X,GR,AREA,
     $ INFO,DRR,DTT,DTR,CG,URC,UDR,TQ,B,F,IPIV)

C      V. CALCULATE THE VISCOSITY FACTOR IF REQUESTED

       IF(IVISCO.EQ.1) CALL VISCO(OUTFNM2,NAR,D1,D2,D3,
     $ CHIPRIME,AREA,X,ICOV,CG,GR,IPIV,INFO,G,H,W,V0)

C      VI. OUTPUT

       CALL OUTPUT(OUTFNM,VER,IAREA,TRUEAREA,AFAC,TRUEVOL,MW,
     $ WM,NAR,CX,CY,CZ,URC,UDR,DTT,DTR,DRR,IVISCO,ICOV,D1,D2,D3,
     $ CHIPRIME,CG,TOTAREA)

       call fdate(datetime3)

       open(12,file=outfnm2,status='old',position='append')
       write(12,*)' '
       write(12,*)'Program finished at:',datetime3
       call cpu_time(tarray4)
       write(12,*)'Elapsed time CPU_TIME:',tarray4-tarray2
       close(12)

       STOP
       END
