       SUBROUTINE VISCO(OUTFNM2,NAR,D1,D2,D3,
     $ CHIPRIME,AREA,X,ICOV,CG,GR,IPIV,INFO,G,H,W,V0)

C      THIS SUBROUTINE CALCULATES THE VISCOSITY FACTOR (S. ALLISON,
C      MACROMOLECULES 1998,31,4464-4474.) AT EITHER THE VISCOSITY CENTER
C      OR THE CENTROID 

       DIMENSION BB(3,3),BBI(3,3),TEMP01(5),TEMP02(5),TEMP03(5),
     $ TEMP11(5),TEMP12(5),TEMP13(5),TEMP21(5),TEMP22(5),TEMP23(5),
     $ TEMP31(5),TEMP32(5),TEMP33(5),GR(3*NAR,3*NAR),AREA(NAR),
     $ G(NAR,5,3),H(NAR,5,3,3),IPIV(3*NAR),X(NAR,3),W(5,3*NAR),
     $ V0(3*NAR) 
       CHARACTER*80 OUTFNM2
       
       PI = 3.1415926535897932385
       NAR3 = 3*NAR

C      DETERMINE THE VISCOSITY CENTER, IF REQUESTED
       IF(ICOV.EQ.0) GOTO 600

       DO N1 = 1,NAR
       I0 = 3*(N1-1)
 
       DO L = 1,5
       TEMP01(L) = 0.0
       TEMP02(L) = 0.0
       TEMP03(L) = 0.0
       TEMP11(L) = 0.0
       TEMP12(L) = 0.0
       TEMP13(L) = 0.0
       TEMP21(L) = 0.0
       TEMP22(L) = 0.0
       TEMP23(L) = 0.0
       TEMP31(L) = 0.0
       TEMP32(L) = 0.0
       TEMP33(L) = 0.0
       ENDDO
 
       DO N2 = 1,NAR
       J0 = 3*(N2-1)
       GRE11 = GR(I0+1,J0+2)
       GRE12 = GR(I0+1,J0+1)
       GRE21 = GR(I0+2,J0+2)
       GRE22 = GR(I0+2,J0+1)
       GRE31 = GR(I0+3,J0+2)
       GRE32 = GR(I0+3,J0+1)
       TEMP01(1) = TEMP01(1) + GRE11*X(N2,1) + GRE12*X(N2,2)
       TEMP02(1) = TEMP02(1) + GRE21*X(N2,1) + GRE22*X(N2,2)
       TEMP03(1) = TEMP03(1) + GRE31*X(N2,1) + GRE32*X(N2,2)
       TEMP11(1) = TEMP11(1) + GRE11
       TEMP12(1) = TEMP12(1) + GRE12
       TEMP21(1) = TEMP21(1) + GRE21
       TEMP22(1) = TEMP22(1) + GRE22
       TEMP31(1) = TEMP31(1) + GRE31
       TEMP32(1) = TEMP32(1) + GRE32
 
       GRE11 = GR(I0+1,J0+3)
       GRE13 = GR(I0+1,J0+1)
       GRE21 = GR(I0+2,J0+3)
       GRE23 = GR(I0+2,J0+1)
       GRE31 = GR(I0+3,J0+3)
       GRE33 = GR(I0+3,J0+1)
       TEMP01(2) = TEMP01(2) + GRE11*X(N2,1) + GRE13*X(N2,3)
       TEMP02(2) = TEMP02(2) + GRE21*X(N2,1) + GRE23*X(N2,3)
       TEMP03(2) = TEMP03(2) + GRE31*X(N2,1) + GRE33*X(N2,3)
       TEMP11(2) = TEMP11(2) + GRE11
       TEMP13(2) = TEMP13(2) + GRE13
       TEMP21(2) = TEMP21(2) + GRE21
       TEMP23(2) = TEMP23(2) + GRE23
       TEMP31(2) = TEMP31(2) + GRE31
       TEMP33(2) = TEMP33(2) + GRE33
 
       GRE12 = GR(I0+1,J0+3)
       GRE13 = GR(I0+1,J0+2)
       GRE22 = GR(I0+2,J0+3)
       GRE23 = GR(I0+2,J0+2)
       GRE32 = GR(I0+3,J0+3)
       GRE33 = GR(I0+3,J0+2)
       TEMP01(3) = TEMP01(3) + GRE12*X(N2,2) + GRE13*X(N2,3)
       TEMP02(3) = TEMP02(3) + GRE22*X(N2,2) + GRE23*X(N2,3)
       TEMP03(3) = TEMP03(3) + GRE32*X(N2,2) + GRE33*X(N2,3)
       TEMP12(3) = TEMP12(3) + GRE12
       TEMP13(3) = TEMP13(3) + GRE13
       TEMP22(3) = TEMP22(3) + GRE22
       TEMP23(3) = TEMP23(3) + GRE23
       TEMP32(3) = TEMP32(3) + GRE32
       TEMP33(3) = TEMP33(3) + GRE33
 
       GRE11 = GR(I0+1,J0+1)
       GRE12 = -GR(I0+1,J0+2)
       GRE21 = GR(I0+2,J0+1)
       GRE22 = -GR(I0+2,J0+2)
       GRE31 = GR(I0+3,J0+1)
       GRE32 = -GR(I0+3,J0+2)
       TEMP01(4) = TEMP01(4) + GRE11*X(N2,1) + GRE12*X(N2,2)
       TEMP02(4) = TEMP02(4) + GRE21*X(N2,1) + GRE22*X(N2,2)
       TEMP03(4) = TEMP03(4) + GRE31*X(N2,1) + GRE32*X(N2,2)
       TEMP11(4) = TEMP11(4) + GRE11
       TEMP12(4) = TEMP12(4) + GRE12
       TEMP21(4) = TEMP21(4) + GRE21
       TEMP22(4) = TEMP22(4) + GRE22
       TEMP31(4) = TEMP31(4) + GRE31
       TEMP32(4) = TEMP32(4) + GRE32
 
       GRE11 = GR(I0+1,J0+1)
       GRE13 = -GR(I0+1,J0+3)
       GRE21 = GR(I0+2,J0+1)
       GRE23 = -GR(I0+2,J0+3)
       GRE31 = GR(I0+3,J0+1)
       GRE33 = -GR(I0+3,J0+3)
       TEMP01(5) = TEMP01(5) + GRE11*X(N2,1) + GRE13*X(N2,3)
       TEMP02(5) = TEMP02(5) + GRE21*X(N2,1) + GRE23*X(N2,3)
       TEMP03(5) = TEMP03(5) + GRE31*X(N2,1) + GRE33*X(N2,3)
       TEMP11(5) = TEMP11(5) + GRE11
       TEMP13(5) = TEMP13(5) + GRE13
       TEMP21(5) = TEMP21(5) + GRE21
       TEMP23(5) = TEMP23(5) + GRE23
       TEMP31(5) = TEMP31(5) + GRE31
       TEMP33(5) = TEMP33(5) + GRE33
       ENDDO                     
 
       DO L = 1,5
       G(N1,L,1) = TEMP01(L)
       G(N1,L,2) = TEMP02(L)
       G(N1,L,3) = TEMP03(L)
       H(N1,L,1,1) = TEMP11(L)
       H(N1,L,1,2) = TEMP12(L)
       H(N1,L,1,3) = TEMP13(L)
       H(N1,L,2,1) = TEMP21(L)
       H(N1,L,2,2) = TEMP22(L)
       H(N1,L,2,3) = TEMP23(L)
       H(N1,L,3,1) = TEMP31(L)
       H(N1,L,3,2) = TEMP32(L)
       H(N1,L,3,3) = TEMP33(L)
       ENDDO
       ENDDO
 
       BB(1,1) = 0.0
       BB(1,2) = 0.0
       BB(1,3) = 0.0
       BB(2,1) = 0.0
       BB(2,2) = 0.0
       BB(2,3) = 0.0
       BB(3,1) = 0.0
       BB(3,2) = 0.0
       BB(3,3) = 0.0
       C1 = 0.0
       C2 = 0.0
       C3 = 0.0
 
       DO N = 1,NAR
       BB(1,1) = BB(1,1) + AREA(N)*(3*(H(N,1,2,1) + H(N,2,3,1)) 
     $ + 2*(H(N,4,1,1) + H(N,5,1,1)))
 
       BB(1,2) = BB(1,2) + AREA(N)*(-4*H(N,4,2,1)
     $ + 3*(H(N,1,1,1) + H(N,1,2,2) + H(N,2,3,2) + H(N,3,3,1))
     $ + 2*(H(N,4,1,2) + H(N,5,1,2) + H(N,5,2,1)))
 
       BB(1,3) = BB(1,3) + AREA(N)*(-4*H(N,5,3,1)
     $ + 3*(H(N,1,2,3) + H(N,2,1,1) + H(N,2,3,3) + H(N,3,2,1))
     $ + 2*(H(N,4,1,3) + H(N,4,3,1) + H(N,5,1,3)))
 
       BB(2,2) = BB(2,2) + AREA(N)*(-4*H(N,4,2,2)
     $ + 3*(H(N,1,1,2) + H(N,3,3,2)) + 2*H(N,5,2,2))
 
       BB(2,3) = BB(2,3) + AREA(N)*(-4*H(N,4,2,3)
     $ + 3*(H(N,3,2,2) + H(N,3,3,3) + H(N,1,1,3) + H(N,2,1,2))
     $ + 2*(H(N,4,3,2) + H(N,5,2,3)))
 
       BB(3,3) = BB(3,3) + AREA(N)*(-4*H(N,5,3,3)
     $ + 3*(H(N,2,1,3) + H(N,3,2,3)) + 2*H(N,4,3,3))
 
       C1 = C1 + AREA(N)*
     $ (-3*G(N,1,2) - 3*G(N,2,3) - 2*G(N,4,1) - 2*G(N,5,1)
     $ -X(N,1)*(3*(H(N,1,2,1)+H(N,2,3,1))+2*(H(N,4,1,1)+H(N,5,1,1)))
     $ -X(N,2)*(3*(H(N,1,1,1)+H(N,3,3,1))-4*H(N,4,2,1)+2*H(N,5,2,1))
     $ -X(N,3)*(3*(H(N,2,1,1)+H(N,3,2,1))-4*H(N,5,3,1)+2*H(N,4,3,1)))
 
       C2 = C2 + AREA(N)*
     $ (-3*G(N,1,1) - 3*G(N,3,3) + 4*G(N,4,2) - 2*G(N,5,2)
     $ -X(N,1)*(3*(H(N,1,2,2)+H(N,2,3,2))+2*(H(N,4,1,2)+H(N,5,1,2)))
     $ -X(N,2)*(3*(H(N,1,1,2)+H(N,3,3,2))-4*H(N,4,2,2)+2*H(N,5,2,2))
     $ -X(N,3)*(3*(H(N,2,1,2)+H(N,3,2,2))-4*H(N,5,3,2)+2*H(N,4,3,2)))
 
       C3 = C3 + AREA(N)*
     $ (-3*G(N,2,1) - 3*G(N,3,2) - 2*G(N,4,3) + 4*G(N,5,3)
     $ -X(N,1)*(3*(H(N,1,2,3)+H(N,2,3,3))+2*(H(N,4,1,3)+H(N,5,1,3)))
     $ -X(N,2)*(3*(H(N,1,1,3)+H(N,3,3,3))-4*H(N,4,2,3)+2*H(N,5,2,3))
     $ -X(N,3)*(3*(H(N,2,1,3)+H(N,3,2,3))-4*H(N,5,3,3)+2*H(N,4,3,3)))
 
       ENDDO
 
       BB(1,1) = BB(1,1)/15.
       BB(1,2) = BB(1,2)/30.
       BB(1,3) = BB(1,3)/30.
       BB(2,1) = BB(1,2)
       BB(2,2) = BB(2,2)/15.
       BB(2,3) = BB(2,3)/30.
       BB(3,1) = BB(1,3)
       BB(3,2) = BB(2,3)
       BB(3,3) = BB(3,3)/15.
       C1 = C1/30.
       C2 = C2/30.
       C3 = C3/30.

       IERR = 0  
       CALL MATIN3(BB,3,BBI,IERR)
       IF(IERR.NE.0) THEN
       OPEN(12,FILE=OUTFNM2,STATUS='OLD')
       WRITE(12,*) 'PROBLEM IN INVERTING BB'
       CLOSE(12)
       ENDIF
 
       D1 = BBI(1,1)*C1 + BBI(1,2)*C2 + BBI(1,3)*C3
       D2 = BBI(2,1)*C1 + BBI(2,2)*C2 + BBI(2,3)*C3
       D3 = BBI(3,1)*C1 + BBI(3,2)*C2 + BBI(3,3)*C3

C      REFER INCENTER COORDINATES TO THE CENTER OF VISCOSITY

       DO N = 1,NAR
       X(N,1) = X(N,1) - D1
       X(N,2) = X(N,2) - D2
       X(N,3) = X(N,3) - D3
       ENDDO

600    IF(ICOV.EQ.1) THEN

C      COMPUTE THE STRESS FORCE AT THE CENTER OF VISCOSITY

       DO N1 = 1,NAR
       I0 = 3*(N1-1)
       TEMPA = 0.0
       TEMPB = 0.0
       TEMPC = 0.0
       TEMPD = 0.0
       TEMPE = 0.0
       TEMPF = 0.0
       TEMPG = 0.0
       TEMPH = 0.0
       TEMPI = 0.0
       TEMPJ = 0.0
       TEMPK = 0.0
       TEMPL = 0.0
       TEMPM = 0.0
       TEMPN = 0.0
       TEMPO = 0.0
       DO N2 = 1,NAR
       J0 = 3*(N2-1)
       TEMPA = TEMPA + (X(N2,2)*GR(I0+1,J0+1) + X(N2,1)*GR(I0+1,J0+2))
       TEMPB = TEMPB + (X(N2,2)*GR(I0+2,J0+1) + X(N2,1)*GR(I0+2,J0+2))
       TEMPC = TEMPC + (X(N2,2)*GR(I0+3,J0+1) + X(N2,1)*GR(I0+3,J0+2))
       TEMPD = TEMPD + (X(N2,3)*GR(I0+1,J0+1) + X(N2,1)*GR(I0+1,J0+3))
       TEMPE = TEMPE + (X(N2,3)*GR(I0+2,J0+1) + X(N2,1)*GR(I0+2,J0+3))
       TEMPF = TEMPF + (X(N2,3)*GR(I0+3,J0+1) + X(N2,1)*GR(I0+3,J0+3))
       TEMPG = TEMPG + (X(N2,3)*GR(I0+1,J0+2) + X(N2,2)*GR(I0+1,J0+3))
       TEMPH = TEMPH + (X(N2,3)*GR(I0+2,J0+2) + X(N2,2)*GR(I0+2,J0+3))
       TEMPI = TEMPI + (X(N2,3)*GR(I0+3,J0+2) + X(N2,2)*GR(I0+3,J0+3))
       TEMPJ = TEMPJ + (X(N2,1)*GR(I0+1,J0+1) - X(N2,2)*GR(I0+1,J0+2))
       TEMPK = TEMPK + (X(N2,1)*GR(I0+2,J0+1) - X(N2,2)*GR(I0+2,J0+2))
       TEMPL = TEMPL + (X(N2,1)*GR(I0+3,J0+1) - X(N2,2)*GR(I0+3,J0+2))
       TEMPM = TEMPM + (X(N2,1)*GR(I0+1,J0+1) - X(N2,3)*GR(I0+1,J0+3))
       TEMPN = TEMPN + (X(N2,1)*GR(I0+2,J0+1) - X(N2,3)*GR(I0+2,J0+3))
       TEMPO = TEMPO + (X(N2,1)*GR(I0+3,J0+1) - X(N2,3)*GR(I0+3,J0+3))
       ENDDO
       W(1,I0+1) = 0.5*TEMPA
       W(1,I0+2) = 0.5*TEMPB
       W(1,I0+3) = 0.5*TEMPC
       W(2,I0+1) = 0.5*TEMPD
       W(2,I0+2) = 0.5*TEMPE
       W(2,I0+3) = 0.5*TEMPF
       W(3,I0+1) = 0.5*TEMPG
       W(3,I0+2) = 0.5*TEMPH
       W(3,I0+3) = 0.5*TEMPI
       W(4,I0+1) = 0.5*TEMPJ
       W(4,I0+2) = 0.5*TEMPK
       W(4,I0+3) = 0.5*TEMPL
       W(5,I0+1) = 0.5*TEMPM
       W(5,I0+2) = 0.5*TEMPN
       W(5,I0+3) = 0.5*TEMPO
       ENDDO
       ELSE

C      COMPUTE THE STRESS FORCE AT THE CENTROID

       DO L = 1,5
       IF(L.EQ.1) THEN
       DO N = 1,NAR
       I0=3*(N-1)
       V0(I0+1) = 0.5*X(N,2)
       V0(I0+2) = 0.5*X(N,1)
       V0(I0+3) = 0.0
       ENDDO 
       ENDIF
       IF(L.EQ.2) THEN
       DO N = 1,NAR
       I0=3*(N-1)
       V0(I0+1) = 0.5*X(N,3)
       V0(I0+2) = 0.0
       V0(I0+3) = 0.5*X(N,1)
       ENDDO
       ENDIF
       IF(L.EQ.3) THEN
       DO N = 1,NAR
       I0=3*(N-1)
       V0(I0+1) = 0.0
       V0(I0+2) = 0.5*X(N,3)
       V0(I0+3) = 0.5*X(N,2)
       ENDDO
       ENDIF
       IF(L.EQ.4) THEN
       DO N = 1,NAR
       I0=3*(N-1)
       V0(I0+1) = 0.5*X(N,1)
       V0(I0+2) =-0.5*X(N,2)
       V0(I0+3) = 0.0
       ENDDO
       ENDIF
       IF(L.EQ.5) THEN
       DO N = 1,NAR
       I0=3*(N-1)
       V0(I0+1) = 0.5*X(N,1)
       V0(I0+2) = 0.0
       V0(I0+3) =-0.5*X(N,3)
       ENDDO
       ENDIF
       CALL DGETRS('N',NAR3,1,GR,NAR3,IPIV,V0,NAR3,INFO)
       DO N = 1,NAR
       I0=3*(N-1)
       W(L,I0+1) = V0(I0+1)
       W(L,I0+2) = V0(I0+2)
       W(L,I0+3) = V0(I0+3)
       ENDDO
       ENDDO
       ENDIF

C      COMPUTE THE VISCOSITY FACTOR

       CHI112 = 0.0
       CHI213 = 0.0 
       CHI323 = 0.0
       CHI411 = 0.0
       CHI422 = 0.0
       CHI433 = 0.0
       CHI511 = 0.0
       CHI522 = 0.0
       CHI533 = 0.0
       DO N = 1,NAR
       I0=3*(N-1)
       CHI112 = CHI112 + AREA(N)*(W(1,I0+1)*X(N,2) + W(1,I0+2)*X(N,1))
       CHI213 = CHI213 + AREA(N)*(W(2,I0+1)*X(N,3) + W(2,I0+3)*X(N,1)) 
       CHI323 = CHI323 + AREA(N)*(W(3,I0+2)*X(N,3) + W(3,I0+3)*X(N,2)) 
       CHI411 = CHI411 + AREA(N)*W(4,I0+1)*X(N,1)
       CHI422 = CHI422 + AREA(N)*W(4,I0+2)*X(N,2)
       CHI433 = CHI433 + AREA(N)*W(4,I0+3)*X(N,3)
       CHI511 = CHI511 + AREA(N)*W(5,I0+1)*X(N,1)
       CHI522 = CHI522 + AREA(N)*W(5,I0+2)*X(N,2)
       CHI533 = CHI533 + AREA(N)*W(5,I0+3)*X(N,3)
       ENDDO
       CHI112 = 0.5*CHI112
       CHI213 = 0.5*CHI213
       CHI323 = 0.5*CHI323 
 
       CHIPRIME = 8*PI*(0.2*(CHI112 + CHI213 + CHI323) +   
     $ (CHI411 - 2*CHI422 + CHI433 + CHI511 + CHI522 - 2*CHI533)/15.0)

C      MULTIPLY BY 0.75 FOR PREAVERAGED CALCULATION
       IF(CG.EQ.0) CHIPRIME = 0.75*CHIPRIME

       RETURN
       END
