       SUBROUTINE VISCO(OUTFNM2,NAR,D1,D2,D3,
     $ CHIPRIME,AREA,X,ICOV,CG,GR,IPIV,INFO,G,H,W,V0)

C      THIS SUBROUTINE CALCULATES THE VISCOSITY FACTOR (S. ALLISON,
C      MACROMOLECULES 1998,31,4464-4474.) AT EITHER THE VISCOSITY CENTER
C      OR THE CENTROID 

       DIMENSION BB(3,3),BBI(3,3),GR(3*NAR,3*NAR),AREA(NAR),
     $ G(NAR,5,3),H(NAR,5,3,3),IPIV(3*NAR),X(NAR,3),W(5,3*NAR),
     $ V0(3*NAR,8) 
       CHARACTER*80 OUTFNM2
       
       PI = 3.1415926535897932385
       NAR3 = 3*NAR

C      DETERMINE THE VISCOSITY CENTER, IF REQUESTED
       IF(ICOV.EQ.0) GOTO 600

C      G AND H HOLD THE PRODUCTS OF THE INVERSE OF GR WITH THE VELOCITIES OF
C      THE FIVE PURE STRAIN FLOWS (TIMES 2) AND OF THE THREE UNIT TRANSLATIONS,
C      AT THE CENTROID.  INSTEAD OF INVERTING GR, SOLVE FOR ALL EIGHT FLOWS AT
C      ONCE WITH THE LU FACTORS OF GR: V0(:,1-5) THE STRAIN FLOWS, V0(:,6-8)
C      THE TRANSLATIONS ALONG X,Y,Z.

       CALL STRAIN(NAR,X,V0)
       DO N = 1,NAR
       I0 = 3*(N-1)
       DO K = 1,3
       V0(I0+K,6) = 0.0
       V0(I0+K,7) = 0.0
       V0(I0+K,8) = 0.0
       ENDDO
       V0(I0+1,6) = 1.0
       V0(I0+2,7) = 1.0
       V0(I0+3,8) = 1.0
       ENDDO
       CALL DGETRS('N',NAR3,8,GR,NAR3,IPIV,V0,NAR3,INFO)

C      G(N,L,K): COMPONENT K AT TRIANGLE N FOR STRAIN FLOW L.
C      H(N,L,K,M): THE CHANGE IN G(N,L,K) PER UNIT SHIFT OF ALL COORDINATES
C      X(:,M); THE DERIVATIVE OF FLOW L'S VELOCITY WITH RESPECT TO X(:,M) IS
C      PLUS OR MINUS A UNIT TRANSLATION, OR ZERO.

       DO N1 = 1,NAR
       I0 = 3*(N1-1)
       DO K = 1,3
       DO L = 1,5
       G(N1,L,K) = 2.0*V0(I0+K,L)
       ENDDO
       E1 = V0(I0+K,6)
       E2 = V0(I0+K,7)
       E3 = V0(I0+K,8)
       H(N1,1,K,1) = E2
       H(N1,1,K,2) = E1
       H(N1,1,K,3) = 0.0
       H(N1,2,K,1) = E3
       H(N1,2,K,2) = 0.0
       H(N1,2,K,3) = E1
       H(N1,3,K,1) = 0.0
       H(N1,3,K,2) = E3
       H(N1,3,K,3) = E2
       H(N1,4,K,1) = E1
       H(N1,4,K,2) = -E2
       H(N1,4,K,3) = 0.0
       H(N1,5,K,1) = E1
       H(N1,5,K,2) = 0.0
       H(N1,5,K,3) = -E3
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

C      COMPUTE THE STRESS FORCE AT THE CENTER OF VISCOSITY (X WAS SHIFTED THERE
C      ABOVE) OR AT THE CENTROID.  SOLVE FOR THE FIVE PURE STRAIN FLOWS AT ONCE
C      WITH THE LU FACTORS OF GR.

600    CALL STRAIN(NAR,X,V0)
       CALL DGETRS('N',NAR3,5,GR,NAR3,IPIV,V0,NAR3,INFO)
       DO L = 1,5
       DO N = 1,NAR
       I0=3*(N-1)
       W(L,I0+1) = V0(I0+1,L)
       W(L,I0+2) = V0(I0+2,L)
       W(L,I0+3) = V0(I0+3,L)
       ENDDO
       ENDDO

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

       SUBROUTINE STRAIN(NAR,X,V0)

C      VELOCITIES OF THE FIVE PURE STRAIN FLOWS AT THE INCENTERS X, IN
C      V0(:,L), L = 1,...,5 (S. ALLISON, MACROMOLECULES 1998,31,4464-4474.)

       DIMENSION X(NAR,3),V0(3*NAR,5)

       DO N = 1,NAR
       I0=3*(N-1)
       V0(I0+1,1) = 0.5*X(N,2)
       V0(I0+2,1) = 0.5*X(N,1)
       V0(I0+3,1) = 0.0
       V0(I0+1,2) = 0.5*X(N,3)
       V0(I0+2,2) = 0.0
       V0(I0+3,2) = 0.5*X(N,1)
       V0(I0+1,3) = 0.0
       V0(I0+2,3) = 0.5*X(N,3)
       V0(I0+3,3) = 0.5*X(N,2)
       V0(I0+1,4) = 0.5*X(N,1)
       V0(I0+2,4) =-0.5*X(N,2)
       V0(I0+3,4) = 0.0
       V0(I0+1,5) = 0.5*X(N,1)
       V0(I0+2,5) = 0.0
       V0(I0+3,5) =-0.5*X(N,3)
       ENDDO

       RETURN
       END
