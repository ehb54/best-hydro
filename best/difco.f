       SUBROUTINE DIFCO(OUTFNM2,NAR,X,GR,AREA,
     $ INFO,DRR,DTT,DTR,CG,URC,UDR,TQ,B,F,IPIV)

c      THE DIFFUSION TENSORS ARE DETERMINED IN THIS SUBROUTINE

       DIMENSION FTOTH(3),TTOTH(3),OMGA(3),D(3),V(3,3),U(3),
     $ DTR(3,3),DRR(3,3),DTT(3,3),UU(3,3),AUX(3,3),URC(3),UDR(3),
     $ GG(3,3),RR(3,3),X(NAR,3),GR(3*NAR,3*NAR),AREA(NAR),
     $ F(NAR,3),B(3*NAR,6),TQ(NAR,3),IPIV(3*NAR)
       CHARACTER*80 OUTFNM2
       REAL KTT(3,3),KTR(3,3),KRT(3,3),KRR(3,3),KTTINV(3,3)

       NAR3 = 3*NAR

C      BE MAIN LOOP: COMPUTE THE FORCES AND TORQUES FOR EACH FLOW
c      The flows correspond to three orthogonal translations, and three orthogonal
c      rotations.  The translation velocity is the U array and the angular velocity
c      is the array OMGA. The velocity has units of A/s, and the angular velocity
c      has units of 1/s. This assumes the triangule coordinates were input in A.

C      B(:,MX) contains the total body velocity for flow MX, composed of the sum
C      of the translational velocity u and the rotational velocity omga X x:
C      unit translations along x,y,z for MX=1-3, unit rotations about x,y,z for MX=4-6.

       DO N=1,NAR
       I0=3*(N-1)
       DO MX=1,6
       B(I0+1,MX)=0.0
       B(I0+2,MX)=0.0
       B(I0+3,MX)=0.0
       ENDDO
       B(I0+1,1)=1.0
       B(I0+2,2)=1.0
       B(I0+3,3)=1.0
       B(I0+2,4)=-x(n,3)
       B(I0+3,4)= x(n,2)
       B(I0+1,5)= x(n,3)
       B(I0+3,5)=-x(n,1)
       B(I0+1,6)=-x(n,2)
       B(I0+2,6)= x(n,1)
       enddo

C      Call the LAPACK computational routine to solve the linear system of equations
C      for all six flows at once, with the LU factors of GR.

       CALL DGETRS('N',NAR3,6,GR,NAR3,IPIV,B,NAR3,INFO)

       OPEN(12,FILE=OUTFNM2,STATUS='OLD',POSITION='APPEND')
       WRITE(12,*) ' **** TRANSLATION/ROTATION FLOWS ****'
       DO 100 MX=1,6
       OMGA(1)=0.0
       OMGA(2)=0.0
       OMGA(3)=0.0
       U(1)=0.0
       U(2)=0.0
       U(3)=0.0
       IF(MX .EQ. 1)U(1)=1.0
       IF(MX .EQ. 2)U(2)=1.0
       IF(MX .EQ. 3)U(3)=1.0
       IF(MX .EQ. 4)OMGA(1)=1.0
       IF(MX .EQ. 5)OMGA(2)=1.0
       IF(MX .EQ. 6)OMGA(3)=1.0
       WRITE(12,404)U(1),U(2),U(3)
404    FORMAT(4X,' U(X-Z) =',3F10.3,' A/SEC')
       WRITE(12,412)OMGA(1),OMGA(2),OMGA(3)
412    FORMAT(1X,' OMGA(X-Z) =',3F10.3,' RAD/S',/)

C      CALCULATE THE STRESS FORCE

       DO N = 1,NAR
       I0 = 3*(N-1)
       F(N,1) = B(I0+1,MX)
       F(N,2) = B(I0+2,MX)
       F(N,3) = B(I0+3,MX)
       TQ(N,1) = (X(N,2)*F(N,3)-X(N,3)*F(N,2))
       TQ(N,2) = (X(N,3)*F(N,1)-X(N,1)*F(N,3))
       TQ(N,3) = (X(N,1)*F(N,2)-X(N,2)*F(N,1))
       ENDDO

C      COMPUTE TOTAL FORCES/TORQUES ON PARTICLE

       FTOTH(1)=0.0
       FTOTH(2)=0.0
       FTOTH(3)=0.0
       TTOTH(1)=0.0
       TTOTH(2)=0.0
       TTOTH(3)=0.0
       DO 76 N=1,NAR
       FTOTH(1)=FTOTH(1)-AREA(N)*F(N,1)
       FTOTH(2)=FTOTH(2)-AREA(N)*F(N,2)
       FTOTH(3)=FTOTH(3)-AREA(N)*F(N,3)
       TTOTH(1)=TTOTH(1)-AREA(N)*TQ(N,1)
       TTOTH(2)=TTOTH(2)-AREA(N)*TQ(N,2)
       TTOTH(3)=TTOTH(3)-AREA(N)*TQ(N,3)
76     CONTINUE
       WRITE(12,78)FTOTH(1),FTOTH(2),FTOTH(3)
78     FORMAT(' FTOTH(X-Z) =',3E15.5)
       WRITE(12,90)TTOTH(1),TTOTH(2),TTOTH(3)
90     FORMAT(' TTOTH(X-Z) =',3E15.5,/)

C      CONSTRUCT COMPONENTS OF Ktt,Ktr,Krt,Krr FRICTION TENSORS

       IOP=0
       I=MX
       IF(MX .GT. 3)GO TO 119
       TEMP=1.0/U(I)
       GO TO 120
119    CONTINUE
       IOP=1
       I=MX-3
       TEMP=(1.0)/OMGA(I)
120    CONTINUE
       DO 122 J=1,3
       IF(IOP .EQ. 1)GO TO 124
       Ktt(J,I)=-FTOTH(J)*TEMP
       Ktr(J,I)=-TTOTH(J)*TEMP
       GO TO 122
124    CONTINUE
       Krt(J,I)=-FTOTH(J)*TEMP
       Krr(J,I)=-TTOTH(J)*TEMP
122    CONTINUE
100    CONTINUE

C      CYCLES COMPLETE
C      Write raw friction tensors

       DO  I=1,3
       DO  J=1,3
       WRITE(12,128)I,J,Ktt(I,J),Ktr(I,J),Krt(I,J),Krr(I,J)
128    FORMAT(' I,J,Ktt,Ktr,Krt,Krr=',2I5,4E12.4)
       enddo
       enddo

C      COMPUTE THE CENTER OF RESISTANCE/TRANSFORM TO CR
C      Fill the Center of Resistance matrix and vector.  Invert the matrix
C      and solve for the CR vector (stored in u).

       uu(1,1)=Ktt(2,2)+Ktt(3,3)
       uu(2,2)=Ktt(1,1)+Ktt(3,3)
       uu(3,3)=Ktt(2,2)+Ktt(1,1)
       uu(1,2)=-Ktt(1,2)
       uu(1,3)=-Ktt(1,3)
       uu(3,2)=-Ktt(2,3)
       uu(2,1)=uu(1,2)
       uu(2,3)=uu(3,2)
       uu(3,1)=uu(1,3)

C      Uniformitize Krt, Ktr (smooths out fluctuations)

       Krt(1,3) = (Krt(1,3)+Ktr(3,1))*0.5
       Krt(3,1) = (Krt(3,1)+Ktr(1,3))*0.5
       Krt(1,2) = (Krt(1,2)+Ktr(2,1))*0.5
       Krt(2,1) = (Krt(2,1)+Ktr(1,2))*0.5
       Krt(2,3) = (Krt(2,3)+Ktr(3,2))*0.5
       Krt(3,2) = (Krt(3,2)+Ktr(2,3))*0.5
       Krt(1,1) = (Krt(1,1)+Ktr(1,1))*0.5
       Krt(2,2) = (Krt(2,2)+Ktr(2,2))*0.5
       Krt(3,3) = (Krt(3,3)+Ktr(3,3))*0.5

       sur1= -Krt(3,2)+Krt(2,3)
       sur2= -Krt(1,3)+Krt(3,1)
       sur3= -Krt(2,1)+Krt(1,2)

       call matin3(uu,3,aux,ierr)
       if(ierr.ne.0) write(12,*)'Problem with C.R. inversion'

c      Vector between Center of Resistance/Centroid

       do  j=1,3
       uRC(j)= (aux(j,1)*sur1+aux(j,2)*sur2+aux(j,3)*sur3)
       enddo

c      Define the elements of antisymmetric matrix E.R 
c	This is where Jerome Fung caught error in 2020. uu matrix was previously set to zero.
C	Needed to use the above vector between Centers of Resistance/Centroid.

       uu(3,1)= -uRC(2)
       uu(3,2)= uRC(1)
       uu(1,2)= -uRC(3)
       uu(2,1)=-uu(1,2)
       uu(2,3)=-uu(3,2)
       uu(1,3)=-uu(3,1)
       uu(1,1)=0.0
       uu(2,2)=0.0
       uu(3,3)=0.0

c      Transform Krr, Krt to Center of Resistance. Krt should be symmetric..
C      Krt(j,i) is force j per unit rotation i.  With r = uRC and uu = [r]x,
C      Krt' = Krt + Ktt uu and Krr' = Krr - uu Ktt uu + Krt^T uu - uu Krt.
C      These used Krt where Krt^T belongs (Krt - uu Ktt, Krt uu - uu Krt^T),
C      which left Krt' unsymmetric.

       do i=1,3
       do j=1,3
       temp=0.0
       do k=1,3
       temp=temp+uu(i,k)*(Ktt(k,1)*uu(1,j)+Ktt(k,2)*uu(2,j)
     $ +Ktt(k,3)*uu(3,j))
       enddo
       temp1 = Krt(1,i)*uu(1,j)+Krt(2,i)*uu(2,j)+Krt(3,i)*uu(3,j)
       temp2 = uu(i,1)*Krt(1,j)+uu(i,2)*Krt(2,j)+uu(i,3)*Krt(3,j)
       Krr(i,j)=Krr(i,j)-temp +temp1-temp2 
       enddo
       enddo

       do i=1,3
       do j=1,3
       Krt(i,j)=Krt(i,j)+Ktt(i,1)*uu(1,j)+Ktt(i,2)*uu(2,j)
     $ +Ktt(i,3)*uu(3,j)
       enddo
       enddo
      
C      PRINT THE FRICTION TENSORS

       write(12,*) ' '
       write(12,*) ' Friction Tensors in Center of Resistance'
       write(12,*) ' '
       write(12,*) ' Ktt Tensor (A)'
       do j=1,3
       write(12,340)Ktt(j,1),Ktt(j,2),Ktt(j,3)
       end do
340    format(3(3x,e12.5))
       write(12,*) ' '
       write(12,*) ' Symmetric Krt Tensor (A^2)'
       do j=1,3
       write(12,340)Krt(j,1),Krt(j,2),Krt(j,3)
       end do
       write(12,*) ' '
       write(12,*) ' Krr Tensor (A^3)'
       do j=1,3
       write(12,340)Krr(j,1),Krr(j,2),Krr(j,3)
       end do

C      COMPUTE DIFFUSION TENSORS/TRANSFORM TO CM
C      Approximate computation of Dtt and Drr by ignoring translation-
C      rotation coupling.
C      Preserve Krr because matin destroys input.
C      Note: Use matin3 because the array sizes are all 3x3.
C      Approximate raw Dtt tensor

       CALL matin3(Ktt,3,Kttinv,ierr)
       if(ierr.ne.0) write(12,*) 'problem in inverting Ktt'
       write(12,*) ' '
       write(12,*) 'Raw Dtt tensor (no rot/trans coupling correction)'
       do j=1,3
       write(12,340) Kttinv(j,1),Kttinv(j,2),Kttinv(j,3)
       end do

C      Approximate raw Drr tensor
C      Save Krr in Ktt array before inverting.

       do j=1,3
       Ktt(1,j)=Krr(1,j)
       Ktt(2,j)=Krr(2,j)
       Ktt(3,j)=Krr(3,j) 
       enddo

       CALL matin3(Ktt,3,aux,ierr)
       if(ierr.ne.0) write(12,*) 'problem in inverting Krr'
       write(12,*) ' '
       write(12,*) 'Raw Drr tensor (no rot/trans coupling correction)'
       do j=1,3
       write(12,340) aux(j,1),aux(j,2),aux(j,3)
       end do

C      COMPUTE EXACT ROTATIONAL DIFFUSION TENSOR
C      Note:  Use Krt (not Ktr) in the remaining computations.  Ktr has not
C      been shifted to the Center of Resistance.
C      Array Ktt is now free for storage. Use it to hold intermediate values.

       DO I=1,3
       DO J=1,3
       temp = 0.0
       DO K=1,3
       temp = temp +Krt(i,k)*(Kttinv(k,1)*Krt(j,1)+
     $ Kttinv(k,2)*Krt(j,2)+Kttinv(k,3)*Krt(j,3))
       end do
       Ktt(I,J)=Krr(I,J)-temp
       end do
       end do

C      Invert the above to get Drr

       CALL matin3(Ktt,3,Drr,ierr)
       if(ierr.ne.0) write(12,*)'problem inverting modified Ktt for Drr'

C      MULTIPLY BY 4/3 FOR PRE-AVERAGED CALCULATION
       IF(CG.EQ.0) THEN
       DO I = 1,3
       DO J = 1,3
       DRR(I,J) = 4.0*DRR(I,J)/3.0
       ENDDO
       ENDDO
       ENDIF

C      Compute Dtr (to be used in Center of Diffusion Transformation)
C      Dtr = -Kttinv Krt^T Drr.  The sum over k used to keep only its
C      last term (k=3).

       do i=1,3
       do j=1,3
       temp = 0.0
       do k=1,3
       temp = temp +Kttinv(i,k)*(Krt(1,k)*Drr(1,j)+
     $ Krt(2,k)*Drr(2,j)+Krt(3,k)*Drr(3,j))
       end do
       dtr(i,j) =-temp
C      MULTIPLY BY 4/3 FOR PRE-AVERAGED CALCULATION
       IF(CG.EQ.0) DTR(I,J) = 4.0*DTR(I,J)/3.0
       end do
       end do

C      COMPUTE DTRANS EXACTLY

       DO  I=1,3
       DO  J=1,3
       temp = 0.0
       DO  K=1,3
C      Correct index order used below. The matrix dtr is not symmetric. 
       temp = temp +Kttinv(i,k)*(Krt(1,k)*dtr(j,1)+
     $ Krt(2,k)*dtr(j,2)+Krt(3,k)*dtr(j,3))
       enddo
       Dtt(I,J)=Kttinv(I,J)-temp
C      MULTIPLY BY 4/3 FOR PRE-AVERAGED CALCULATION
       IF(CG.EQ.0) DTT(I,J) = 4.0*DTT(I,J)/3.0
       enddo 
       enddo

       write(12,*) ' '
       write(12,*) 'Raw Dtt tensor, before CD transformation'
       do j=1,3
       write(12,340)Dtt(j,1),Dtt(j,2),Dtt(j,3)
       end do

C      We must now transform dtr, dtt to the center of diffusion. First find
C      the center of diffusion and then transform.  Use equations of
C      H. Brenner, J. Colloid Interface Science, 23, 407 (1967).
C      Fill the Center of Diffusion matrix and vector.  Invert the matrix
C      and solve for the CD vector (stored in u).

       uu(1,1)=drr(2,2)+drr(3,3)
       uu(2,2)=drr(1,1)+drr(3,3)
       uu(3,3)=drr(2,2)+drr(1,1)
       uu(1,2)=-drr(1,2)
       uu(1,3)=-drr(1,3)
       uu(3,2)=-drr(3,2)
       uu(2,1)=uu(1,2)
       uu(2,3)=uu(3,2)
       uu(3,1)=uu(1,3)

       sur1= dtr(3,2)-dtr(2,3)
       sur2= dtr(1,3)-dtr(3,1)
       sur3= dtr(2,1)-dtr(1,2)

       call matin3(uu,3,aux,ierr)
       if(ierr.ne.0) write(12,*)'Problem with C.D. inversion'

C      Vector between Centers of Diffusion/Resistance

       do  j=1,3
       uDR(j)= (aux(j,1)*sur1+aux(j,2)*sur2+aux(j,3)*sur3)
       enddo

C      Define the elements of antisymmetric matrix E.R
C	This is where Jerome Fung found error  in 2020. uu matrix was being set to zero.
C	Needed to use the above vector between Centers of Diffusion/Resistance.

       uu(3,1)= -uDR(2)
       uu(3,2)= uDR(1)
       uu(1,2)= -uDR(3)
       uu(2,1)=-uu(1,2)
       uu(2,3)=-uu(3,2)
       uu(1,3)=-uu(3,1)
       uu(1,1)=0.0
       uu(2,2)=0.0
       uu(3,3)=0.0

C      Transform dtr, dtt to center of diffusion.

       do i=1,3
       do j=1,3
       temp=0.0
       do k=1,3
       temp=temp+uu(i,k)*drr(k,1)*uu(1,j)+uu(i,k)*drr(k,2)*uu(2,j)
     $ +uu(i,k)*drr(k,3)*uu(3,j)
       enddo
       temp1 = Dtr(i,1)*uu(1,j)+Dtr(i,2)*uu(2,j)+Dtr(i,3)*uu(3,j)
       temp2 =-uu(i,1)*Dtr(j,1)-uu(i,2)*Dtr(j,2)-uu(i,3)*Dtr(j,3)
       dtt(i,j)=dtt(i,j)-temp+temp1+temp2
       enddo
       enddo

       do i=1,3
       do j=1,3
       dtr(i,j)=dtr(i,j)-uu(i,1)*drr(1,j)-uu(i,2)*drr(2,j)
     $ -uu(i,3)*drr(3,j)
       enddo
       enddo

       CLOSE(12)

       RETURN
       END
