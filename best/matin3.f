      SUBROUTINE MATIN3(A,N,Y,IERR)
c Single precision version. S. Aragon 10-19-98.
C LOWER UPPER DECOMPOSITION METHOD OF INVERSE from Numerical Recipes.
C A is the input matrix which is destroyed, N is the dimension of A,
C Y is the inverse matrix and IERR is an error flag.

      PARAMETER (M=1)
      DIMENSION A(3*M,3*M),Y(3*M,3*M),INDX(3*M)
      DO 12 I = 1,N
      DO 11 J = 1,N
      Y(I,J) = 0.0
11    CONTINUE
      Y(I,I) = 1.0
12    CONTINUE

      CALL LUDCMP3(A,N,INDX,D,IERR)
      IF (IERR.EQ.1) RETURN
      DO 13 J = 1,N
      CALL LUBKSB3(A,N,INDX,Y(1,J))
13    CONTINUE
      RETURN
      END

      SUBROUTINE LUDCMP3(A,N,INDX,D,IERR)
      PARAMETER (M=1,TINY=1.D-20)
      DIMENSION A(3*M,3*M),INDX(3*M),VV(3*M)
      IERR = 0
      D = 1.0
      DO 12 I = 1,N
      AAMAX = 0.0
      DO 11 J = 1,N
      IF (DABS(A(I,J)).GT.AAMAX)  AAMAX = DABS(A(I,J))
11    CONTINUE
      IF (AAMAX.EQ.0.0) THEN
        IERR = 1
        AAMAX = 1.0
      ENDIF
      VV(I) = 1.0/AAMAX
12    CONTINUE
      IF (IERR.EQ.1) RETURN
      DO 19 J = 1,N
      IF (J.GT.1) THEN
        DO 14 I = 1,J-1
        SUM1 = A(I,J)
        IF (I.GT.1) THEN
          DO 13 K = 1,I-1
          SUM1 = SUM1 - A(I,K)*A(K,J)
13        CONTINUE
          A(I,J) = SUM1
        ENDIF
14      CONTINUE
      ENDIF
      AAMAX = 0.0
      DO 16 I = J,N
      SUM1 = A(I,J)
      IF (J.GT.1) THEN
        DO 15 K = 1,J-1
        SUM1 = SUM1 - A(I,K)*A(K,J)
15      CONTINUE
        A(I,J) = SUM1
      ENDIF
      DUM = VV(I)*DABS(SUM1)
      IF (DUM.GE.AAMAX) THEN
        IMAX = I
        AAMAX = DUM
      ENDIF
16    CONTINUE
      IF (J.NE.IMAX) THEN
        DO 17 K = 1,N
        DUM = A(IMAX,K)
        A(IMAX,K) = A(J,K)
        A(J,K) = DUM
17      CONTINUE
      D = -D
      VV(IMAX) = VV(J)
      ENDIF
      INDX(J) = IMAX
      IF (J.NE.N) THEN
        IF(A(J,J).EQ.0.0) A(J,J) = TINY
        DUM = 1.0/A(J,J)
        DO 18 I = J+1,N
        A(I,J) = A(I,J)*DUM
18      CONTINUE
      ENDIF
19    CONTINUE
      IF (A(N,N).EQ.0.0) A(N,N) = TINY
      RETURN
      END

      SUBROUTINE LUBKSB3(A,N,INDX,B)
      PARAMETER (M=1)
      DIMENSION A(3*M,3*M),INDX(3*M),B(3*M)
      II = 0
      DO 12 I = 1,N
      LL = INDX(I)
      SUM1 = B(LL)
      B(LL) = B(I)
      IF (II.NE.0) THEN
        DO 11 J = II,I-1
        SUM1 = SUM1 - A(I,J)*B(J)
11      CONTINUE
      ELSEIF (SUM1.NE.0.0) THEN
        II = 1
      ENDIF
      B(I) = SUM1
12    CONTINUE
      DO 14 I = N,1,-1
      SUM1 = B(I)
      IF (I.LT.N) THEN
         DO 13 J = I+1,N
         SUM1 = SUM1 - A(I,J)*B(J)
13       CONTINUE
      ENDIF
      B(I) = SUM1/A(I,I)
14    CONTINUE
      RETURN
      END

