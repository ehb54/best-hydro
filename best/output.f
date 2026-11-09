       SUBROUTINE OUTPUT(OUTFNM,VER,IAREA,TRUEAREA,AFAC,TRUEVOL,MW,
     $ WM,NAR,CX,CY,CZ,URC,UDR,DTT,DTR,DRR,IVISCO,ICOV,D1,D2,D3,
     $ CHIPRIME,CG,TOTAREA)

C      OUTPUT TO .be

       DIMENSION DTR(3,3),DRR(3,3),DTT(3,3),URC(3),UDR(3),D(3),V(3,3)
       character*80 outfnm
       character*16 ver

       IOUNIT = 9

C      PRINT AREA AND ASSUMED FLUID PARAMETERS (water, 20 C)
c      Values of physical constants. "cgs" units are used. Distances in A, assumed.

       temperature = 293.15
       eta = 0.01002
       pi = 3.141592653
       tunit = (1.38066E-16)*temperature*1.0E8/(8.0*eta*pi)
       runit = tunit*1.0e16
       rtunit = tunit*1.0e8

c      The units of "unit" quantities are: tunit: cm^2 A/s, rtunit: cm A^2/s, runit: A^3/s.
c      The internal units of the diffusion tensors are: Dtt:1/A, Drt: 1/A^2, Drr: 1/A^3.
c      The "cgs" units of the diffusion tensors are: Dtt: cm^2/s, Drt: cm/s, Drr: 1/s.

       ftunit = 8.0*eta*pi*1.0e8
       ftrunit = ftunit*1.0e8
       frunit = ftunit*1.0e16

c      The units of "unit" quantities are: ftunit: g/A/s, frtunit: g cm/A^2/s, frunit: g cm^2/A^3/s.
c      The internal units of the friction tensors are: Ktt: A, Krt: A^2, Krr: A^3.
c      The "cgs" units of the friction tensors are: Ktt: g/s, Krt: g cm/s, Krr: g cm^2/s.

ccccc  OPEN(9,FILE=OUTFNM,Status='old',position='append')
       OPEN(9,FILE=OUTFNM,Status='old')

       WRITE(iounit,*)'BEST Ver. ',ver 
       WRITE(iounit,400)
400    FORMAT(/,' BE - MACROMOLECULE TRANSPORT TENSORS',/)
       WRITE(iounit,403)totarea
403    FORMAT(4X,' Area (A^2)            ',E12.5)
       if (iarea.eq.1)  then
       write(iounit,414) truearea
414    format(4x,' True Area (A^2)       ',E12.5)
       write(iounit,415) afac
415    format(4x,' Areafactor                 ',f7.4)
       WRITE(IOUNIT,416) TRUEVOL
416    FORMAT(4X,' True Volume (A^3)     ',E12.5)
       endif
       IF (MW.EQ.1) THEN
       WRITE(IOUNIT,420) WM
420    FORMAT(4X,' Molecular Weight (g/mol)  ',F8.1)
       ENDIF
       WRITE(iounit,401)temperature
401    FORMAT(4X,' Temperature (K)           ',F8.3)
       WRITE(iounit,402)eta
402    FORMAT(4X,' Solvent Viscosity (Poise) ',F8.5)
       WRITE(iounit,410)NAR
410    FORMAT(4X,' Number of Surface Patches    ',I5)

       write(iounit,*) ' '                                    
       IF(CG.EQ.0) THEN                                       
       WRITE(IOUNIT,*) " PRE-AVERAGING USED FOR THIS CALCULATION"
       ENDIF                                                  

C      PRINT THE DIFFUSION TENSORS

       write(iounit,*) ' '
       write(iounit,*) " Centroid Coordinates (A)"
       write(iounit,340) CX,CY,CZ
       write(iounit,*) " Vector between Center of Resistance/Centroid "
       write(iounit,340) URC(1),URC(2),URC(3)
       write(iounit,*) " Vector between Centers of Diffusion/Resistance"
       write(iounit,340) UDR(1),UDR(2),UDR(3)

       write(iounit,*) ' '
       write(iounit,*) ' Diffusion Tensors in Center of Mobility'
       write(iounit,*) ' '
       write(iounit,*) ' Dtt Tensor (1/A)'
       do j=1,3
       write(iounit,340)dtt(j,1),dtt(j,2),dtt(j,3)
340    format(3(3x,e12.5))
       enddo
       write(iounit,*) ' '
       write(iounit,*) ' Symmetric Dtr Tensor (1/A^2)'
       do j=1,3
       write(iounit,340)dtr(j,1),dtr(j,2),dtr(j,3)
       enddo
       write(iounit,*) ' '
       write(iounit,*) '  Drr TENSOR (1/A^3)'
       do j=1,3
       write(iounit,340)drr(j,1),drr(j,2),drr(j,3)
       end do

c      Diagonalize the tensors.

       call jacobi(drr,3,3,d,v,nrot)
       write(iounit,*) ' '
       write(iounit,*) 'EIGENVALUES OF Drr TENSOR (1/A^3)'
       write(iounit,340)d(1),d(2),d(3)
       trace=(d(1)+d(2)+d(3))/3.
       write(iounit,350)trace
       del=DSQRT(DABS(d(1)*d(1)+d(2)*d(2)+d(3)*d(3)
     $ -d(1)*d(2)-d(1)*d(3)-d(2)*d(3)))
       write(iounit,352)del
       write(iounit,*) ' '
       write(iounit,*) 'EIGENVALUES OF Drr TENSOR (1/s)'
       write(iounit,340)runit*d(1),runit*d(2),runit*d(3)
       DR=runit*trace
       WRITE(iounit,182)DR
182    FORMAT('  1/3 Trace DRR(EXACT) =',E15.5,' 1/s',/)
350    format(2x,'1/3 Trace = ',e12.5)
352    format(2x,'Anisotropy=',e12.5)

       call jacobi(Dtt,3,3,d,v,nrot)
       write(iounit,*) 'EIGENVALUES OF Dtt TENSOR (1/A)'
       write(iounit,340)d(1),d(2),d(3)
       trace=(d(1)+d(2)+d(3))/3.
       write(iounit,350)trace
       del=DSQRT(DABS(d(1)*d(1)+d(2)*d(2)+d(3)*d(3)
     $ -d(1)*d(2)-d(1)*d(3)-d(2)*d(3)))
       write(iounit,352) del
       write(iounit,*) ' '
       write(iounit,*) 'EIGENVALUES OF Dtt TENSOR (cm^2/s)'
       write(iounit,340)tunit*d(1),tunit*d(2),tunit*d(3)
       DTRANS=tunit*trace
       WRITE(iounit,134)DTRANS
134    FORMAT('  1/3 Trace DTT(EXACT) =',E15.5,' cm^2/s',/)

       IF(IVISCO.EQ.0) GOTO 1000

C      PRINT THE INTRINSIC VISCOSITY

       IF(ICOV.EQ.1) THEN
       IF(MW.EQ.1) THEN
       WRITE(IOUNIT,*) "INTRINSIC VISCOSITY AT THE CENTER OF VISCOSITY"
       ELSE
       WRITE(IOUNIT,*) "VISCOSITY FACTOR AT THE CENTER OF VISCOSITY"
       ENDIF
       WRITE(IOUNIT,*) " Vector between Center of Viscosity/Centroid"
       WRITE(IOUNIT,980)D1,D2,D3
980    FORMAT(3X,'  Dx = ',E12.5,'  Dy = ',E12.5,'  Dz = ',E12.5)
       ELSE
       IF(MW.EQ.1) THEN
       WRITE(IOUNIT,*) "INTRINSIC VISCOSITY AT THE CENTROID"
       ELSE
       WRITE(IOUNIT,*) "VISCOSITY FACTOR AT THE CENTROID"
       ENDIF
       ENDIF

       IF(IAREA.EQ.1) THEN
       CHI = CHIPRIME/TRUEVOL
       WRITE(IOUNIT,990)CHI
990    FORMAT(/,'     CHI               ',E16.6)
       WRITE(IOUNIT,991)CHIPRIME
991    FORMAT('     CHI*VOLUME (A^3)  ',E16.6,/)
       IF(MW.EQ.1) THEN
       ETA = CHIPRIME*0.602214/WM
       WRITE(IOUNIT,993)ETA
993    FORMAT('     ETA (cm^3/g)      ',E16.6,/)
       ENDIF
       ELSE
       WRITE(IOUNIT,997)CHIPRIME
997    FORMAT(/,'     CHI*VOLUME (A^3)  ',E16.6,/)
       IF(MW.EQ.1) THEN
       ETA = CHIPRIME*0.602214/WM
       WRITE(IOUNIT,999)ETA
999    FORMAT('     ETA (cm^3/g)      ',E16.6,/)
       ENDIF
       ENDIF

1000   CLOSE(9)
 
       RETURN
       END
