       subroutine input(infnm,iarea,truearea,TRUEVOL,cg,
     $ IVISCO,ICOV,MW,WM,IHM,HEPS,IHLU)

c      This subroutine reads a triangulation input file.  Only the number of triangles and the vertex 
c      coordinates of each triangle are used as input by the BEST program. The user should change this 
c      subroutine to take the format of the input data into account.  Retain the /com1/ common statement. 
c      S. Aragon  Sum 2001

c      Subroutine arguments:
c      ver = version number of BEST program.  Assigned in MAIN. Input.
c      outfnm = output filename. Output
c      iounit = Fortran unit number for writing the output. Output.
c      nar = number of triangles. Output.
c      iarea = flag to control area correction.  iarea= 1 implies do area correction.
c      truearea = analytically computed area of molecule.  It is larger than the polyhedron area.
c      if iarea = 1, then a value of truearea is requested at input.
c      The vertex coordinate arrays are filled in common: v1,v2,v3.
c      The incenter array, x, and the area array are calculated elsewhere.
c      ihm = solver: 1 H-matrix (-hm), -1 dense (-dn), 0 chosen by the number of triangles.
c      heps = H-matrix accuracy (-he, default 1e-4); ihlu = 1: direct H-LU instead of GMRES (-hlu).
c      E. H. Brookes 11/13/2026: -hm, -dn, -he, -hlu.

       character*80 inline,infnm,outfnm2,infnm3,arg

       NARG=iargc()

       if(narg.lt.2) then
       print *,
     $ 'Usage : best -f input_filename [-na] [-p] [-v] [-vc] [-mw]',
     $ ' [-hm | -dn] [-he eps] [-hlu]'
       stop
       endif

       IAREA = 1
       CG = 1.0
       IVISCO = 0
       ICOV = 0
       MW = 0
       IHM = 0
       HEPS = 1.0E-4
       IHLU = 0
       DO K = 1,NARG
       CALL GETARG(K,ARG)
       IF(ARG(1:2).EQ.'-f') THEN
       CALL GETARG(K+1,ARG)
       INFNM = ARG
       IFILE = 1
       ENDIF
       IF (ARG(1:3).EQ.'-na') IAREA = 0
       IF (ARG(1:2).EQ.'-p') CG = 0.0
       IF (ARG(1:2).EQ.'-v') IVISCO = 1
       IF (ARG(1:3).EQ.'-vc') ICOV = 1
       IF (ARG(1:3).EQ.'-mw') THEN
       CALL GETARG(K+1,ARG)
       READ(ARG,*) WM
       MW = 1
       ENDIF
       IF (ARG(1:3).EQ.'-hm') IHM = 1
       IF (ARG(1:3).EQ.'-dn') IHM = -1
       IF (ARG(1:4).EQ.'-hlu') IHLU = 1
       IF (ARG(1:3).EQ.'-he') THEN
       CALL GETARG(K+1,ARG)
       READ(ARG,*) HEPS
       ENDIF
       ENDDO

       return
       end
