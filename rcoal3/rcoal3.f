	program rcoal
c This program uses the routine Coalesce (see below) to produce a set
c of input files for the program BEST.  The purpose of using rcoal 
c instead of Coalesce is to specify the number of triangles one desires
c in the BEST input files instead of using "edgfrac", the control para-
c meter in Coalesce.  In order to relate the number of triangles and
c e (edgfrac), we have determined that the functional form of this 
c relation can be approximately expressed by
c	N(e) = Nc /(1 + a e^p)
c where Nc is the number of triangles in the Connolly .c3p file.
c This program runs two standard calculations with e = .7 and e = 1.0
c in order to determine the values of the constants a and p.  Then
c  these values are used to determine the e value for a requested
c number of triangles. No output is written to disk for the initial two calcs. 

c Performance:  The relationship typically underestimates the number of
c triangles requested by 5%.  Thus, it is recommended that when setting
c the values of the triangle range you desire, that you increase those
c by 5% at the outset.  For example, if you desire 5 files with triangle 
c numbers varying between 6000 and 3000, then run rcoal as follows:
c   rcoal -f test.c3p -nmax 6300 -nmin 3150 -n 5

c To request a single file with a given number ntriang of triangles, set
c nmax = nmin = ntriang, n =1.

c To use the control parameter directly, specify only option "-e value" on the command line.
c In this case, no apriori estimate of the number of triangles will be made.
c This is useful when analyzing a molecular dynamics trajectory where one wants
c to use the same number of triangles for each run.  In that case we run rcoal
c once to get the edgefrac value that corresponds to the desired number of
c triangles, typically about 3500, and then call rcoal with that value of 
c edgefrac.  A typical command in this case is:
c	rcoal -f test.c3p -e 0.82
c Note that values of any other parameters are ignored if the -e option is used.

c Aragon 10/15/07

c  Addition of msms support.  Aragon 3/9/2012
c  For Sanner msms files, the new flag "-s" is included in the command line.
c  The input file given is the molecule.face from the msms output.  The program will
c  automatically also use the molecule.vert file produced.
c  To request various files within a triangle range:
c   rcoal -f molecule.face -nmax 6300 -nmin 3150 -n 4 -s

c  To request a triangulation on the basis of "edgefrac":
c	rcoal -f molecule.face -e 0.82 -s


	dimension nt(10),nn(10),e(10)
c     Set by Coalesce: the number of triangles of zero area in its last
c     triangulation, and the files written so far (two for each).
      common /zarea/ nzero
      common /outlst/ nout
      character*80 outf
      common /outnam/ outf(20)
	CHARACTER*80 tfile, text
	logical sanner, connolly 

        nmax = -1
        nmin = -1
        n = -1
	edgefrac = -1.0

	print*,"Rcoal Ver. 3.0"
      nzero = 0
      nout = 0
	print*, "SFSU, Aragon Lab."

C READ INPUT FROM COMMAND LINE
        CALL INPUT(tfile,NMAX,NMIN,N,WM,edgefrac,sanner)

c  Set Connolly or Sanner case: Connolly is default.
	connolly = .true.
	if (sanner) then
	connolly = .false.
	endif

c Read number of triangles from input file
	connolly1_case: if (connolly) then
	open(10,file=tfile,status='old',iostat=ierr)
	read(10,*)dum,dum,ntriang
	close(10)
	endif connolly1_case

	sanner1_case: if (sanner) then
c  Make sure tfile is of the form molecule.face

	do j=1,28
        k=j
        if(k+3.gt.28) then
        print *, 'Input file is not a Sanner .face file '
        stop
        endif
        if(tfile(k:k+3).eq.'.fac') exit
	enddo

	open(10,file=tfile,status='old',iostat=ierr)
	read(10,*) text
	read(10,*) text
	read(10,*)ntriang, dum,dum,dum
	close(10)
	endif sanner1_case

c If value of edgefrac was not specified upon input, then we assume that
c the user wants files with specific number of triangles.
	if(edgefrac.gt.0.0) goto 110

        if(nmax.eq.-1.OR.nmin.le.3.OR.n.eq.0.OR.n.eq.-1) then
	print*,"Check input parameters: n=0?, nmin <3,n,nmax missing?"
        goto 220
        else if(nmax.lt.nmin) then
	print*,"nmax is less than nmin"
        goto 220
        endif

        if(ntriang.lt.nmax) then
        print*,'Nmax must be less than',ntriang+1,'in this instance.'
        goto 220
        endif  

c If value of edgefrac was not specified upon input, then we assume that
c the user wants files with specific number of triangles.
	if(edgefrac.gt.0.0) goto 110
	if(nmax.eq.nmin)then
	 n = 1
	goto 100
	endif
	if(n.gt.1) ndelta = (nmax-nmin)/(n-1)

100     print *, ' '
	print *, 'Target values of n'
	do j=1,n
	nt(j) = nmin + (j-1)*ndelta
 	print *, '             ', nt(j)
	end do
        print *, ' '

c Make two standard runs of coalesce to collect data of N(e).
	call scoalesce(tfile,0.7,n1,-1,sanner)
	print *, "Just finished first scoalesce call"
	call scoalesce(tfile,1.0,n2,-1,sanner)
	print *, 'initial n ',n1,n2

c Generate two parameters of the formula N(e) = N/(1+ a*e^p): a, p
	y2 = log(float(ntriang)/float(n2) -1.0)
	y1 = log(float(ntriang)/float(n1) -1.0)
	x2 = log(1.0)
	x1 = log(0.7)
	p = (y2-y1)/(x2-x1)
	a = exp(y2-p*x2)
	print *, ' p = ', p
	print *, ' a = ', a

c Generate list of edgfrac values from formula & call scoalesce.
	do j=1,n
	e(j) = 1.000*((float(ntriang)/float(nt(j)) -1.0)/a)**(1/p)
	print *, ' e = ',e(j)
	call scoalesce(tfile,e(j),nn(j),1,sanner)
      if(nzero.gt.0) goto 230
	print *, 'Target n         ',nt(j)
        print *, 'Actual coalesce n',nn(j)
        print *, ' '
	end do
	goto 220
110	continue
	call scoalesce(tfile,edgefrac,n3,1,sanner)
      if(nzero.gt.0) goto 230
        print *, 'Actual coalesce n',n3
c     List the files written, as Coalesce did after each triangulation.
220   do i=1,nout-1,2
      print *,' '
      print *,' Output files written:'
      print *,'                      ',outf(i)
      print *,'                      ',outf(i+1)
      print *,' '
      end do
      stop
c     A triangulation has a triangle with two vertices at the same
c     position (msroll can give them where its surface folds onto
c     itself), which BEST cannot use.  Remove the files already
c     written and fail, so that the caller can try another input.
230   do i=1,nout
      open(16,file=outf(i),status='old',iostat=ierr)
      if(ierr.eq.0) close(16,status='delete')
      end do
      write(0,*) 'rcoal: a triangulation has triangles of zero area'
      write(0,*) '(two vertices at the same position), which BEST'
      write(0,*) 'cannot use.  No triangulation written;'
      write(0,*) 'try a coarser input triangulation.'
      flush(0)
      stop 1
	end program

       SUBROUTINE INPUT(tfile,NMAX,NMIN,N,WM,edgefrac,sanner)

       CHARACTER*80 tfile ,ARG
	logical sanner

	sanner = .false.

       NARG = IARGC()
       IF(NARG.LT.2) then
      PRINT*,'Usage: rcoal -f filename -nmax val -nmin val '
      PRINT*,'-n numb [-s]'
      PRINT*,'Usage: rcoal -f filename -e val [-s]'
	STOP
       ENDIF

c By default, rcoal assumes the file is a Connolly .c3p file.  If the flag -s
c is present, then rcoal will pass msms sanner = .true.  to coalesce.

       DO K=1,NARG
       CALL GETARG(K,ARG)
       IF (ARG(1:5).EQ.'-f') THEN
       CALL GETARG(K+1,ARG)
       READ(ARG,*) tfile 
       END IF
       IF (ARG(1:5).EQ.'-nmax') THEN
       CALL GETARG(K+1,ARG)
       READ(ARG,*) NMAX
       ENDIF
       IF (ARG(1:5).EQ.'-nmin') THEN
       CALL GETARG(K+1,ARG)
       READ(ARG,*) NMIN
       ENDIF
       IF (ARG(1:2).EQ.'-n') THEN
       CALL GETARG(K+1,ARG)
       READ(ARG,*) N
       ENDIF
       IF (ARG(1:3).EQ.'-mw') THEN                             
       CALL GETARG(K+1,ARG)                                   
       READ(ARG,*) WM
       ENDIF                                                  
       IF (ARG(1:2).EQ.'-e') THEN                             
       CALL GETARG(K+1,ARG)                                   
       READ(ARG,*) edgefrac 
       ENDIF                                                  
       IF (ARG(1:2).EQ.'-s') THEN                             
	sanner = .true.
       ENDIF                                                  
       ENDDO

       RETURN
       END

