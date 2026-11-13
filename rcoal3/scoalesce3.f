C     PROGRAM coalesce.f  This is the subroutine version scoalesce.f.

c	AUTHOR:  Sergio Aragon
c			Professor of Chemistry
c			San Francisco State University
c			aragons@sfsu.edu

c	DESCRIPTION

C	COALESCE modifies triangulated surfaces to contain only triangles 
c	suitable for use in Boundary Element Calculations. COALESCE can also be
c	used to coarsen a triangulation by reducing the triangle number. 
c	Triangulations given by Connolly's MSROLL and Sanner's MSMS are 
c	supported. (Zauhar's SMART is in progress).  Coding is in Fortran 90.

C	PROBLEMS THAT COALESCE FIXES
c	Triangulation programs, such as the above, produce copious amounts of
c	slender triangles of small areas that can cause problems in BE calculations.
c	These slender triangles are of two types:
c		Slender base triangles that have one very acute angle.
c		Slender obtuse triangles that have two very acute angles.
c	Infrequently, one may also find Degenerate triangles (a line) which have zero area.
c     Coalesce removes degenerate, slender base, and slender obtuse triangles from
c	the triangulation.
c	If a triangulation contains too many triangles, coalesce can be used to
c	reduce the number by coalescing smaller triangles into larger ones.
c	COALESCE tries to preserve the shape as much as possible during this process.
 
c	METHODS
c	First, slender base triangles are removed. The average edgelength is computed.
c	The control variable is the fraction, edgfrac, of the average 
c	edge length that is considered to be small. Triangles with an edge smaller
c	than edgfrac*average_length are removed. The other triangles sharing the vertices
c	of triangles to be removed, are moved to a weighted intermediate point between the
c	original vertices. The default value is edgfrac =0.11.
c	Because coalescing the vertices of a small edge can create new triangles which
c	are small, a second pass of coalescence is done. Typically, 90% of small
c	triangles are eliminated on the first pass, and the rest in the second pass.
c	If a value greater than 1 is given to edgefrac, COALESCE will eventually
c	yield a triangulation with few triangles irrespective of the value of edgfrac.
c	Some triangles that are smaller than the cutoff size will survive even the
c	second pass and these will never be eliminated.  Thus a minimum number of
c	triangles will always be produced, of order 100 when the starting value is
c	several thousand.

c	Second, after slender base elimination, the slender obtuse triangles are
c	identified.  There are many fewer of this type in typical files. The slender
c	obtuse triangle and the triangle joined to it opposite the obtuse vertex
c	are repartioned across their opposing vertices, generating two new triangles.
c	The triangle count stays the same in this process. Infrequently, the new
c	triangles produced may have small areas.  These are checked and if necessary,
c	a coalescence is done to eliminate them.
c
c	Finally, upon user request, triangles that have one or two vertices not shared
c	with any other triangle ("dangling triangles") are eliminated.  Selecting
c	the -d option accomplishes this.  Such triangles can ocurr often for PDB files
c	containing waters of hydration.  If these are not removed, from the input
c	file, they may result in a dangling triangle.  The hydrodynamic interaction
c	of will affect the overall calculation somewhat, depending on the size of
c	the molecule.  For a large protein, the effect is minimal and it is best to
c	delete them to decrease the number of triangles.  Dangling triangles may
c	also appear when a pointy feature exists in a molecule and the triangulation
c	is coarsened greatly. In this case, the dangling triangle should not be
c	removed.

c	OUTPUT
c	The program produces two output files, named after the input file [input].xyz,with 
c	the 'xyz' ending striped off. Processing statistics are provided on the
c	console as well as in the output histogram file.
 
c	[input].his contains a histogram distributions around 0-.5 and .5-4.0 A**2
c        of the patch areas.
c	[input] contains the vertex coordinates of the vertices of each triangle
c	  see format below. This format is used in the BEST hydrodynamics program.

c	USAGE : typing "coalesce" without arguments produces the usage statement:
c	'Usage : coalesce -c connolly_filename.c3p '
c	'         -s sanner_filename.vert sanner_filename.face '
c	'			-z smart_filename.sm '
c	'         [-e edgfrac]'
c	'            [-d] [-v] [-o ncom]'
c	 The -c, -s, -z options select the type of input file as indicated.
c	 The -e option preceeds the value of edgfrac chosen.
c	 The -d option requests that dangling triangles be eliminated.
c	 The -v option provides verbose information.
c TR   The -o option provides for choosing a Connolly component surface. The outer surface
c	  is usually component 1 (ncom =1, the default). See Connolly for details.
c           
c	HISTORY
c	8/2/04 Sergio Aragon
c	To preserve shape, added instructions to move both vertices to a position in between
c	the original vertices, weighted by the number of triangles on each vertex.
c     8/25/01 Sergio Aragon
c	Slender base triangle deletion and coalescence algorithm developed.
C	Connolly's msroll output (*.c3p) implemented. Version 1.0
c	9/15/01 Sergio Aragon
c	Program modified to read Sanner input files as an option. The
c	distinct edge list is computed from the Sanner molecule.face file.
c	11/30/01 Aragon
c	Obtuse slender triangles are found and reshaped, after the slender base
c	triangles are removed. This removes small area triangles without requiring
c	a large value of edgfrac! Version 2.0
c	4/17/02 Aragon
c	 Program logic simplified in slender_del routine. Added option to remove
c	dangling triangles. Removed count that required that triangles be removed
c	in pairs. Added verbose option. Version 2.5
c       4/22/02 T. Rosales
c       Added -o option to choose a Connolly component surface. Default is component 1
c		Version 2.6

c  Operation of the subroutine version:
c	The subroutine takes a triangulation file, tfile, as input.  It can be either a 
c	Connolly .c3p file or a Sanner .face file.
c	By default, scoalesce expects a .c3p file. Otherwise, the sanner logical variable
c	must be set to .true. .  In that case, the program will check that tfile is of the form
c	basename.face .  The program will assume that if this is true, then a file basename.vert
c	also exists and will read it.  The program will stop with an error is this is not the case.
c	The subroutine will output the number of triangles it produces in the variable npoly.
c	The subroutine requires the control variable edgefrac to coarsen the input triangulation.
c	edgefrac is described above.  
c	The np=-1 control variable causes no output to be printed.  This is useful when scoalesce
c	is called from rcoal the first couple of times to determine the parameters to produce
c	the number of triangles desired. If np is not -1, output will be produced.

	subroutine SCOALESCE(tfile,edgfrac,npoly,np,sanner)
	parameter(ver=3.5)
	integer,dimension(:,:),allocatable :: vn,edgn,sedg
      real,dimension(:,:),allocatable :: v1,v2,v3,x
      real,dimension(:),allocatable :: area,edglen
	integer, dimension(:), allocatable :: indx,nsut,nsutc,iobv

	dimension::bins(50),bins2(50),aux(3)
	dimension::ndist1(50),ndist2(50),oarea(200)
	dimension::vx(3),vy(3),vz(3)
      integer first,second
	logical connolly, sanner, smart,verb,dangle
      logical lflip
      common /lnkcnt/ nlinkskip
      common /zarea/ nzero
      common /outlst/ nout
      character*80 outf
      common /outnam/ outf(20)
      character*51 cv1,cv2,cv3
       CHARACTER*80 tfile,CONC3V,OUTPUT,OUTPUTHIS, basename,
     $ smar,CON3b,inputf, outputf, vfile, text1,text2

	connolly = .true.
      nlinkskip = 0
      nflipskip = 0
	smart = .false.
	verb = .false.
	dangle = .false.
	ncom=1
	if (sanner) then
        connolly = .false.
        endif

	print*, " "
	print*, "Begin Coalesce Processing"
	print*, " "
 
c	Determine the basename for histogram and output files.
      do  i=1,28 
		L=i                                                           
      if (tfile(L:L+3).eq.'.c3p'.or. tfile(L:L+3).eq.'.fac') goto 20 
	basename = tfile(1:L)
	end do
20	continue
	print *, "basename = ", basename
c  Note that due to the above loop structure, the length of basename is L-1
c  upon exiting the loop. 
C	READ INPUT FILES

c	CONNOLLY CASE
c	Read data from input file.
	connolly_case: if(connolly) then
	open(10,file=tfile,status='old',iostat=ierr)

c	Read number of vertices, edges, and triangles
	read(10,*)nvert,nedge,ntriang
       allocate(x(nvert,3),nsut(ntriang),area(ntriang))
c	Read vertex coordinate records 
	
        do n=1,nvert 
         read(10,*)x(n,1),x(n,2),x(n,3),dum,dum,
     $dum,dum,dum,dum,dum,dum,dum
	end do
	

c     Read edge records
      allocate(edgn(nedge,2),vn(ntriang,3))
        do n=1,nedge 
	    read(10,*)edgn(n,1),edgn(n,2),dum,dum,dum
	  end do

c	Read triangle records & compute area.
 	totarea = 0.0
       do n=1,ntriang
       read(10,*)dum,dum,dum,vn(n,1),vn(n,2),vn(n,3),nsut(n),dum,dum

c	Calculate the total area of the triangles.
		do k=1,3 
		vx(k) = x(vn(n,1),k)
		vy(k) = x(vn(n,2),k)
		vz(k) = x(vn(n,3),k)
		end do
	call areas(vx,vy,vz,area(n))
	totarea = totarea+area(n)
	end do
	beginarea = totarea
	deallocate(area)
	close(10)

	conc3v = tfile(1:L-1)//'.c3v'
       open(14,file=conc3v)
       do n=1,12
       read(14,'(a80)')inline
       end do
       read(14,*)i1,x2,x3,x4,TRUEVOL,truearea,x5
       close(unit=14)       

	end if connolly_case

c	SMART CASE
c	Read data from input file.
	smart_case: if(smart) then
	tfile = smar
	open(10,file=smar,status='old',iostat=ierr)

c	Read number of vertices, edges, and triangles
	read(10,*)nvert,nedge,ntriang
       allocate(x(nvert,3),nsut(ntriang))
c	Read vertex coordinate records
        do n=1,nvert 
         read(10,*)x(n,1),x(n,2),x(n,3),dum,dum,
     $dum,dum,dum,dum,dum,dum,dum
	end do

c     Read edge records
      allocate(edgn(nedge,2),vn(ntriang,3))
        do n=1,nedge 
	    read(10,*)edgn(n,1),edgn(n,2),dum,dum,dum
	  end do

c	Read triangle records. 
       do n=1,ntriang
      read(10,*)dum,dum,dum,vn(n,1),vn(n,2),vn(n,3),dum,dum,dum
         nsut(n)=1
        end do

	end if smart_case

c	SANNER CASE.
	sanner_case: if(sanner) then

c	First make sure tfile is of the form basename.face . 
	
	if(tfile(L:L+3).ne.'.fac') then
	print *, "tfile piece = ", tfile(L:L+3)
	print *, 'No valid Sanner input files'
	stop
	endif
	vfile = basename(1:L-1)//'.vert'
c	 print *, "vertex file = ", vfile
	nvunit =11
	nfunit = 10
	open(nfunit,file=tfile,status='old',iostat=ierr)
	open(nvunit,file=vfile,status='old',iostat=ierr)
c	Read number of vertices and triangles. Skip Comment lines.
	read(nvunit,*)
	read(nvunit,*)
	read(nvunit,*)nvert,ndum,dum,dum
       allocate(x(nvert,3))
c	Read vertex coordinate records
      do n=1,nvert 
      read(nvunit,*)x(n,1),x(n,2),x(n,3),dum,dum,dum,ndum,ndum
	end do
	close(nvunit)
c	print *, "Just read vertices"
c	Read triangle records. Skip comment lines.
	read(nfunit,*)
	read(nfunit,*)
 	read(nfunit,*)ntriang,ndum,dum,dum
	print *, "Number triangles = ", ntriang
       allocate(vn(ntriang,3),nsut(ntriang),area(ntriang))
	totarea = 0.0
       do n=1,ntriang
       read(nfunit,*)vn(n,1),vn(n,2),vn(n,3),ndum,ndum
       nsut(n)=1
c	Calculate the total area of the triangles.
		do k=1,3 
		vx(k) = x(vn(n,1),k)
		vy(k) = x(vn(n,2),k)
		vz(k) = x(vn(n,3),k)
		end do
	call areas(vx,vy,vz,area(n))
	totarea = totarea+area(n)
	enddo
	beginarea = totarea
	deallocate(area)
	close(nfunit)
c	Generate the edge records from the triangle records. Sanner case.
c	Collect all edges, 3 from each triangle, with many repetitions.
	allocate(sedg(3*ntriang,2))
	k = 1
	do n =1, ntriang

	sedg(k,1) = vn(n,1)
	sedg(k,2) = vn(n,2)
	sedg(k+1,1) = vn(n,2)
	sedg(k+1,2) = vn(n,3)
	sedg(k+2,1) = vn(n,3)
	sedg(k+2,2) = vn(n,1)
	k = k+3
	end do

	if(k-1.ne.3*ntriang) then
	print *,'Edge count mismatch'
	print *,'total # edges= ',k-1,' and ',3*ntriang
	else
c	print *,'Edge count correct= ',k-1
	continue
	end if

c  Now mark edges that are repeated by setting the vertex number to -1.
	do n = 1, 3*ntriang

	if(sedg(n,1).ne.-1) then
	first = sedg(n,1)
	second = sedg(n,2)
		do k=1, 3*ntriang
		if(k.ne.n) then

		if(sedg(k,1).eq.first.and.sedg(k,2).eq.second) then
		sedg(k,1) = -1
		sedg(k,2) = -1
		else if(sedg(k,2).eq.first.and.sedg(k,1).eq.second) then
		sedg(k,1) = -1
		sedg(k,2) = -1
		end if
		endif
		end do
	end if
	end do

c	Now assemble all unique edges
	allocate(edgn(3*ntriang,2))
	nj=0.
	do n = 1, 3*ntriang
	if(sedg(n,1).ne.-1) then
	nj=nj+1
	edgn(nj,1) = sedg(n,1)
	edgn(nj,2) = sedg(n,2)
c	print *,'Edge= ',nj,' vertices= ',edgn(nj,1),edgn(nj,2)
c	if(modulo(nj,20).eq.1) pause
	end if
	end do
	deallocate(sedg)
	nedge = nj
	print *, nj, ' unique edges found'

c	Read msms.log file to obtain areas and volume.

	inputf = basename(1:L-1)//'.msms.log'
	print *, "Input file = ", inputf
	open(12,file=inputf,status='old',iostat=ierr)
	 do j=1,200
        read(12,*,iostat=ierr) text1, text2
        if (ierr.lt.0) exit
	if (text1(1:10).eq.'ANALYTICAL'.and.text2(1:7).eq.'SURFACE') then
                read(12,*)
                read(12,*) dum,dum,dum,dum,dum,truearea,dum
        endif
        if (text1(1:9).eq.'NUMERICAL') then
                read(12,*)
                read(12,*) dum,dum, polyvol, polyarea
        endif
        enddo
	close(12)

c	MSMS does not compute an analytical volume.  Approximate it from
c	the corresponding ratio of the areas.
c	truevol = polyvol*(truearea/polyarea)**1.5

c	Potential problem:  msms produces a polygonal volume that is larger
c	than the analytical volume produced by msroll.  Furtheremore,
c	msms produces an analytical area that is slightly larger than that
c	of msroll (for 1crn.pdb, for example), but a much larger polyvol
c	than produced by msroll.  Must test this to see who is right.
c	Meanwhile, we just set truevol = polyvol for the sanner case.
	truevol = polyvol
	print *, "Areas: True, Poly = ", truearea, polyarea
	print *, "Volumes: True, Poly =", truevol, polyvol
	end if sanner_case

c	SMALL TRIANGLE ELIMINATION

C	The basic information in building a valid polygon is the presence of
c	a valid vertex number record in the array vn.

c	DEGENERATE triangle elimination
c	Degenerate triangles have zero area, which implies at least two
c	identical vertex numbers. We flag the bad triangle by setting its vertices to -1.
	 idegen = 0
	do m=1,ntriang
	if(vn(m,1).eq.vn(m,2).or.vn(m,1).eq.vn(m,3).or.
     $   vn(m,2).eq.vn(m,3) .or. nsut(m) .ne. ncom) then
		idegen = idegen +1
		vn(m,1) = -1
		vn(m,2) = -1
		vn(m,3) = -1
	end if
        end do

c	SMALL EDGE detection.  Compute the edge lengths.
	ipass = 1
100	continue
	allocate(edglen(nedge))
	avelength = 0.0
	do j = 1, nedge
	aux(1) = x(edgn(j,1),1)-x(edgn(j,2),1)
	aux(2) = x(edgn(j,1),2)-x(edgn(j,2),2)
	aux(3) = x(edgn(j,1),3)-x(edgn(j,2),3)
	edglen(j) = sqrt(dot(aux,aux))
	avelength = avelength +edglen(j)
	end do
	cutofflength = edgfrac*avelength/nedge

c	Construct an index array to indicate the order of the edge lengths.
	allocate(indx(nedge))
	call index(nedge,edglen,indx)

c	Now search the edge list in order of smallest to largest.
	idelete = 0

c	SLENDER BASE TRIANGLE ELIMINATION
	do j = 1, nedge

	if(edglen(indx(j)).ge.cutofflength) exit
c	A small edge has been found.  Now find the triangles that share
c	 those edge vertices, delete the slender base triangles and coalesce.
	
	call slender_del(edgn(indx(j),1),edgn(indx(j),2),x,nvert,vn,
	1 ntriang,idelete)
	
	end do 

 	print *, 'Finished deleting small triangles: Pass= ',ipass
	print *," "
	deallocate(edgn,edglen,indx)


	if(ipass.eq.2) then
	idang = 0
c	DELETE DANGLING TRIANGLES IF REQUESTED
	if(dangle) call dangle_del(ntriang,vn,verb,idang)
	endif
		
	if(ipass.eq.1) npolygon = ntriang-idegen-idelete

c	RESHAPING SLENDER OBTUSE TRIANGLES
	if(ipass.eq.1) goto 125
	slender = 0.97
	slenderfrac = slender/(1+slender)
	nobtuse = 0
c     One entry per triangle at most, so the list cannot overflow.
	allocate(indx(ntriang),iobv(ntriang))
	slenderloop: do j=1, ntriang
	if(vn(j,1).ne.-1) then
c	Compute the edge lengths
	do n =1,3
	vx(n) = x(vn(j,1),n)
	vy(n) = x(vn(j,2),n)
	vz(n) = x(vn(j,3),n)
	enddo
	call slength(vy,vx,slen21)
	call slength(vz,vx,slen31)
	call slength(vz,vy,slen32)
c	Find the largest edge length and test for slenderness.

	slenmax = max(slen21,slen31,slen32)
	perim = slen21+slen31+slen32
	if(slenmax.ge.slenderfrac*perim) then
c	Triangle is either slender obtuse or slender base
c	A slender obtuse triangle has two very acute angles.
c	A slender base has only one very acute angle.
c	If triangle is slender base, then it has a small edge.
c	The small edges have already been eliminated.
c	Test for slender obtuse: min/max greater than 1/3.
	slenmin = min(slen21,slen31,slen32)

		if(slenmin/slenmax.ge.0.33) then
c	A slender obtuse triangle has been found. Store the triangle index.
	nobtuse = nobtuse +1
	indx(nobtuse) = j
c	Determine which vertex is the obtuse one. Store it in iobv.
c     Break ties in a fixed order, so that iobv is always set.
      if(slen32.ge.slen21.and.slen32.ge.slen31) then
      iobv(nobtuse) = 1
      else if(slen31.ge.slen21.and.slen31.ge.slen32) then
      iobv(nobtuse) = 2
      else
      iobv(nobtuse) = 3
      endif

		endif
	endif
	endif
	enddo slenderloop
	print *, nobtuse, " Slender obtuse triangles found!"

c	Now, reshape the slender obtuse triangles.
c	Start by finding the triangle opposite to the obtuse vertex.
	icon = 0
	do j = 1, nobtuse
	L1= iobv(j)-1
	L2= iobv(j)+1
	if(L1.eq.0) L1 = 3
	if(L2.eq.4) L2 = 1
	iflag = 0
	icon = icon + 1
	do k = 1 ,ntriang
	if(k.ne.indx(j).and.vn(k,1).ne.-1.and.iflag.eq.0) then
c	We simply need to find the other triangle that shares the acute
c	vertices with the obtuse triangle. iflag ensures we use this loop only once per j.

	if(vn(indx(j),L1).eq.vn(k,1).or.vn(indx(j),L1).eq.vn(k,2).
     $  or.vn(indx(j),L1).eq.vn(k,3)) then
c This triangle contains one acute vertex.
		if(vn(indx(j),L2).eq.vn(k,1).or.vn(indx(j),L2).eq.vn(k,2).
     $	 or.vn(indx(j),L2).eq.vn(k,3)) then
c  This triangle also contains the second acute vertex. They share a side.
c  Now make two new triangles: the obtuse vertex is given to the opposite triangle.
c  One acute vertex is replaced by the middle vertex of the other triangle.	
	
c     Skip the flip if its new edge, from the obtuse vertex to the other
c     triangle's third vertex, already exists: flipping would make that
c     edge non-manifold.  Skipped flips are counted in nflipskip.
      ia = vn(indx(j),iobv(j))
      id = 0
      do m = 1,3
      if(vn(k,m).ne.vn(indx(j),L1).and.vn(k,m).ne.vn(indx(j),L2))
     $ id = vn(k,m)
      enddo
      lflip = ia.ne.id
      do kk = 1,ntriang
      if(lflip.and.kk.ne.k.and.kk.ne.indx(j).and.vn(kk,1).ne.-1) then
      if((vn(kk,1).eq.ia.or.vn(kk,2).eq.ia.or.vn(kk,3).eq.ia).and.
     $   (vn(kk,1).eq.id.or.vn(kk,2).eq.id.or.vn(kk,3).eq.id))
     $ lflip = .false.
      endif
      enddo
      if(.not.lflip) nflipskip = nflipskip+1
      if(lflip) then
	if(vn(indx(j),L1).eq.vn(k,1).and.vn(indx(j),L2).eq.vn(k,2)) then 
			vn(indx(j),L1) = vn(k,3)
			vn(k,2) = vn(indx(j),iobv(j))
	else if(vn(indx(j),L1).eq.vn(k,1).and.vn(indx(j),L2).eq.vn(k,3)) 
	1then 
			vn(indx(j),L1) = vn(k,2)
			vn(k,3) = vn(indx(j),iobv(j))
							
	else if(vn(indx(j),L1).eq.vn(k,2).and.vn(indx(j),L2).eq.vn(k,3))
	1 then 
			vn(indx(j),L1) = vn(k,1)
			vn(k,3) = vn(indx(j),iobv(j))
	else if(vn(indx(j),L1).eq.vn(k,2).and.vn(indx(j),L2).eq.vn(k,1))
	1 then 
			vn(indx(j),L1) = vn(k,3)
			vn(k,1) = vn(indx(j),iobv(j))

	else if(vn(indx(j),L1).eq.vn(k,3).and.vn(indx(j),L2).eq.vn(k,1))
	1 then 
			vn(indx(j),L1) = vn(k,2)
			vn(k,1) = vn(indx(j),iobv(j))
	else if(vn(indx(j),L1).eq.vn(k,3).and.vn(indx(j),L2).eq.vn(k,2))
	1 then 
			vn(indx(j),L1) = vn(k,1)
			vn(k,2) = vn(indx(j),iobv(j))
	end if
c	print *,indx(j),' New ob vrtx',vn(indx(j),1),vn(indx(j),2),
c	1vn(indx(j),3)
c	print *,k,' New oppst vrtx',vn(k,1),vn(k,2),vn(k,3)
	
c	We have now replaced the slender obtuse triangle and its opposite triangle,
c	by two new triangles with decent areas.  However, if the opposite triangle
c	was also slender obtuse, the two new triangles MAY  have small
c	areas, and will share a new small edge.  This must be fixed.
c	Test that the new areas are ok:
         do n=1,3
            vx(n)=x(vn(k,1),n)
            vy(n)=x(vn(k,2),n)
            vz(n)=x(vn(k,3),n)
         end do
		call areas(vx,vy,vz,areanew)
c	print *,'obt=',vn(indx(j),L1),vn(indx(j),L2),vn(indx(j),
c	1 iobv(j)),' area=',areanew
		if(areanew.lt.0.06) then
		if(verb) print *,'Opposite:new area too small, will coalesce'

		
c     Pass copies: slender_del writes -1 into vn, so passing elements of
c     vn would let it overwrite its own arguments.
      iv1 = vn(k,1)
      iv2 = vn(k,2)
      call slender_del(iv1,iv2,x,nvert,vn,ntriang,idelete)
		endif
c     The coalescence above may have deleted this triangle too.
      if(vn(indx(j),1).ne.-1) then
         do n=1,3
            vx(n)=x(vn(indx(j),1),n)
            vy(n)=x(vn(indx(j),2),n)
            vz(n)=x(vn(indx(j),3),n)
         end do
		call areas(vx,vy,vz,areanew)
c	print *,'match:',vn(k,1),vn(k,2),vn(k,3),' area=',areanew
		if(areanew.lt.0.06) then
		if (verb) print *,'Obtuse: new area too small, will coalesce'
	
c     Pass copies, as above.
      iv1 = vn(indx(j),L1)
      iv2 = vn(indx(j),L2)
      call slender_del(iv1,iv2,x,nvert,vn,ntriang,idelete)
		endif
      endif
      endif
	iflag = 1
		endif
		
	endif
	endif
	enddo

	enddo
	npolygon = npolygon-idang-idelete
	print *,'Processed ',icon,' obtuse slender triangles.'

125	continue

C	BUILD NEW TRIANGLES

c     Read the vertex indices and get the coordinates of each vertex.
c     The vertex coordinates are in arrays v1,v2,v3. The array index
c	icount should reach up to npolygon.

	if(ipass.eq.1) then
	allocate(v1(npolygon,3),v2(npolygon,3),v3(npolygon,3))
	allocate(area(npolygon))
	end if

	totarea=0.0
	icount = 0
         do k=1,ntriang
	if(vn(k,1).ne.-1) then
		icount= icount+1
         do n=1,3
            v1(icount,n)=x(vn(k,1),n)
            v2(icount,n)=x(vn(k,2),n)
            v3(icount,n)=x(vn(k,3),n)
		  vx(n) = v1(icount,n)
		  vy(n) = v2(icount,n)
		  vz(n) = v3(icount,n)
         end do
c	Calculate the total area of the triangles. 

	call areas(vx,vy,vz,area(icount))
	totarea = totarea+area(icount)

c	Print out characteristics of any remaining small triangles after pass 2.     	
		if(area(icount).lt.0.06.and.ipass.eq.2.and.verb) then
	call slength(vy,vx,slen21)
	call slength(vz,vx,slen31)
	call slength(vz,vy,slen32)
	print *,'Small Triangle; lengths: ', slen21,slen31,slen32
	print *,icount,' triangle ',vn(k,1),vn(k,2),vn(k,3),area(icount)
	if(modulo(icount,25).eq.1) pause
		end if
	end if
        end do
	npoly = icount
	
	if(npoly.ne.npolygon) print *,"Triangles counted in pass",ipass,":"
	1 ,npoly," computed npolygon= ",npolygon



	if(ipass.eq.1) then
	idelete1 = idelete
	npolygon1 = npolygon
	totarea1 = totarea
	end if

C	CHECK NEW EDGE LIST FOR NEW SMALL TRIANGLES
	if(ipass.eq.2) goto 200
c	This section is skipped on the second pass.
	ipass = 2
	print *,' '
	print *, '******Now checking new edge list*******'
	allocate(sedg(3*npolygon,2))
	k = 1
	do n =1, ntriang
	if(vn(n,1).ne.-1) then
	sedg(k,1) = vn(n,1)
	sedg(k,2) = vn(n,2)
	sedg(k+1,1) = vn(n,2)
	sedg(k+1,2) = vn(n,3)
	sedg(k+2,1) = vn(n,3)
	sedg(k+2,2) = vn(n,1)
	k = k+3
	end if
	end do

c  Now mark edges that are repeated by setting the vertex number to -1.
	do n = 1, 3*npolygon

	if(sedg(n,1).ne.-1) then

	first = sedg(n,1)
	second = sedg(n,2)
		do k=1, 3*npolygon
		if(k.ne.n) then
		if(sedg(k,1).eq.first.and.sedg(k,2).eq.second) then
		sedg(k,1) = -1
		sedg(k,2) = -1
		else if(sedg(k,2).eq.first.and.sedg(k,1).eq.second) then
		sedg(k,1) = -1
		sedg(k,2) = -1
		end if
		end if
		end do
	end if

	end do

c	Now assemble all new unique edges
	allocate(edgn(3*npolygon,2))
	nj=0.
	do n = 1, 3*npolygon
	if(sedg(n,1).ne.-1) then
	nj=nj+1
	edgn(nj,1) = sedg(n,1)
	edgn(nj,2) = sedg(n,2)
c	print *,'Edge= ',nj,' vertices= ',edgn(nj,1),edgn(nj,2)
c	if(modulo(nj,20).eq.1) pause
	end if
	end do
	deallocate(sedg)
	nedge = nj
	print *, nj, ' New unique edges found'
	print *, " "
	go to 100
200	continue

	avarea = totarea/npoly
	print *, '  '
	write(6,*) "COALESCE PROCESS RESULTS -Version= ",ver
	write(6,*) " "
		if(sanner) then
		write(6,*) 'Input File (Sanner): ', tfile 
		else if(connolly) then
		write(6,*) 'Input File (Connolly): ', tfile
		else if(smart) then
		write(6,*) 'Input File (SMART): ', tfile
		end if
	print *, " Edge fraction = ", edgfrac
	print *, " "
	print *, " Beginning polygon Area =  ",beginarea
	print *, 'Number of triangle records, pass 1= ',ntriang
	print *, 'Number of degenerate triangles, pass 1= ',idegen
	print *, 'Number of small triangles deleted, pass 1= ',idelete1
	print *, 'Number of polygons generated, pass 1= ', npolygon1
	print *,'Total area, pass 1 = ',totarea1
	print *,' '
	print *, 'Number of small triangles deleted, pass 2= ',idelete
	print *, 'Number of polygons generated, pass 2= ', npoly
	print *,'Total area, pass 2 = ',totarea
	print *, 'Average Triangle Area (A^2)= ', avarea
	print *,'Average edgelength (A)= ', avelength/nedge
	print *, 'Cut-off length (A) = ', cutofflength
	print *,'Found',nobtuse,' slender obtuse triangles, processed',icon
      print *,'Collapses skipped by the topology check: ',nlinkskip
      print *,'Flips skipped by the topology check: ',nflipskip


c	OUTPUT AND HISTOGRAM
                
C      APPEND NUMBER OF POLYGONS TO END OF OUTPUTFILE
       NPOLYX = NPOLY
       CALL APPEND(basename,outputf,NPOLYX,L)
       output = outputf(1:L+5) 
       outputhis = outputf(1:L+5)//'.his'  



c	Output is suppressed if the np flag has value -1
	if(np.eq.-1) goto 250
c     A triangle with two vertices at the same position has zero area,
c     and BEST gives NaN for it.  msroll can give such vertices where
c     its surface folds onto itself.  A triangulation with one is not
c     written; rcoal then reports an error.  The positions are compared
c     as written to the output file, which is what BEST reads.
      nzero = 0
      do m=1,npoly
      write(cv1,'(3e17.8)') (v1(m,n),n=1,3)
      write(cv2,'(3e17.8)') (v2(m,n),n=1,3)
      write(cv3,'(3e17.8)') (v3(m,n),n=1,3)
      if(cv1.eq.cv2.or.cv2.eq.cv3.or.cv1.eq.cv3) nzero = nzero+1
      end do
      if(nzero.gt.0) then
      print *,'Triangles with two vertices at the same position: ',
     $ nzero,'; not written'
      goto 250
      endif

       open(15,file=outputhis)
c      open(14,file=output)
c	 open(15,file=outputhis)

	write(15,*) "COALESCE PROCESS RESULTS -Version= ",ver
	write(15,*) " "
		if(sanner) then
		write(15,*) 'Input File (Sanner): ',tfile 
		else if(connolly) then
		write(15,*) 'Input File (Connolly): ',tfile 
		else if(smart) then
		write(15,*) 'Input File (SMART): ',tfile 
		end if

	write(15,*) " Edge fraction = ", edgfrac
	write(15,*) " "
	write(15,*) " Beginning polygon Area =  ",beginarea
	write(15,*) 'Initial number of triangle records= ',ntriang
	write(15,*) 'Number of degenerate triangles, pass 1= ',idegen
	write(15,*) 'Number of small triangles deleted, pass 1= ',idelete1
	write(15,*) 'Number of polygons generated, pass 1= ', npolygon1
	write(15,*) 'total area, pass 1 = ',totarea1
	write(15,*) ' '
	write(15,*) 'Number of small triangles deleted, pass 2= ',idelete
	write(15,*) 'Number of polygons generated, pass 2= ', npoly
	write(15,*) 'total area, pass 2 = ',totarea
	write(15,*) 'Average Triangle Area (A^2)= ', avarea
	write(15,*) 'Average edgelength (A)= ', avelength/nedge
	write(15,*) 'Cut-off length (A) = ', cutofflength
	write(15,*) 'Found',nobtuse,' slender obtuse triangles, 
	1processed',icon
	write(15,*) ' '

c      Call bin program and dstrb to create histogram
c       This creates a bimodal distribution of points centered
c      around 10% of average triangle area and 1.0   A**2

      write(15,*)"Histogram bins and count follows:"
c     The counts start from zero; they were not set before.
      do n=1,50
      ndist1(n) = 0
      ndist2(n) = 0
      end do
      numbin=12
      call bin(0.1*avarea,0.5,bins,numbin+1)
c     do n=1,numbin
c	write(15,*)bins(n)
c	end do
      do m=1,npoly
	if (area(m) <= 0.2*avarea ) then
       call dstrb(area(m),ndist1,bins,numbin)
	end if
	end do
	write(15,*)"Area distribution around",.1*avarea,' A^2'
	nsmallarea = 0
      do n=2,numbin
	write(15,*)bins(n),ndist1(n)
	nsmallarea = nsmallarea + ndist1(n)
	end do
	if(nsmallarea.gt.0) then
	print *," "
	print *, "************************"
	print *, 'N =  ',npoly
	print *,'ALERT! ',nsmallarea,' Triangles with Area<
	1 20% average area Remain!'
	print *,'Review histogram.' 
	print *, 'Edgefrac = ', edgfrac
	print *, "************************"
	endif

	numbin=22
      call bin(avarea,0.2,bins2,numbin+1)
c     do n=1,numbin
c	write(15,*)bins2(n)
c	end do
      do m=1,npoly
	if (area(m) > 0.2*avarea ) then
       call dstrb(area(m),ndist2,bins2,numbin)
	end if
	end do
      write(15,*)"Area distribution around", avarea,'A^2'
      do n=2,numbin
	write(15,*)bins2(n),ndist2(n)
	end do
	close(15)


c     Write output file for the hydrodynamics program:
c
        open(14,file=output)
	write(14,*) npoly,truearea,TRUEVOL
       do n=1,npoly
      write(14,12)v1(n,1),v1(n,2),v1(n,3),v2(n,1),v2(n,2),v2(n,3),
     $    v3(n,1),v3(n,2),v3(n,3)
12      format(9e17.8)
        end do
      close(14)
c     rcoal lists the files at the end, once all are written.
      nout = nout+1
      outf(nout) = output
      nout = nout+1
      outf(nout) = outputhis

250	continue

	return

      END subroutine


       subroutine Bin(e,s,bins,numbin)
c      This program produces a set of values, placed in the array bins, which serve
c      as separators in producing a distribution (see subroutine Dstrb).
c      The value e specifies the value of the independent variable at the mid point of
c      the distribution, while the value s defines the spacing between the separators
c      or bins. 

c       implicit double precision (a-h,o-z)
         
        dimension bins(numbin)
        do  j=1,(numbin-1)/2
        bins(j)=e-((numbin+1)/2-j)*s*e/2.5
        bins(j+(numbin+1)/2)=e+j*s*e/2.5
        end do
        bins((numbin+1)/2)=e
        return
       end subroutine


        Subroutine Dstrb(z,ndist,bins,numbin)
c     The quantity z is placed in the histogram ndist, given the
c     set of  comparison points contained in array bins.
c     The actual bin width is bin(j+1) - bin(j).  The bin separators
c      are produced in program bin.

c       implicit double precision (a-h,o-z)

        dimension ndist(numbin),bins(numbin)
        nc=0
        do  10 j=1,(numbin-2)/2
       if (z.le.bins(-j+1+(numbin/2)).and.z.gt.bins((numbin/2)-j))then
        ndist((numbin/2)+1-j)=ndist((numbin/2)+1-j)+1
        nc=1
       goto 15
        endif
        if(z.gt.bins((numbin/2)-1+j).and.z.le.bins((numbin/2)+j))then
        ndist((numbin/2)-1+j)=ndist((numbin/2)-1+j)+1
        nc=1
	  goto 15
       endif
10     continue
        if(nc.eq.0) then
        if (z.le.bins(1)) then
        ndist(1)=ndist(1)+1
        else if  (z.gt.bins(numbin)) then
        ndist(numbin)=ndist(numbin)+1
	else if (z.le.bins(numbin).and.z.gt.bins(numbin-1)) then
		ndist(numbin-1) = ndist(numbin-1)+1
        endif
        endif
 15     continue

        return
        end subroutine

	function dot(a,b)
	dimension a(3),b(3)
	dot = a(1)*b(1)+a(2)*b(2)+a(3)*b(3)
	return
	end function
	
	subroutine areas(v1,v2,v3,area)
	dimension v1(3),v2(3),v3(3)
	v21x = v2(1)-v1(1)
	v21y = v2(2)-v1(2)
	v21z = v2(3)-v1(3)
	v31x = v3(1)-v1(1)
	v31y = v3(2)-v1(2)
	v31z = v3(3)-v1(3)
	area = Sqrt((v31y*v21z-v31z*v21y)**2 + 
     $ (v31z*v21x-v31x*v21z)**2+(v31x*v21y-v31y*v21x)**2)/2.
	return
	end

	subroutine slength(v1,v2,slen)
	dimension v1(3),v2(3)
	v21x = v2(1)-v1(1)
	v21y = v2(2)-v1(2)
	v21z = v2(3)-v1(3)
	slen = sqrt(v21x*v21x+v21y*v21y+v21z*v21z)
	return
	end

      subroutine slender_del(v1,v2,x,nvx,vn,ntrg,idelete)
c     Collapses the small edge v1-v2: deletes the triangles on the edge
c     (setting their vertex numbers to -1) and moves the vertex shared by
c     fewer triangles onto the other.  x holds the nvx vertex coordinates.
c     The collapse is skipped, and counted in /lnkcnt/, unless exactly two
c     triangles share the edge and v1, v2 have no common neighbours other
c     than those two triangles' third vertices.  A collapse that fails
c     this test glues two triangles together, leaving duplicate triangles
c     and non-manifold edges, which later deletions can turn into holes.
c     A zero-length edge is always collapsed (see below).
      integer v1,v2,first,second,vn,tag1,tag2,nb1,nb2,iapex
      dimension vn(ntrg,3),x(nvx,3),tag1(ntrg),tag2(ntrg)
      dimension nb1(2*ntrg),nb2(2*ntrg),iapex(2)
      logical has1,has2,dupl,found
      common /lnkcnt/ nlinkskip

c     Find the triangles on v1 and on v2, their other vertices, and the
c     third vertices of the triangles on the edge.
      first = 0
      second = 0
      nboth = 0
      n1 = 0
      n2 = 0
      do k = 1,ntrg
      if(vn(k,1).ne.-1) then
      has1 = v1.eq.vn(k,1).or.v1.eq.vn(k,2).or.v1.eq.vn(k,3)
      has2 = v2.eq.vn(k,1).or.v2.eq.vn(k,2).or.v2.eq.vn(k,3)
      if(has1) then
      first = first+1
      tag1(first) = k
      do m = 1,3
      if(vn(k,m).ne.v1.and.vn(k,m).ne.v2) then
      n1 = n1+1
      nb1(n1) = vn(k,m)
      endif
      enddo
      endif
      if(has2) then
      second = second+1
      tag2(second) = k
      do m = 1,3
      if(vn(k,m).ne.v1.and.vn(k,m).ne.v2) then
      n2 = n2+1
      nb2(n2) = vn(k,m)
      endif
      enddo
      endif
      if(has1.and.has2) then
      nboth = nboth+1
      if(nboth.le.2) then
      do m = 1,3
      if(vn(k,m).ne.v1.and.vn(k,m).ne.v2) iapex(nboth) = vn(k,m)
      enddo
      endif
      endif
      endif
      enddo

c     A vertex with no triangles left: nothing to do.
      if(first.eq.0.or.second.eq.0) return

c     The topology check, for edges of nonzero length.  msroll can give
c     two vertex numbers the same position; the edge between them is
c     always collapsed, as before: that changes no geometry, while
c     skipping it would leave triangles of zero area.
      d2 = (x(v1,1)-x(v2,1))**2+(x(v1,2)-x(v2,2))**2
     $   +(x(v1,3)-x(v2,3))**2
      if(d2.gt.1.0e-10) then
      if(nboth.ne.2) goto 900
      if(iapex(1).eq.iapex(2)) goto 900
      ncommon = 0
      do i = 1,n1
      dupl = .false.
      do j = 1,i-1
      if(nb1(j).eq.nb1(i)) then
      dupl = .true.
      exit
      endif
      enddo
      if(.not.dupl) then
      found = .false.
      do j = 1,n2
      if(nb2(j).eq.nb1(i)) then
      found = .true.
      exit
      endif
      enddo
      if(found) ncommon = ncommon+1
      endif
      enddo
      if(ncommon.ne.2) goto 900
      endif

c     Delete the triangles on the edge.
      do i = 1,first
      k = tag1(i)
      if(v2.eq.vn(k,1).or.v2.eq.vn(k,2).or.v2.eq.vn(k,3)) then
      idelete = idelete+1
      vn(k,1) = -1
      vn(k,2) = -1
      vn(k,3) = -1
      endif
      enddo

c     Coalesce in the direction of the vertex shared by more triangles.
      if(second.ge.first) then
      do i = 1,first
      k = tag1(i)
      if(vn(k,1).eq.v1) then
      vn(k,1) = v2
      else if(vn(k,2).eq.v1) then
      vn(k,2) = v2
      else if(vn(k,3).eq.v1) then
      vn(k,3) = v2
      endif
      enddo
      else
      do i = 1,second
      k = tag2(i)
      if(vn(k,1).eq.v2) then
      vn(k,1) = v1
      else if(vn(k,2).eq.v2) then
      vn(k,2) = v1
      else if(vn(k,3).eq.v2) then
      vn(k,3) = v1
      endif
      enddo
      endif
      return

900   nlinkskip = nlinkskip+1
      return
      end

	subroutine dangle_del(ntrg,vn,verb,idang)
c	This subroutine identifies the triangles with vertices shared with two
c	or less triangles. Deletes the triangles with no shared vertices = isolated.
c	Deletes triangles which have singly shared vertices = dangling.
c	Provides statistics on the average number of triangles sharing a vertex.

	integer vn,v1,v2,v3
	dimension vn(ntrg,3)
	logical verb
			
	isodel = 0
	idang = 0
	nshare1t = 0
	nshare2t = 0
	nshare3t = 0
	npoly = 0

	do j = 1,ntrg
		if(vn(j,1).ne.-1) then
		npoly = npoly +1
		v1 = vn(j,1)
		v2 = vn(j,2)
		v3 = vn(j,3)
		nshare1 = 0
		nshare2 = 0
		nshare3 = 0
		
c	Find the number of times each vertex is shared with some other triangle.
		do k=1,ntrg
		if(k.ne.j.and.v1.ne.-1) then

		if(v1.eq.vn(k,1).or.v1.eq.vn(k,2).or.v1.eq.vn(k,3)) then
		nshare1 = nshare1 +1
		endif
		if(v2.eq.vn(k,1).or.v2.eq.vn(k,2).or.v2.eq.vn(k,3)) then
		nshare2 = nshare2 +1
		endif
		if(v3.eq.vn(k,1).or.v3.eq.vn(k,2).or.v3.eq.vn(k,3)) then
		nshare3 = nshare3 +1
		endif

		endif
		enddo
c	Analyse the nshare values and take action.

c	If any vertices are not shared, the triangle is dangling or isolated
c	and will be deleted.
		if(nshare1.eq.0.or.nshare2.eq.0.or.nshare3.eq.0) then
		vn(j,1) = -1
		vn(j,2) = -1
		vn(j,3) = -1
		idang = idang +1
		if(verb) print *, "Dangling or isolated triangle deleted", j
		endif
c	If all nsharej numbers are zero, the triangle is not connected
		if(nshare1.eq.0.and.nshare2.eq.0.and.nshare3.eq.0) then
		isodel = isodel +1
		if(verb) print *, "The deleted triangle was Isolated", j
		endif

c	If nsharej= 1, then there is a two-triangle fin.  Delete.
		if(nshare1.eq.1.or.nshare2.eq.1.or.nshare3.eq.1) then
		vn(j,1) = -1
		vn(j,2) = -1
		vn(j,3) = -1
		idang = idang +1
		if(verb) print *, "Double-Dangling triangle deleted", j
		endif
			
		nshare1t = nshare1t + nshare1
		nshare2t = nshare2t + nshare2
		nshare3t = nshare3t + nshare3
		endif				
	enddo
		share1t = nshare1t/npoly
		share2t = nshare2t/npoly
		share3t = nshare3t/npoly
	print *, "Average shared vertices=",share1t,share2t,share3t
	print *, " Deleted", idang," dangling triangles."
	print *, "Of these,",isodel," where isolated."
	print *, " "
		
	return
	end

       SUBROUTINE APPEND(basename,outputf,NPOLYX,L)

       CHARACTER*80 STRING,basename,outputf

       DO I = 0,4

       IREMAINDER = MOD(NPOLYX, 10)

       IF (IREMAINDER == 0) THEN
       STRING (5-I:5-I) = '0'
       ELSE IF (IREMAINDER == 1) THEN
       STRING (5-I:5-I) = '1'
       ELSE IF (IREMAINDER == 2) THEN
       STRING (5-I:5-I) = '2'
       ELSE IF (IREMAINDER == 3) THEN
       STRING (5-I:5-I) = '3'
       ELSE IF (IREMAINDER == 4) THEN
       STRING (5-I:5-I) = '4'
       ELSE IF (IREMAINDER == 5) THEN
       STRING (5-I:5-I) = '5'
       ELSE IF (IREMAINDER == 6) THEN
       STRING (5-I:5-I) = '6'
       ELSE IF (IREMAINDER == 7) THEN
       STRING (5-I:5-I) = '7'
       ELSE IF (IREMAINDER == 8) THEN
       STRING (5-I:5-I) = '8'
       ELSE
       STRING (5-I:5-I) = '9'
       ENDIF

       NPOLYX = NPOLYX/10                       
                                               
       ENDDO                                    
                                                 
      outputf = basename(1:L-1)//'_'//STRING                

       RETURN           
       END               

