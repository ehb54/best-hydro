! BEST's boundary-element matrix GR and the solves with it.
!
! Dense (the default up to NHMAT triangles, or with -dn): GR (3*NAR x 3*NAR)
! is filled by TGETGR, deflated and factored into LU form by DGETRF in the
! main program, and BESOLV solves with those factors (DGETRS).
!
! H-matrix (above NHMAT triangles, or with -hm): HMPREP sets up the solver in
! hmatrix/ (C++, called through bind(C)), which never forms GR.  It computes
! the same entries, compresses the blocks of GR between well-separated groups
! of triangles to low rank at a relative accuracy eps (-he, default 1e-4), and
! solves by GMRES with a coarse H-LU preconditioner, or with a direct H-LU
! factorization (-hlu).  Its memory grows about as N log N instead of 72 N^2
! bytes.  BEST's results then agree with the dense ones to about 1e-5.  With
! -hs it stores the blocks in single precision, in about half the memory.
!
! Either way DIFCO and VISCO call BESOLV, for several right-hand sides at once.
!
! E. H. Brookes, 11/13/2026.

module besolve
   use, intrinsic :: iso_c_binding
   implicit none
   private
   public :: gr, ipiv, besolv, hmavail, hmprep, nhmat

   ! GR and, after DGETRF, its LU factors and pivots (dense only).
   real, allocatable :: gr(:,:)
   integer, allocatable :: ipiv(:)

   ! Above this many triangles BEST uses the H-matrix solver unless -dn is given.
   integer, parameter :: nhmat = 10000

   ! The GMRES preconditioner's accuracy.
   real, parameter :: epspre = 0.1

   ! Set by HMPREP: BESOLV then uses the H-matrix solver.  Its messages go to logfile.
   logical :: hm = .false.
   character(len=80) :: logfile = ' '

   ! With -v or -vc (H-matrix solver only): VISCO first solves for the five strain flows about
   ! the centroid, STRAIN(X).  BESOLV solves for them already with DIFCO's six flows, in one
   ! GMRES run of 11 right-hand sides (its time is set by the preconditioner's triangular
   ! solves, about the same for 11 as for 6), and keeps them: srhs the right-hand sides, ssol
   ! the solutions.  hmx is X.
   logical :: hmvisc = .false.
   real, allocatable :: hmx(:,:), srhs(:,:), ssol(:,:)

   interface
      integer(c_int) function bhm_available() bind(C, name='bhm_available')
         import :: c_int
      end function bhm_available
      integer(c_int) function bhm_create(n, v1, v2, v3, area, afac, cg, threads) &
            bind(C, name='bhm_create')
         import :: c_int, c_double
         integer(c_int), value :: n, threads
         real(c_double), intent(in) :: v1(*), v2(*), v3(*), area(*)
         real(c_double), value :: afac, cg
      end function bhm_create
      subroutine bhm_incenters(x) bind(C, name='bhm_incenters')
         import :: c_double
         real(c_double), intent(out) :: x(*)
      end subroutine bhm_incenters
      subroutine bhm_entries(i0, m, j0, nc, a, lda) bind(C, name='bhm_entries')
         import :: c_int, c_double
         integer(c_int), value :: i0, m, j0, nc, lda
         real(c_double), intent(out) :: a(*)
      end subroutine bhm_entries
      real(c_double) function bhm_diag_mean() bind(C, name='bhm_diag_mean')
         import :: c_double
      end function bhm_diag_mean
      integer(c_int) function bhm_build(eps, direct, eps_pre, gdef, xmu, single) bind(C, name='bhm_build')
         import :: c_int, c_double
         real(c_double), value :: eps, eps_pre, gdef
         integer(c_int), value :: direct, single
         real(c_double), intent(in) :: xmu(*)
      end function bhm_build
      integer(c_int) function bhm_solve(nrhs, b) bind(C, name='bhm_solve')
         import :: c_int, c_double
         integer(c_int), value :: nrhs
         real(c_double), intent(inout) :: b(*)
      end function bhm_solve
      integer(c_int) function bhm_log_line(buf, len) bind(C, name='bhm_log_line')
         import :: c_int, c_char
         character(kind=c_char), intent(out) :: buf(*)
         integer(c_int), value :: len
      end function bhm_log_line
   end interface

contains

   ! B(3*NAR, NRHS) holds the right-hand sides on entry and the solutions on
   ! return.  INFO is DGETRS's, or the H-matrix solver's (which stops BEST).
   subroutine besolv(nrhs, b, info)
      integer, intent(in) :: nrhs
      real, intent(inout) :: b(*)
      integer, intent(out) :: info
      integer :: n, nar
      real, allocatable :: bb(:,:)

      if (hm) then
         n = 3*size(hmx, 1)
         nar = size(hmx, 1)
         if (hmvisc .and. nrhs == 6 .and. .not. allocated(srhs)) then
            allocate(bb(n, 11), srhs(n, 5), ssol(n, 5))
            bb(:, 1:6) = reshape(b(1:6*n), [n, 6])
            call strain(nar, hmx, srhs)
            bb(:, 7:11) = srhs
            info = bhm_solve(11, bb)
            b(1:6*n) = reshape(bb(:, 1:6), [6*n])
            ssol = bb(:, 7:11)
            deallocate(bb)
         else if (nrhs == 5 .and. allocated(ssol)) then
            if (all(b(1:5*n) == reshape(srhs, [5*n]))) then
               b(1:5*n) = reshape(ssol, [5*n])
               info = 0
               open(12, file=logfile, status='old', position='append')
               write(12,'(1x,a)') 'solve: the 5 strain flows were solved with the 6 above'
               close(12)
            else
               info = bhm_solve(nrhs, b)
            end if
            deallocate(srhs, ssol)
         else
            info = bhm_solve(nrhs, b)
         end if
         call hmlog()
         if (info /= 0) call hmstop('the H-matrix solve failed')
         return
      end if
      n = size(gr, 1)
      call dgetrs('N', n, nrhs, gr, n, ipiv, b, n, info)
   end subroutine besolv

   ! True if BEST was built with the H-matrix solver.
   logical function hmavail()
      hmavail = bhm_available() /= 0
   end function hmavail

   ! Sets up the H-matrix solver: what best.f does from TGETGR to DGETRF on
   ! the dense path.  X receives the incenters (TGETGR's X).  The log gets
   ! the same lines about GR and its deflation as on the dense path, then the
   ! solver's own.  HEPS is the H-matrix accuracy; IHLU = 1 asks for the
   ! direct H-LU solver instead of GMRES; IVISCO = 1 (-v, -vc) has BESOLV solve
   ! VISCO's strain flows with DIFCO's; IHS = 1 (-hs) stores the blocks in
   ! single precision.
   subroutine hmprep(nar, afac, cg, v1, v2, v3, x, area, outfnm2, heps, ihlu, ivisco, ihs)
      integer, intent(in) :: nar, ihlu, ivisco, ihs
      real, intent(in) :: afac, cg, heps
      real, intent(in) :: v1(nar,3), v2(nar,3), v3(nar,3), area(nar)
      real, intent(out) :: x(nar,3)
      character(len=*), intent(in) :: outfnm2
      real, allocatable :: xmu(:)
      real :: a(9,9), xmu2, diagm, gdef
      real :: v21x, v21y, v21z, v31x, v31y, v31z, xnx, xny, xnz, xnl
      integer :: n, j, m, info
      character(len=8) :: val
      integer :: st

      logfile = outfnm2
      info = bhm_create(nar, v1, v2, v3, area, afac, cg, 0)
      call hmlog()
      if (info /= 0) call hmstop('the H-matrix solver could not be set up')
      call bhm_incenters(x)

      ! For the tests: with BEST_HMATRIX_CHECK set, also fill GR with TGETGR
      ! and compare.
      call get_environment_variable('BEST_HMATRIX_CHECK', val, status=st)
      if (st == 0 .and. len_trim(val) > 0 .and. val /= '0') &
         call hmcheck(nar, afac, cg, v1, v2, v3, x, area)

      open(12, file=logfile, status='old', position='append')
      m = min(9, 3*nar)
      call bhm_entries(0, m, 0, m, a, 9)
      write(12,*) 'first few rows of gr matrix'
      do j = 1, m
         write(12,*) a(j,1), a(j,2), a(j,3), a(j,4), a(j,5)
      end do
      write(12,*) 'first few columns of gr matrix'
      do j = 1, m
         write(12,*) a(1,j), a(2,j), a(3,j), a(4,j), a(5,j)
      end do
      write(12,*) ' '

      ! The deflation of the normal field, as in best.f (see the comments
      ! there); the H-matrix solver adds GDEF*XMU*XMU^T to every block.
      allocate(xmu(3*nar))
      xmu2 = 0.0
      do n = 1, nar
         v21x = v2(n,1) - v1(n,1)
         v21y = v2(n,2) - v1(n,2)
         v21z = v2(n,3) - v1(n,3)
         v31x = v3(n,1) - v1(n,1)
         v31y = v3(n,2) - v1(n,2)
         v31z = v3(n,3) - v1(n,3)
         xnx = v31y*v21z - v31z*v21y
         xny = v31z*v21x - v31x*v21z
         xnz = v31x*v21y - v31y*v21x
         xnl = sqrt(xnx**2 + xny**2 + xnz**2)
         xmu(3*n-2) = area(n)*xnx/xnl
         xmu(3*n-1) = area(n)*xny/xnl
         xmu(3*n) = area(n)*xnz/xnl
         xmu2 = xmu2 + xmu(3*n-2)**2 + xmu(3*n-1)**2 + xmu(3*n)**2
      end do
      diagm = bhm_diag_mean()
      gdef = diagm/xmu2
      write(12,*) 'normal-field deflation: mean diagonal of GR = ', diagm
      write(12,*) ' added GDEF*XMU*XMU^T to GR, GDEF = ', gdef
      write(12,*) ' '
      if (ihlu == 1) then
         write(12,*) 'H-matrix solver, direct H-LU, eps = ', heps
      else
         write(12,*) 'H-matrix solver, GMRES, eps = ', heps
      end if
      if (ihs == 1) write(12,*) 'H-matrix stored in single precision (-hs)'
      close(12)

      info = bhm_build(heps, ihlu, epspre, gdef, xmu, ihs)
      call hmlog()
      if (info /= 0) call hmstop('the H-matrix could not be built')
      deallocate(xmu)
      hm = .true.
      hmvisc = ivisco == 1
      hmx = x
   end subroutine hmprep

   ! Fills GR with TGETGR, as on the dense path, and logs how many of its
   ! entries, and of TGETGR's incenters, differ in any bit from the H-matrix
   ! solver's.  They should all be the same.
   subroutine hmcheck(nar, afac, cg, v1, v2, v3, x, area)
      integer, intent(in) :: nar
      real, intent(in) :: afac, cg, x(nar,3)
      real, intent(in) :: v1(nar,3), v2(nar,3), v3(nar,3), area(nar)
      real, allocatable :: g(:,:), x2(:,:), c(:,:)
      integer :: n3, nc, j0, i, j
      integer(8) :: ngr, nx

      n3 = 3*nar
      nc = 3*min(nar, 256)
      allocate(g(n3,n3), x2(nar,3), c(n3,nc))
      call tgetgr(nar, afac, cg, v1, v2, v3, x2, area, g)
      nx = 0
      do j = 1, 3
         do i = 1, nar
            if (transfer(x(i,j), 0_8) /= transfer(x2(i,j), 0_8)) nx = nx + 1
         end do
      end do
      ngr = 0
      do j0 = 0, n3 - 1, nc
         call bhm_entries(0, n3, j0, min(nc, n3 - j0), c, n3)
         do j = 1, min(nc, n3 - j0)
            do i = 1, n3
               if (transfer(c(i,j), 0_8) /= transfer(g(i,j0+j), 0_8)) ngr = ngr + 1
            end do
         end do
      end do
      open(12, file=logfile, status='old', position='append')
      write(12,'(1x,a,i0,a,i0,a,i0,a,i0,a)') 'GR check: ', ngr, ' of ', int(n3,8)**2, &
         ' entries and ', nx, ' of ', 3*nar, ' incenter coordinates differ from TGETGR''s'
      close(12)
      deallocate(g, x2, c)
   end subroutine hmcheck

   ! Appends the H-matrix solver's messages to the log.  Unit 12 must not be
   ! open in the caller.
   subroutine hmlog()
      character(kind=c_char, len=512) :: line
      integer :: n

      open(12, file=logfile, status='old', position='append')
      do
         n = bhm_log_line(line, len(line))
         if (n < 0) exit
         write(12,'(1x,a)') line(1:n)
      end do
      close(12)
   end subroutine hmlog

   subroutine hmstop(msg)
      character(len=*), intent(in) :: msg

      print *, 'BEST: ', msg, '; see ', trim(logfile)
      stop 1
   end subroutine hmstop

end module besolve
