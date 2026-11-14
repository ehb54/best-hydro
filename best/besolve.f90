! BEST's boundary-element matrix GR and the solves with it.
!
! GR (3*NAR x 3*NAR) is filled by TGETGR, deflated and factored into LU form
! by DGETRF in the main program.  DIFCO and VISCO then call BESOLV, which
! solves GR X = B for several right-hand sides at once with those factors.
! Keeping the matrix here means the subroutines that need solutions do not
! hold it themselves.
!
! E. H. Brookes, 11/13/2026.

module besolve
   implicit none
   private
   public :: gr, ipiv, besolv

   ! GR and, after DGETRF, its LU factors and pivots.
   real, allocatable :: gr(:,:)
   integer, allocatable :: ipiv(:)

contains

   ! B(3*NAR, NRHS) holds the right-hand sides on entry and the solutions on
   ! return.  INFO is DGETRS's.
   subroutine besolv(nrhs, b, info)
      integer, intent(in) :: nrhs
      real, intent(inout) :: b(*)
      integer, intent(out) :: info
      integer :: n

      n = size(gr, 1)
      call dgetrs('N', n, nrhs, gr, n, ipiv, b, n, info)
   end subroutine besolv

end module besolve
