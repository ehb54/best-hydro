! Stand-ins for the H-matrix solver (hmatrix/best_hmatrix.h), linked instead
! of it when BEST is built without C++ (make HMATRIX=).  BEST then uses the
! dense solver at every size, and -hm stops with a message.
!
! E. H. Brookes, 11/13/2026.

integer(c_int) function bhm_available() bind(C, name='bhm_available')
   use, intrinsic :: iso_c_binding
   implicit none
   bhm_available = 0
end function bhm_available

integer(c_int) function bhm_create(n, v1, v2, v3, area, afac, cg, threads) bind(C, name='bhm_create')
   use, intrinsic :: iso_c_binding
   implicit none
   integer(c_int), value :: n, threads
   real(c_double), intent(in) :: v1(*), v2(*), v3(*), area(*)
   real(c_double), value :: afac, cg
   bhm_create = 1
end function bhm_create

subroutine bhm_incenters(x) bind(C, name='bhm_incenters')
   use, intrinsic :: iso_c_binding
   implicit none
   real(c_double), intent(out) :: x(*)
end subroutine bhm_incenters

subroutine bhm_entries(i0, m, j0, nc, a, lda) bind(C, name='bhm_entries')
   use, intrinsic :: iso_c_binding
   implicit none
   integer(c_int), value :: i0, m, j0, nc, lda
   real(c_double), intent(out) :: a(*)
end subroutine bhm_entries

real(c_double) function bhm_diag_mean() bind(C, name='bhm_diag_mean')
   use, intrinsic :: iso_c_binding
   implicit none
   bhm_diag_mean = 0
end function bhm_diag_mean

integer(c_int) function bhm_build(eps, direct, eps_pre, gdef, xmu, single) bind(C, name='bhm_build')
   use, intrinsic :: iso_c_binding
   implicit none
   real(c_double), value :: eps, eps_pre, gdef
   integer(c_int), value :: direct, single
   real(c_double), intent(in) :: xmu(*)
   bhm_build = 1
end function bhm_build

integer(c_int) function bhm_solve(nrhs, b) bind(C, name='bhm_solve')
   use, intrinsic :: iso_c_binding
   implicit none
   integer(c_int), value :: nrhs
   real(c_double), intent(inout) :: b(*)
   bhm_solve = 1
end function bhm_solve

! One message, once.
integer(c_int) function bhm_log_line(buf, len) bind(C, name='bhm_log_line')
   use, intrinsic :: iso_c_binding
   implicit none
   character(kind=c_char), intent(out) :: buf(*)
   integer(c_int), value :: len
   character(len=*), parameter :: msg = 'error: this BEST was built without the H-matrix solver (make HMATRIX=)'
   logical, save :: said = .false.
   integer :: i

   bhm_log_line = -1
   if (said) return
   said = .true.
   bhm_log_line = min(len, len_trim(msg))
   do i = 1, bhm_log_line
      buf(i) = msg(i:i)
   end do
end function bhm_log_line

subroutine bhm_destroy() bind(C, name='bhm_destroy')
end subroutine bhm_destroy
