! The connected surfaces of BEST's triangulation, for the normal-field
! deflation (best.f, besolve.f90).
!
! On every closed surface the Stokes single layer of that surface's normal
! field vanishes, so a triangulation of K separate closed surfaces gives GR K
! near-null vectors, one per surface.  One deflation term over all triangles
! removes only their sum; the other K-1 stay, and GR is then singular to the
! discretization error.  On a probe-0 fragment of lysozyme with a second
! surface of 130 triangles the 1-norm condition estimate was 4e11 (2e6 for the
! main surface alone), and the dense D_t and [eta] were off by 3% and 27%.
!
! COMPON finds the surfaces: triangles that share a vertex (the same three
! coordinates, as read) belong to one.  ICOMP(N) is the surface of triangle
! N, numbered 1..NCOMP in the order of their first triangles.  NOPEN counts
! the edges not shared by exactly two triangles (0 for closed surfaces).
!
! E. H. Brookes, 11/15/2026.

subroutine compon(nar, v1, v2, v3, icomp, ncomp, nopen)
   implicit none
   integer, intent(in) :: nar
   real, intent(in) :: v1(nar,3), v2(nar,3), v3(nar,3)
   integer, intent(out) :: icomp(nar), ncomp, nopen
   real, allocatable :: key(:,:)
   integer, allocatable :: idx(:), vid(:), root(:), lab(:)
   integer :: nv, k, j, n, c, r, a, b, nid

   ! Vertices: k = 3*(n-1) + c is corner c of triangle n.
   nv = 3*nar
   allocate(key(3,nv), idx(nv), vid(nv))
   do n = 1, nar
      key(:,3*n-2) = v1(n,:)
      key(:,3*n-1) = v2(n,:)
      key(:,3*n) = v3(n,:)
   end do
   call msort3(nv, key, idx)
   nid = 0
   do j = 1, nv
      if (j == 1) then
         nid = 1
      else if (any(key(:,idx(j)) /= key(:,idx(j-1)))) then
         nid = nid + 1
      end if
      vid(idx(j)) = nid
   end do

   ! Union-find over triangles that share a vertex.
   allocate(root(nar))
   do n = 1, nar
      root(n) = n
   end do
   do j = 2, nv
      if (vid(idx(j)) == vid(idx(j-1))) call unite((idx(j)+2)/3, (idx(j-1)+2)/3)
   end do
   allocate(lab(nar))
   lab = 0
   ncomp = 0
   do n = 1, nar
      r = find(n)
      if (lab(r) == 0) then
         ncomp = ncomp + 1
         lab(r) = ncomp
      end if
      icomp(n) = lab(r)
   end do

   ! Edges: (smaller, larger) vertex id of each side; closed means each pair twice.
   do n = 1, nar
      do c = 1, 3
         a = vid(3*(n-1) + c)
         b = vid(3*(n-1) + mod(c, 3) + 1)
         key(1,3*(n-1)+c) = min(a, b)
         key(2,3*(n-1)+c) = max(a, b)
         key(3,3*(n-1)+c) = 0.0
      end do
   end do
   call msort3(nv, key, idx)
   nopen = 0
   k = 1
   do j = 2, nv + 1
      if (j > nv) then
         if (j - k /= 2) nopen = nopen + 1
      else if (any(key(:,idx(j)) /= key(:,idx(k)))) then
         if (j - k /= 2) nopen = nopen + 1
         k = j
      end if
   end do
   deallocate(key, idx, vid, root, lab)

contains

   integer function find(i)
      integer, intent(in) :: i
      integer :: p, q, t
      p = i
      do while (root(p) /= p)
         p = root(p)
      end do
      q = i
      do while (root(q) /= p)
         t = root(q)
         root(q) = p
         q = t
      end do
      find = p
   end function find

   subroutine unite(i, j)
      integer, intent(in) :: i, j
      integer :: ri, rj
      ri = find(i)
      rj = find(j)
      if (ri /= rj) root(max(ri, rj)) = min(ri, rj)
   end subroutine unite

end subroutine compon

! IDX(1..M) such that the columns KEY(:,IDX(J)) are in lexicographic order
! (a stable bottom-up merge sort).
subroutine msort3(m, key, idx)
   implicit none
   integer, intent(in) :: m
   real, intent(in) :: key(3,m)
   integer, intent(out) :: idx(m)
   integer, allocatable :: tmp(:)
   integer :: w, lo, mid, hi, i, j, k

   do i = 1, m
      idx(i) = i
   end do
   allocate(tmp(m))
   w = 1
   do while (w < m)
      lo = 1
      do while (lo <= m)
         mid = min(lo + w - 1, m)
         hi = min(lo + 2*w - 1, m)
         i = lo
         j = mid + 1
         k = lo
         do while (i <= mid .and. j <= hi)
            if (before(key(:,idx(j)), key(:,idx(i)))) then
               tmp(k) = idx(j)
               j = j + 1
            else
               tmp(k) = idx(i)
               i = i + 1
            end if
            k = k + 1
         end do
         do while (i <= mid)
            tmp(k) = idx(i)
            i = i + 1
            k = k + 1
         end do
         do while (j <= hi)
            tmp(k) = idx(j)
            j = j + 1
            k = k + 1
         end do
         lo = lo + 2*w
      end do
      idx = tmp
      w = 2*w
   end do
   deallocate(tmp)

contains

   logical function before(a, b)
      real, intent(in) :: a(3), b(3)
      if (a(1) /= b(1)) then
         before = a(1) < b(1)
      else if (a(2) /= b(2)) then
         before = a(2) < b(2)
      else
         before = a(3) < b(3)
      end if
   end function before

end subroutine msort3

! The deflation weights for NDEF surfaces: GDEFK(K)*|XMU(K)|^2 = DIAGM, with
! XMU(K) the part of XMU on the triangles of surface K (ICOMP).
subroutine defgam(nar, ndef, icomp, xmu, diagm, gdefk)
   implicit none
   integer, intent(in) :: nar, ndef, icomp(nar)
   real, intent(in) :: xmu(3*nar), diagm
   real, intent(out) :: gdefk(ndef)
   integer :: n

   gdefk = 0.0
   do n = 1, nar
      gdefk(icomp(n)) = gdefk(icomp(n)) + xmu(3*n-2)**2 + xmu(3*n-1)**2 + xmu(3*n)**2
   end do
   gdefk = diagm/gdefk
end subroutine defgam

! True if all of A(1..N) are finite (no NaN or Inf).
logical function allfin(n, a)
   implicit none
   integer, intent(in) :: n
   real, intent(in) :: a(n)
   integer :: i

   allfin = .true.
   do i = 1, n
      if (.not. (abs(a(i)) <= huge(a(i)))) allfin = .false.
   end do
end function allfin

! Stops BEST (exit code 1) because WHAT is not finite, with a line in the log
! LOGFILE and on standard output.
subroutine notfin(logfile, what)
   implicit none
   character(len=*), intent(in) :: logfile, what

   open(12, file=logfile, status='old', position='append')
   write(12,'(1x,4a)') 'error: ', what, ' are not finite (NaN or Inf); ', &
      'nothing was written to the .be file'
   close(12)
   print *, 'BEST: ', what, ' are not finite (NaN or Inf); see ', trim(logfile)
   stop 1
end subroutine notfin
