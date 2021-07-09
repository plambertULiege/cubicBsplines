!=======================================================================
subroutine tpower(x, t, p, ans)
!=======================================================================
! Input: x, t are double ; p is integer
!
! Output: (x - t) ^ p  is (x > t)
!
  implicit none

  integer p
  double precision x, t, ans
  ans = 0.
  if (x .gt. t) then
     ans = exp(p * log(x-t))
  end if
end subroutine tpower

!=======================================================================
double precision function dif1(x, j, nt, t)
!=======================================================================
! Compute first order divided difference of
!     fun(.) = max(([.] - t(j))^3, 0.)
! Thus, it yields
!  [t(j-1), t(j)] fun(.) = (fun(t(j)) - fun(t(j-1))) / (t(j) - t(j-1))
!
  integer, intent(in) :: nt, j
  double precision, intent(in) :: x, t(nt)
  double precision :: g1, g2
  g2 = 0.d0 ; g1 = 0.d0
  if (x > t(j)) g2 = (x-t(j))**3
  if (x > t(j-1)) g1 = (x-t(j-1))**3
  dif1 = (g2 - g1) / (t(j) - t(j-1))
end function dif1

!=======================================================================
double precision function integrated_dif1(x, j, nt, t)
!=======================================================================
! Compute the integral of the first order divided difference of
!     fun(.) = ([.] - t(j))^4/4  if ([.] > t(j))
! Thus, it yields the integral of
!  [t(j-1), t(j)] fun(.) = (fun(t(j)) - fun(t(j-1))) / (t(j) - t(j-1))
!
  integer, intent(in) :: nt, j
  double precision, intent(in) :: x, t(nt)
  double precision :: g1, g2
  g2 = 0.d0 ; g1 = 0.d0
  if (x > t(j)) g2 = ((x-t(j))**4)/4.d0
  if (x > t(j-1)) g1 = ((x-t(j-1))**4)/4.d0
  integrated_dif1 = (g2 - g1) / (t(j) - t(j-1))
end function integrated_dif1

!=======================================================================
double precision function D1_dif1(x, j, nt, t)
!=======================================================================
! Compute the 1st derivative of the first order divided difference of
!     fun(.) = max(([.] - t(j))^3, 0.)
! Thus, it yields the 1st derivative of
!  [t(j-1), t(j)] fun(.) = (fun(t(j)) - fun(t(j-1))) / (t(j) - t(j-1))
!
  integer, intent(in) :: nt, j
  double precision, intent(in) :: x, t(nt)
  double precision :: g1, g2
  g2 = 0.d0 ; g1 = 0.d0
  if (x > t(j)) g2 = 3.d0 * (x-t(j))**2
  if (x > t(j-1)) g1 = 3.d0 * (x-t(j-1))**2
  D1_dif1 = (g2 - g1) / (t(j) - t(j-1))
end function D1_dif1

!=======================================================================
double precision function D2_dif1(x, j, nt, t)
!=======================================================================
! Compute the 2ns derivative of the first order divided difference of
!     fun(.) = max(([.] - t(j))^3, 0.)
! Thus, it yields the 2nd derivative of
!  [t(j-1), t(j)] fun(.) = (fun(t(j)) - fun(t(j-1))) / (t(j) - t(j-1))
!
  integer, intent(in) :: nt, j
  double precision, intent(in) :: x, t(nt)
  double precision :: g1, g2
  g2 = 0.d0 ; g1 = 0.d0
  if (x > t(j)) g2 = 6.d0 * (x-t(j))
  if (x > t(j-1)) g1 = 6.d0 * (x-t(j-1))
  D2_dif1 = (g2 - g1) / (t(j) - t(j-1))
end function D2_dif1

!=======================================================================
subroutine cubicBsplines_general(nx, x, nknots, knots, B)
!=======================================================================
! Input: x(nx): values where to compute the cubic B-splines basis ;
!        knots(nknots): possibly unequidistant (but distinct) knots
!
! Output: B(nx, ndx+3)
!
! By P. Lambert (ULiege 2012) from de De Boor formula relating
!  cubic B-splines to truncated powers :
!    Bj(t) = (t(j+4) - t(j)) [t(j),...,t(j+4)] max(0,(.-t)^3)
!  with divided differences
!

  implicit none
  integer, intent(in) :: nx, nknots
  double precision, intent(in) :: x(nx), knots(nknots)
  double precision, intent(out) :: B(nx,nknots+2)

  integer :: nknots2, i, j
  double precision :: knots2(nknots+6), ff(3), dd(4)
  double precision :: eps1, eps2, dif1, temp

  ! Put 2x3 extra knots on the left and on the right of the knots range
  nknots2 = nknots + 6
  eps1 = knots(2) - knots(1)
  eps2 = knots(nknots) - knots(nknots-1)
  knots2 = (/ knots(1)-3.d0*eps1, knots(1)-2.d0*eps1, knots(1)-eps1, &
              knots(1:nknots), &
              knots(nknots)+eps2, knots(nknots)+2.d0*eps2, knots(nknots)+3.d0*eps2 /)
  B(:,:) = 0.d0
  do i = 1, nx
     do j = 1, nknots2-4
        temp = 0.
        if (x(i) .lt. knots2(j)) exit
        ! First order divided differences
        dd(1) = dif1(x(i), j+1, nknots2, knots2)
        dd(2) = dif1(x(i), j+2, nknots2, knots2)
        dd(3) = dif1(x(i), j+3, nknots2, knots2)
        dd(4) = dif1(x(i), j+4, nknots2, knots2)
        ! Second order divided differences
        ff(1) = (dd(2) - dd(1)) / (knots2(j+2) - knots2(j))
        ff(2) = (dd(3) - dd(2)) / (knots2(j+3) - knots2(j+1))
        ff(3) = (dd(4) - dd(3)) / (knots2(j+4) - knots2(j+2))
        ! Computation of cubic B-splines
        temp = (ff(3) - ff(2)) / (knots2(j+4) - knots2(j+1))
        temp = temp - (ff(2) - ff(1)) / (knots2(j+3) - knots2(j))
        if (abs(temp) .lt. 1.d-10) temp = 0.d0
        B(i,j) = temp
     end do
  end do
end subroutine cubicBsplines_general

!=======================================================================
subroutine integrated_cubicBsplines_general(t0, nx, x, nknots, knots, IB)
!=======================================================================
! Input: x(nx): values where to compute the integrated cubic
!               B-splines basis ;
!        knots(nknots): possibly unequidistant (but distinct) knots
!
! Output: IB(nx, ndx+3)
!
! By P. Lambert (ULiege 2012) from de De Boor formula relating
!  cubic B-splines to truncated powers :
!    Bj(t) = (t(j+4) - t(j)) [t(j),...,t(j+4)] max(0,(.-t)^3)
!  with divided differences
!

  implicit none
  integer, intent(in) :: nx, nknots
  double precision, intent(in) :: t0, x(nx), knots(nknots)
  double precision, intent(out) :: IB(nx,nknots+2)

  integer :: nknots2, i, j
  double precision :: knots2(nknots+6), ff(3), dd(4)
  double precision :: eps1, eps2, integrated_dif1, x2(nx+1), IB0(nknots+2), temp

  nknots2 = nknots + 6
  eps1 = knots(2) - knots(1)
  eps2 = knots(nknots) - knots(nknots-1)
  knots2 = (/ knots(1)-3.d0*eps1, knots(1)-2.d0*eps1, knots(1)-eps1, &
              knots(1:nknots), &
              knots(nknots)+eps2, knots(nknots)+2.d0*eps2, knots(nknots)+3.d0*eps2 /)
  ! Put 2x3 extra knots on the left and on the right of the knots range

  x2 = (/ t0, x(:) /)
  IB(:,:) = 0.d0 ; IB0(:) = 0.d0
  do i = 1, nx + 1
     do j = 1, nknots2-4
        if (x2(i) .lt. knots2(j)) exit
        ! First order divided differences
        dd(1) = integrated_dif1(x2(i), j+1, nknots2, knots2)
        dd(2) = integrated_dif1(x2(i), j+2, nknots2, knots2)
        dd(3) = integrated_dif1(x2(i), j+3, nknots2, knots2)
        dd(4) = integrated_dif1(x2(i), j+4, nknots2, knots2)
        ! Second order divided differences
        ff(1) = (dd(2) - dd(1)) / (knots2(j+2) - knots2(j))
        ff(2) = (dd(3) - dd(2)) / (knots2(j+3) - knots2(j+1))
        ff(3) = (dd(4) - dd(3)) / (knots2(j+4) - knots2(j+2))
        ! Computation of integrated cubic B-splines
        temp = (ff(3) - ff(2)) / (knots2(j+4) - knots2(j+1))
        temp = temp - (ff(2) - ff(1)) / (knots2(j+3) - knots2(j))
        if (abs(temp) .lt. 1.d-10) temp = 0.d0
        if (i .eq. 1) then
           IB0(j) = temp
        else
           IB(i-1,j) = temp
        endif
     end do
  end do
  do i = 1, nx
     do j = 1, nknots2-4
        temp = IB(i,j) - IB0(j)
!        if (abs(temp) .lt. 1.d-10) temp = 0.d0
        IB(i,j) = temp
     end do
  end do
end subroutine integrated_cubicBsplines_general

!=======================================================================
subroutine D1_cubicBsplines_general(nx, x, nknots, knots, dB)
!=======================================================================
! Input: x(nx): values where to compute the 1st derivative of the cubic
!               B-splines basis ;
!        knots(nknots): possibly unequidistant (but distinct) knots
!
! Output: dB(nx, ndx+3)
!
! By P. Lambert (2012) from de De Boor formula relating cubic B-splines
!  to truncated powers :
!    Bj(t) = (t(j+4) - t(j)) [t(j),...,t(j+4)] max(0,(.-t)^3)
!  with divided differences
!

  implicit none
  integer, intent(in) :: nx, nknots
  double precision, intent(in) :: x(nx), knots(nknots)
  double precision, intent(out) :: dB(nx,nknots+2)

  integer :: nknots2, i, j
  double precision :: knots2(nknots+6), ff(3), dd(4)
  double precision :: eps1, eps2, D1_dif1, temp

  ! Put 2x3 extra knots on the left and on the right of the knots range
  nknots2 = nknots + 6
  eps1 = knots(2) - knots(1)
  eps2 = knots(nknots) - knots(nknots-1)
  knots2 = (/ knots(1)-3.d0*eps1, knots(1)-2.d0*eps1, knots(1)-eps1, &
              knots(1:nknots), &
              knots(nknots)+eps2, knots(nknots)+2.d0*eps2, knots(nknots)+3.d0*eps2 /)
  dB(:,:) = 0.d0
  do i = 1, nx
     do j = 1, nknots2-4
        if (x(i) .lt. knots2(j)) exit
        temp = 0.
        ! First order divided differences
        dd(1) = D1_dif1(x(i), j+1, nknots2, knots2)
        dd(2) = D1_dif1(x(i), j+2, nknots2, knots2)
        dd(3) = D1_dif1(x(i), j+3, nknots2, knots2)
        dd(4) = D1_dif1(x(i), j+4, nknots2, knots2)
        ! Second order divided differences
        ff(1) = (dd(2) - dd(1)) / (knots2(j+2) - knots2(j))
        ff(2) = (dd(3) - dd(2)) / (knots2(j+3) - knots2(j+1))
        ff(3) = (dd(4) - dd(3)) / (knots2(j+4) - knots2(j+2))
        ! Computation of cubic B-splines
        temp = (ff(3) - ff(2)) / (knots2(j+4) - knots2(j+1))
        temp = temp - (ff(2) - ff(1)) / (knots2(j+3) - knots2(j))
        if (abs(temp) .lt. 1.d-10) temp = 0.d0
        dB(i,j) = temp
     end do
  end do
end subroutine D1_cubicBsplines_general

!=======================================================================
subroutine D2_cubicBsplines_general(nx, x, nknots, knots, d2B)
!=======================================================================
! Input: x(nx): values where to compute the 1st derivative of the cubic
!               B-splines basis ;
!        knots(nknots): possibly unequidistant (but distinct) knots
!
! Output: dB(nx, ndx+3)
!
! By P. Lambert (2012) from de De Boor formula relating cubic B-splines
!  to truncated powers :
!    Bj(t) = (t(j+4) - t(j)) [t(j),...,t(j+4)] max(0,(.-t)^3)
!  with divided differences
!

  implicit none
  integer, intent(in) :: nx, nknots
  double precision, intent(in) :: x(nx), knots(nknots)
  double precision, intent(out) :: d2B(nx,nknots+2)

  integer :: nknots2, i, j
  double precision :: knots2(nknots+6), ff(3), dd(4)
  double precision :: eps1, eps2, D2_dif1, temp

  ! Put 2x3 extra knots on the left and on the right of the knots range
  nknots2 = nknots + 6
  eps1 = knots(2) - knots(1)
  eps2 = knots(nknots) - knots(nknots-1)
  knots2 = (/ knots(1)-3.d0*eps1, knots(1)-2.d0*eps1, knots(1)-eps1, &
              knots(1:nknots), &
              knots(nknots)+eps2, knots(nknots)+2.d0*eps2, knots(nknots)+3.d0*eps2 /)
  d2B(:,:) = 0.d0
  do i = 1, nx
     do j = 1, nknots2-4
        if (x(i) .lt. knots2(j)) exit
        temp = 0.
        ! First order divided differences
        dd(1) = D2_dif1(x(i), j+1, nknots2, knots2)
        dd(2) = D2_dif1(x(i), j+2, nknots2, knots2)
        dd(3) = D2_dif1(x(i), j+3, nknots2, knots2)
        dd(4) = D2_dif1(x(i), j+4, nknots2, knots2)
        ! Second order divided differences
        ff(1) = (dd(2) - dd(1)) / (knots2(j+2) - knots2(j))
        ff(2) = (dd(3) - dd(2)) / (knots2(j+3) - knots2(j+1))
        ff(3) = (dd(4) - dd(3)) / (knots2(j+4) - knots2(j+2))
        ! Computation of cubic B-splines
        temp = (ff(3) - ff(2)) / (knots2(j+4) - knots2(j+1))
        temp = temp - (ff(2) - ff(1)) / (knots2(j+3) - knots2(j))
        if (abs(temp) .lt. 1.d-10) temp = 0.d0
        d2B(i,j) = temp
     end do
  end do
end subroutine D2_cubicBsplines_general


!=======================================================================
subroutine cubicBsplines(x, nx, xl, xr, ndx, B)
!=======================================================================
! Input: x(nx): values where to compute the cubic B-splines basis ;
!        xl, xr, ndx: typically knots = seq(xl, xr, by=(xr-xl)/ndx)
!
! Output: B(nx, ndx+3)
!
! By P. Lambert (2010) inspired by Paul Eilers R code.
!

  implicit none
  integer nx, ndx
  double precision x(nx), xl, xr, B(nx,ndx+3)

  integer nknots, i, j
  double precision knots(ndx+7), dx, temp, cub

  nknots = ndx + 7
  dx = (xr - xl) / ndx

  knots(1) = xl - 3.*dx
  do i= 2, nknots
     knots(i) = knots(i-1) + dx
  end do

  do i = 1, nx
     do j = 1, nknots-4
        temp = 0.
        cub = x(i) - knots(j)
        if (cub .gt. 0.) then
           temp = temp + cub*cub*cub
           cub = x(i) - knots(j+1)
           if (cub .gt. 0.) then
              temp = temp - 4.*cub*cub*cub
              cub = x(i) - knots(j+2)
              if (cub .gt. 0.) then
                 temp = temp + 6.*cub*cub*cub
                 cub = x(i) - knots(j+3)
                 if (cub .gt. 0.) then
                    temp = temp - 4.*cub*cub*cub
                    cub = x(i) - knots(j+4)
                    if (cub .gt. 0.) then
                       temp = temp + cub*cub*cub
                    end if
                 end if
              end if
           end if
        end if
        B(i,j) = temp / (6.*dx*dx*dx)
        if (abs(B(i,j)) .lt. 1e-10) B(i,j)= 0.
     end do
  end do
end subroutine cubicBsplines

!=======================================================================
subroutine DcubicBsplines(x, nx, xl, xr, ndx, B)
!=======================================================================
! Input: x(nx): values where to compute the cubic B-splines basis ;
!        xl, xr, ndx: typically knots = seq(xl, xr, by=(xr-xl)/ndx)
!
! Output: 1st derivative of B(nx, ndx+3)
!
! By P. Lambert (2010) inspired by Paul Eilers R code.
!

  implicit none
  integer nx, ndx
  double precision x(nx), xl, xr, B(nx,ndx+3)

  integer nknots, i, j
  double precision knots(ndx+7), dx, temp, cub

  nknots = ndx + 7
  dx = (xr - xl) / ndx

  knots(1) = xl - 3.*dx
  do i= 2, nknots
     knots(i) = knots(i-1) + dx
  end do

  do i = 1, nx
     do j = 1, nknots-4
        temp = 0.
        cub = x(i) - knots(j)
        if (cub .gt. 0.) then
           temp = temp + cub**2
           cub = x(i) - knots(j+1)
           if (cub .gt. 0.) then
              temp = temp - 4.*cub**2
              cub = x(i) - knots(j+2)
              if (cub .gt. 0.) then
                 temp = temp + 6.*cub**2
                 cub = x(i) - knots(j+3)
                 if (cub .gt. 0.) then
                    temp = temp - 4.*cub**2
                    cub = x(i) - knots(j+4)
                    if (cub .gt. 0.) then
                       temp = temp + cub**2
                    end if
                 end if
              end if
           end if
        end if
        B(i,j) = 3. * temp / (6.*dx**3)
        if (abs(B(i,j)) .lt. 1e-10) B(i,j)= 0.
     end do
  end do
end subroutine DcubicBsplines

!=======================================================================
subroutine D2cubicBsplines(x, nx, xl, xr, ndx, B)
!=======================================================================
! Input: x(nx): values where to compute the cubic B-splines basis ;
!        xl, xr, ndx: typically knots = seq(xl, xr, by=(xr-xl)/ndx)
!
! Output: 1st derivative of B(nx, ndx+3)
!
! By P. Lambert (2010) inspired by Paul Eilers R code.
!

  implicit none
  integer nx, ndx
  double precision x(nx), xl, xr, B(nx,ndx+3)

  integer nknots, i, j
  double precision knots(ndx+7), dx, temp, cub

  nknots = ndx + 7
  dx = (xr - xl) / ndx

  knots(1) = xl - 3.*dx
  do i= 2, nknots
     knots(i) = knots(i-1) + dx
  end do

  do i = 1, nx
     do j = 1, nknots-4
        temp = 0.
        cub = x(i) - knots(j)
        if (cub .gt. 0.) then
           temp = temp + cub
           cub = x(i) - knots(j+1)
           if (cub .gt. 0.) then
              temp = temp - 4.*cub
              cub = x(i) - knots(j+2)
              if (cub .gt. 0.) then
                 temp = temp + 6.*cub
                 cub = x(i) - knots(j+3)
                 if (cub .gt. 0.) then
                    temp = temp - 4.*cub
                    cub = x(i) - knots(j+4)
                    if (cub .gt. 0.) then
                       temp = temp + cub
                    end if
                 end if
              end if
           end if
        end if
        B(i,j) = 6. * temp / (6.*dx**3)
        if (abs(B(i,j)) .lt. 1e-10) B(i,j)= 0.
     end do
  end do
end subroutine D2cubicBsplines

!=======================================================================
subroutine integratedcubicBsplines(t0, x, nx, xl, xr, ndx, IB)
!=======================================================================
! Input: t0: origin in the computation of  int_t0^x B(u) du
!        x(nx): values where to compute the integrated cubic B-splines basis ;
!        xl, xr, ndx: typically knots = seq(xl, xr, by=(xr-xl)/ndx)
!
! Output: IB(nx, ndx+3): :  int_t0^x B(u) du
!
! By P. Lambert (2011) inspired by Paul Eilers R code.
!

  implicit none
  ! input
  integer nx, ndx
  double precision t0, x(nx), xl, xr, IB(nx,ndx+3)

  ! within
  integer nknots, i, j, xmin_idx(1)
  double precision knots(ndx+7), dx, temp, cub, x2(nx+1), IB0(ndx+3)

  nknots = ndx + 7
  dx = (xr - xl) / ndx

  knots(1) = xl - 3.*dx
  do i= 2, nknots
     knots(i) = knots(i-1) + dx
  end do

  x2 = (/ t0, x(:) /)
  do i = 1, nx + 1
     do j = 1, nknots-4
        temp = 0.
        cub = x2(i) - knots(j)
        if (cub .gt. 0.) then
           temp = temp + cub**4
           cub = x2(i) - knots(j+1)
           if (cub .gt. 0.) then
              temp = temp - 4.*cub**4
              cub = x2(i) - knots(j+2)
              if (cub .gt. 0.) then
                 temp = temp + 6.*cub**4
                 cub = x2(i) - knots(j+3)
                 if (cub .gt. 0.) then
                    temp = temp - 4.*cub**4
                    cub = x2(i) - knots(j+4)
                    if (cub .gt. 0.) then
                       temp = temp + cub**4
                    end if
                 end if
              end if
           end if
        end if
        if (i .eq. 1) then
           IB0(j) = temp / (6.*dx**3) / 4.
        else
           IB(i-1,j) = temp / (6.*dx**3) / 4.
        endif
     end do
  end do
  do i = 1, nx
     IB(i,:) = IB(i,:) - IB0(:)
  end do
end subroutine IntegratedcubicBsplines

