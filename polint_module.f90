module polint_module
  implicit none
  private

! polynomial interpolation
! Source: Numerical Recipies

! Author: Takeshi Enomoto
! History:
! 2007-11-13 translated from Pascal

  real(8), dimension(:), allocatable, private :: c, d, ytmp

  public :: polint, polin2

contains

  subroutine polint(xa, ya, x, y, dy)
    implicit none

    real(8), dimension(:), intent(in) :: xa, ya
    real(8), intent(in) :: x
    real(8), intent(out) :: y, dy

    integer(8) :: n, ns, i, m
    real(8) :: w, hp, ho, dift, dif, den

    n = size(xa)
! find index ns of the closest table entry 
    ns = 1
    dif = abs(x-xa(1))
    do i=1, n
      dift = abs(x-xa(i))
      if (dift < dif) then
        ns = i
        dif = dift
      end if
      c(i) = ya(i) ! initialization
      d(i) = ya(i)
    end do
    y = ya(ns) ! initial approximation
    ns = ns - 1
    do m=1, n-1 
      do i=1, n-m
        ho = xa(i) - x
        hp = xa(i+m) - x
        w = c(i+1) - d(i)
        den = ho - hp
        if (den==0) then
          print *, "error in POLINT"
          stop
        end if
        den = w/den
        d(i) = hp*den
        c(i) = ho*den
      end do
      if (2*ns<n-m) then
        dy = c(ns+1)
      else
        dy = d(ns)
        ns = ns - 1
      end if
      y = y + dy
    end do
      
  end subroutine polint

  subroutine polin2(x1a, x2a, ya, x1, x2, y, dy)
    implicit none

    real(8), dimension(:), intent(in) :: x1a, x2a
    real(8), dimension(:,:), intent(in) :: ya
    real(8), intent(in) :: x1, x2
    real(8), intent(out) :: y, dy

    integer(8) :: j, n

    n = size(x2a)
    do j=1, n
      call polint(x1a, ya(:,j), x1, ytmp(j), dy)
    end do
    call polint(x2a, ytmp, x2, y, dy)

  end subroutine polin2

end module polint_module