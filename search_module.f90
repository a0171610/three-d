module search_module
  implicit none
  private

  public :: search_linear, search_bisection

contains
  function search_linear(y, x, i00) result(i)
    implicit none

    real(8), dimension(:), intent(in) :: y
    real(8) :: x
    integer(8), intent(in), optional :: i00

    integer(8) :: i

    integer(8) :: n, i0
    real(8) :: s

    n = size(y)
! ascending y: s =  1
! decending y: s = -1
    s = sign(1.0d0, y(n)-y(1))

    if (present(i00)) then
      i0 = i00
    else
      i0 = 1
    end if
! out of range
    if (s*(x-y(1))<0.0d0) then
      i = -1
      return
    else if (s*(x-y(n))>=0.0d0) then
      i = n + 1
      return
    end if
    i0 = max(1,min(n-1,i0))

    if (s*(x-y(i0+1))<0.0d0) then
      if (s*(x-y(i0))>=0.0d0) then
! ascending y: y(i0) <= x < y(i0+1)
! descending y: y(i0) > x >= y(i0+1)
        i = i0
        return
      else
! search with decreasing i
        do i=i0, 1, -1
          if (s*(x-y(i))>=0.0d0) then
            exit
          end if
        end do
      end if
    else
! search with increasing i
      do i=i0+1, n-1
        if (s*(x-y(i+1))<0.0d0) then
          exit
        end if
      end do
    end if

  end function search_linear

  function search_bisection(y, x) result(i)
    implicit none

    real(8), dimension(:), intent(in) :: y
    real(8) :: x

    integer(8) :: i

    integer(8) :: n, il, im, iu
    logical :: lascend

    n = size(y)
    lascend = y(n) > y(1)
    il = 0
    iu = n+1
    do
      if (iu-il==1) then
        exit
      end if
      im = (iu+il)/2
      if ((x > y(im)).eqv.lascend) then
        il = im
      else
        iu = im
      end if
    end do
    i = il

  end function search_bisection

end module search_module
