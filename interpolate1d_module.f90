module interpolate1d_module
  use grid_module, only: dh, height
  implicit none
  
  real(8), allocatable, private :: ff(:)
contains

  subroutine interpolate1d_init(f)
    real(8), intent(in) :: f(:)
    integer(8) :: sz

    sz = size(f, 1)
    allocate(ff(sz))
  end subroutine interpolate1d_init

  subroutine interpolate1d_set(f)
    real(8), intent(in) :: f(:)

    ff(:) = f(:)

  end subroutine interpolate1d_set

  subroutine interpolate1d_linear(h, ans)
    real(8), intent(in) :: h
    real(8), intent(out) :: ans
    integer(8) :: id
    real(8) :: ratio
    real(8), parameter :: eps = 1.0d-7

    id = int(h / dh) + 1
    ratio = (h - height(id)) / dh
    !write(*,*) h, dh, id
    ans = ratio * ff(id + 1) + (1.0d0 - ratio) * ff(id)

  end subroutine interpolate1d_linear
end module interpolate1d_module

