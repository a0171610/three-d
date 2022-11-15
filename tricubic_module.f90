module tricubic_module
  ! tricubic interpolation
    implicit none
    private
  
    real(8), dimension(4, 4, 4), private :: c
    real(8), dimension(64, 64), private :: wt
    public :: tricof, trint, tricubic_init
  
  contains

    subroutine tricubic_init
      implicit none
      integer(8) :: i

      open(10, file="matrix.csv")
      do i = 1, 64
        read(10, *) wt(i, :)
      end do
      close(10)

    end subroutine tricubic_init

    subroutine tricof(f, fx, fy, fz, fxy, fxz, fyz, fxyz, dx, dy, dz)
      implicit none
  
      real(8), dimension(8), intent(in) :: f, fx, fy, fz, fxy, fxz, fyz, fxyz
      real(8), intent(in) :: dx, dy, dz
  
      integer(8) :: i, j, k, l
      real(8), dimension(64) :: x, cl
      real(8) :: dxdy
  
      dxdy = dx*dy
      do i = 1, 8
        x(i) = f(i)
        x(i + 8) = fx(i) * dx
        x(i + 16) = fy(i) * dy
        x(i + 24) = fz(i) * dz
        x(i + 32) = fxy(i) * dx * dy
        x(i + 40) = fxz(i) * dx * dz
        x(i + 48) = fyz(i) * dy * dz
        x(i + 56) = fxyz(i) * dx * dy * dz
      end do
      cl(:) = matmul(wt, x)
      l = 0
      do k = 1, 4
        do j = 1, 4
          do i = 1, 4
            l = l + 1
            c(i, j, k) = cl(l)
          enddo
        end do
      end do
  
    end subroutine tricof
  
    function trint(t, u, v) result(fi)
      implicit none
  
      real(8), intent(in) :: t, u, v
      real(8) :: fi
      real(8) :: a, b, d
      integer(8) :: i, j, k
  
      fi = 0.0d0
      do i = 1, 4
        do j = 1, 4
          do k = 1, 4
            a = t ** (i - 1)
            b = u ** (j - 1)
            d = v ** (k - 1)
            fi = fi + c(i,j,k) * a * b * d
          enddo
        enddo
      end do
    end function trint

end module tricubic_module
