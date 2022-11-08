module mass_module
  implicit none
  private

  public :: mass_correct, local_mass_record, local_mass_correct

contains

  subroutine mass_correct(f, fold, gmax, gmin, w)
    implicit none

    real(8), dimension(:, :), intent(inout) :: f
    real(8), dimension(:, :), intent(in) :: fold, gmax, gmin, w

    real(8), dimension(size(f,1),size(f,2)) :: f1, f2, f3
    real(8) :: df, a, b, s2, s3, s2r, f1f2w, f1f1f2f2w, s3s2f3, b1, b2, f1f2

    integer(8) :: nx, ny, i, j

    nx = size(f, 1)
    ny = size(f, 2)
    df = 0.0d0
    s2 = 0.0d0
    s3 = 0.0d0
    do j = 1, ny
      do i = 1, nx
        df = df + (fold(i,j)-f(i,j)) * w(i,j)
        f1(i,j) = gmax(i,j) - f(i,j)
        f2(i,j) = gmin(i,j) - f(i,j)
        f3(i,j) = 0.5d0 * (gmax(i,j)+gmin(i,j)) - f(i,j)
        f1f2w = f1(i,j) * f2(i,j) * w(i,j)
        s2 = s2 + f1f2w
        s3 = s3 + f1f2w * f3(i,j)
      end do
    end do

! Calculate a and b
    s2r = 1.0d0 / s2
    b1 = 0.0d0
    b2 = 0.0d0
    do j = 1, ny
      do i = 1, nx
        f1f1f2f2w = f1(i,j) * f1(i,j) * f2(i,j) * f2(i,j) * w(i,j)
        s3s2f3 = s3*s2r - f3(i,j)
        b1 = b1 + df * s2r * s3s2f3 * f1f1f2f2w
        b2 = b2 + s3s2f3 * s3s2f3 * f1f1f2f2w
      end do
    end do
    b = b1 / b2
    a = (df - b*s3) * s2r

! Correct
    do j = 1, ny
      do i = 1, nx
        f1f2 = f1(i,j) * f2(i,j)
        f(i,j) = f(i,j) + a * f1f2 + b * f1f2 * f3(i,j)
      end do
    end do
    
  end subroutine mass_correct

  subroutine local_mass_record(deplon, deplat, A, B, C, D, w_record)
    use grid_module, only: nlon, nlat, pole_regrid
    use interpolate_module, only: find_stencil_
    implicit none

    integer(8) :: i, j, k
    real(8), intent(out) :: w_record(nlon, nlat)
    real(8), intent(in) :: deplon(nlon, nlat), deplat(nlon, nlat)
    real(8), intent(in) :: A(nlon, nlat), B(nlon, nlat), C(nlon, nlat), D(nlon, nlat)
    integer(8), dimension(4) :: is, js

    w_record(:, :) = 0.0d0
    do i = 1, nlon
      do j = 1, nlat
        call find_stencil_(deplon(i, j), deplat(i, j), is, js)
        do k = 1, 4
          call pole_regrid(is(k), js(k))
        end do
        w_record(is(1), js(1)) = w_record(is(1), js(1)) + A(i, j)
        w_record(is(2), js(2)) = w_record(is(2), js(2)) + B(i, j)
        w_record(is(3), js(3)) = w_record(is(3), js(3)) + C(i, j)
        w_record(is(4), js(4)) = w_record(is(4), js(4)) + D(i, j)
      end do
    end do

  end subroutine local_mass_record

  subroutine local_mass_correct(deplon, deplat, A, B, C, D, w_record)
    use grid_module, only: nlon, nlat, pole_regrid, wgt
    use interpolate_module, only: find_stencil_
    implicit none

    integer(8) :: i, j, k
    real(8), intent(in) :: deplon(nlon, nlat), deplat(nlon, nlat), w_record(nlon, nlat)
    real(8), intent(inout) :: A(nlon, nlat), B(nlon, nlat), C(nlon, nlat), D(nlon, nlat)
    integer(8), dimension(4) :: is, js

    do i = 1, nlon
      do j = 1, nlat
        call find_stencil_(deplon(i, j), deplat(i, j), is, js)
        do k = 1, 4
          call pole_regrid(is(k), js(k))
        end do
        A(i, j) = A(i, j) * wgt(js(1)) / (w_record(is(1), js(1)) * wgt(j))
        B(i, j) = B(i, j) * wgt(js(2)) / (w_record(is(2), js(2)) * wgt(j))
        C(i, j) = C(i, j) * wgt(js(3)) / (w_record(is(3), js(3)) * wgt(j))
        D(i, j) = D(i, j) * wgt(js(4)) / (w_record(is(4), js(4)) * wgt(j))
      end do
    end do
  end subroutine local_mass_correct

end module
