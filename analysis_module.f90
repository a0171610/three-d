module analysis_module
  use grid_module, only: gphi, gphi_initial, lat, wgt, lon, ntrunc, Umax, nz
  use legendre_transform_module, only: legendre_analysis
  use math_module, only: math_pi
  use planet_module, only: planet_radius
  use time_module, only: deltat
  implicit none

contains

  subroutine error_log()
    implicit none
    integer(8) :: i, j, k, nlon, nlat
    real(8), allocatable :: w(:, :)
    real(8) :: dq, dqp, rmse
    real(8) :: sum_g1, sum_g2

    nlat = size(gphi, 2)
    nlon = size(gphi, 1)
    allocate(w(nlon, nlat))

    do j=1, nlat
      w(:, j) = wgt(j)
    end do

    open(10, file="log.txt")
    do i = 1, nlon
      do j = 1, nlat
    !    write(10,*) lon(i), lat(j), gphi(i, j)
      enddo
    enddo
    close(10)

    open(12, file="error.txt")
    do i = 1, nlon
        do j = 1, nlat
    !      write(12,*) lon(i), lat(j), gphi_initial(i, j) - gphi(i, j)
        end do
    end do
    close(12)

    open(14, file = "error_equator.txt")
    do i = 1, nlat
    !  write(14,*) lat(i), gphi(1, i) - gphi_initial(1, i)
    end do

    dq = 0.0d0
    dqp = 0.0d0
    do i = 1, nlon
        do j = 1, nlat
          do k = 1, nz
            dq = dqp + (gphi(i, j, k) - gphi_initial(i, j, k)) * wgt(j)
            dq = dq / dble(nlon * nlat)
            if(gphi(i, j, k) > gphi_initial(i, j, k)) then
              dqp = dqp + (gphi(i, j, k) - gphi_initial(i, j, k)) * wgt(j)
              dqp = dqp / dble(nlon * nlat)
            endif
          enddo
        end do
    end do

    rmse = 0.0d0
    do i = 1, nlon
      do j = 1, nlat
        do k = 1, nz
          rmse = rmse + (gphi(i, j, k) - gphi_initial(i, j, k)) ** 2
          rmse = rmse / dble(nlon * nlat)
        enddo
      end do
    end do
    rmse = sqrt(rmse)

    write(*,*) 'error = ', dq, 'positive error = ', dqp, 'RMSE = ', rmse

    sum_g1 = 0.0d0
    sum_g2 = 0.0d0

    do i = 1, nlon
      do j = 1, nlat
        do k = 1, nz
          sum_g1 = sum_g1 + gphi_initial(i, j, k) * w(i, j)
          sum_g2 = sum_g2 + gphi(i, j, k) * w(i, j)
        enddo
      end do
    end do
    write(*,*) "initial global mass sum", sum_g1, "final global mass sum", sum_g2

  ! l1ノルム
    sum_g1 = 0.0d0
    sum_g2 = 0.0d0

    do i = 1, nlon
      do j = 1, nlat
        do k = 1, nz
          sum_g1 = sum_g1 + ((gphi(i, j, k) - gphi_initial(i, j, k)) * wgt(j))
          sum_g2 = sum_g2 + (gphi_initial(i, j, k) * wgt(j))
        enddo
      enddo
    enddo

    write(*,*) "l1 norm = ", sum_g1 / sum_g2, sum_g1, sum_g2

    ! l2ノルムを求める
    sum_g1 = 0.0d0
    sum_g2 = 0.0d0

    do i = 1, nlon
      do j = 1, nlat
        do k = 1, nz
          sum_g1 = sum_g1 + ((gphi(i, j, k) - gphi_initial(i, j, k)) ** 2) * wgt(j)
          sum_g2 = sum_g2 + (gphi_initial(i, j, k) ** 2) * wgt(j)
        enddo
      enddo
    enddo

    write(*,*) "l2 norm = ", sqrt(sum_g1 / sum_g2), sum_g1, sum_g2

    write(*,*) 'l inf norm = ', maxval(gphi - gphi_initial)

    write(*,*) 'nlon = ', nlon

  end subroutine error_log
end module analysis_module
