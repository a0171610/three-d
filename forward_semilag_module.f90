module forward_semilag_module

  use grid_module, only: nlon, nlat, nz, ntrunc, gphi, height, &
    gu, gv, gw, gphi, sphi_old, sphi, latitudes=>lat, lon, coslatr, wgt
  private
  
  real(8), dimension(:, :, :), allocatable, private :: deplon, deplat, deph
  real(8), dimension(:, :, :), allocatable, private :: gphi_old, gphi_initial
  real(8), dimension(:, :, :), allocatable, private :: intersect_lon, intersect_lat, intersect

  private :: update
  public :: semilag_forward_init, semilag_forward_timeint, semilag_forward_clean

contains

  subroutine semilag_forward_init()
    use interpolate3d_module, only: interpolate_init
    use legendre_transform_module, only: legendre_synthesis
    use cascade_interpolation_module, only: cascade_interpolate_init
    implicit none

    integer(8) :: i,j

    allocate(gphi_old(nlon, nlat, nz), deplon(nlon,nlat,nz),deplat(nlon,nlat,nz), &
      gphi_initial(nlon, nlat, nz), deph(nlon, nlat, nz))
    allocate(intersect_lon(nlon, nlat, nz), intersect_lat(nlon, nlat, nz), intersect(nlon, nlat, nz))
    call interpolate_init(gphi)
    call cascade_interpolate_init

    gphi_old(:, :, :) = gphi(:, :, :)
    gphi_initial(:, :, :) = gphi(:, :, :)

    open(11, file="animation.txt")
    do i = 1, nlat
      do j = 1, nz
        write(11,*) latitudes(i), height(j), gphi(nlon/2, i, j)
      end do
    end do

  end subroutine semilag_forward_init

  subroutine semilag_forward_clean()
    use interpolate3d_module, only: interpolate_clean
    implicit none

    deallocate(gphi_old, deplon, deplat)
    call interpolate_clean()

  end subroutine semilag_forward_clean

  subroutine semilag_forward_timeint()
    use math_module, only: pi=>math_pi
    use time_module, only: nstep, deltat, hstep
    implicit none

    integer(8) :: i, j, k

    do i = 1, nstep
      call update((i-0.5d0)*deltat, deltat)
      write(*, *) "step=", i, "maxval = ", maxval(gphi), 'minval = ', minval(gphi)
      if ( mod(i, hstep) == 0 ) then
        do j = 1, nlat
          do k = 1, nz
            write(11,*) latitudes(j), height(k), gphi(nlon/2, j, k)
          end do
        end do
      endif
    end do
    close(11)
    open(10, file="log.txt")
    do i = 1, nlat
      do j = 1, nz
        write(10,*) latitudes(i), height(j), gphi(nlon/2, i, j)
      enddo
    enddo
    close(10)
  end subroutine semilag_forward_timeint

  subroutine update(t, dt)
    use math_module, only: &
      pi=>math_pi, pir=>math_pir, pih=>math_pih
    use upstream_forward_module, only: find_forward_points
    use time_module, only: case
    use uv_module, only: uv_div
    use uv_hadley_module, only: uv_hadley
    use spline_interpolate_module, only: interpolate_spline_1index
    use search_module, only: search_bisection
    use cascade_interpolation_module, only: cascade_interpolate
    implicit none

    real(8), intent(in) :: t, dt
    real(8) :: target_h, h1, h2, h_dist(nz), h_loc(nz), f(nz)

    integer(8) :: i, j, k, m, n, id

    select case(case)
      case('hadley')
        call uv_hadley(t, lon, latitudes, height, gu, gv, gw)
      case('div')
        call uv_div(t, lon, latitudes, height, gu, gv, gw)
      case default
        print *, "No matching initial field"
      stop
    end select

    call find_forward_points(t, dt, deplon, deplat, deph)

    do k = 1, nz
      target_h = height(k)
      do i = 1, nlon
        do j = 1, nlat
          id = search_bisection(deph(i, j, :), target_h)
          if (id == 0) then
            intersect_lon(i, j, k) = deplon(i, j, 1)
            intersect_lat(i, j, k) = deplat(i, j, 1)
          else
            h1 = deph(i, j, id)
            h2 = deph(i, j, id + 1)
            intersect_lon(i, j, k) = deplon(i, j, id) + ((target_h - h1) / (h2 - h1)) * (deplon(i, j, id+1) - deplon(i, j, id))
            intersect_lat(i, j, k) = deplat(i, j, id) + ((target_h - h1) / (h2 - h1)) * (deplat(i, j, id+1) - deplat(i, j, id))
          endif
        end do
      end do
    end do

    do i = 1, nlon
      do j = 1, nlat
        h_dist(:) = 0.0d0
        do k = 2, nz
          h_dist(k) = h_dist(k - 1) + cartesian_dist(deplon(i,j,k),deplat(i,j,k),deph(i,j,k), &
            deplon(i,j,k-1),deplat(i,j,k-1),deph(i,j,k-1))
        end do
        do k = 1, nz
          id = search_bisection(deph(i, j, :), height(k))
          if (id == 0) then
            h_loc(k) = h_dist(1) - cartesian_dist(intersect_lon(i,j,k),intersect_lat(i,j,k),height(k), &
              deplon(i,j,1),deplat(i,j,1),deph(i,j,1))
          else
            h_loc(k) = h_dist(id) + cartesian_dist(intersect_lon(i,j,k),intersect_lat(i,j,k),height(k), &
              deplon(i,j,id),deplat(i,j,id),deph(i,j,id))
          endif
        end do
        f(:) = gphi_old(i, j, :)
        call interpolate_spline_1index(nz, h_dist, f, h_loc, intersect(i, j, :))
      end do
    end do

    do k = 1, nz
      call cascade_interpolate(deplon(:, :, k), deplat(:, :, k), intersect(:, :, k), gphi(:, :, k))
    end do

    gphi_old = gphi
  end subroutine update

  function cartesian_dist(x1, y1, z1, x2, y2, z2) result(l)
    implicit none
    real(8), intent(in) :: x1, y1, z1, x2, y2, z2
    real(8) :: l

    l = sqrt((x1 - x2)**2 + (y1 - y2)**2 + (z1 - z2)**2)
  end function cartesian_dist

end module forward_semilag_module