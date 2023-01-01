module upstream_forward_module
  use grid_module, only: longitudes=>lon, latitudes=>lat, height
  use uv_hadley_module, only: hadley_u=>calc_u, hadley_v=>calc_v, hadley_w=>calc_w
  use uv_module, only: calc_ua, calc_ud, calc_va, calc_w
  use sphere_module, only: lonlat2xyz, uv2xyz
  use math_module, only: pi2=>math_pi2
  use time_module, only: case
  implicit none

contains

  subroutine find_forward_points(t, dt, deplon, deplat, deph)
    implicit none
    real(8), intent(in) :: t, dt
    real(8), dimension(:, :, :), intent(out) :: deplon, deplat, deph

    integer(8) :: i, j, k
    integer(8) :: nx, ny, nz
    real(8) :: lon, lat, h
    real(8) :: xg, yg, zg
    real(8) :: midx, midy, midz, midh
    real(8) :: k1_x, k1_y, k1_z, k1_h, k2_x, k2_y, k2_z, k2_h, k3_x, k3_y, k3_z, k3_h
    real(8) :: k1_u, k1_v, k1_w, k2_u, k2_v, k2_w, k3_u, k3_v, k3_w
    real(8) :: k1_xdot, k1_ydot, k1_zdot
    real(8) :: k2_xdot, k2_ydot, k2_zdot
    real(8) :: k3_xdot, k3_ydot, k3_zdot
    real(8) :: lon_tmp, lat_tmp
    real(8) :: x_tmp, y_tmp, z_tmp, h_tmp

    nx = size(deplon, 1)
    ny = size(deplon, 2)
    nz = size(deplon, 3)

    do i = 1, nx
      do j = 1, ny
        do k = 1, nz
          lon = longitudes(i)
          lat = latitudes(j)
          h = height(k)

          call lonlat2xyz(lon, lat, xg, yg, zg)

          select case(case)
          case ('hadley')
            k1_u = hadley_u(lat)
            k1_v = hadley_v(lat, h, t)
            k1_w = hadley_w(lat, h, t)
          case ('div')
            k1_u = calc_ua(lon, lat,t) + calc_ud(lon, lat, h, t)
            k1_v = calc_va(lon, lat, t)
            k1_w = calc_w(lon, lat, h, t)
          case default
            print *, 'no matchig velocity'
          end select

          call uv2xyz(k1_u, k1_v, lon, lat, k1_xdot, k1_ydot, k1_zdot)
          k1_x = dt * k1_xdot
          k1_y = dt * k1_ydot
          k1_z = dt * k1_zdot
          k1_h = dt * k1_w

          x_tmp = xg + k1_x / 3.0d0
          y_tmp = yg + k1_y / 3.0d0
          z_tmp = zg + k1_z / 3.0d0
          lon_tmp = modulo(atan2(y_tmp, x_tmp)+pi2, pi2)
          lat_tmp = asin(z_tmp / sqrt(x_tmp**2 + y_tmp**2 + z_tmp**2))
          h_tmp = h + k1_h / 3.0d0
          select case(case)
          case('hadley')
            k2_u = hadley_u(lat_tmp)
            k2_v = hadley_v(lat_tmp, h_tmp, t + dt / 3.0d0)
            k2_w = hadley_w(lat_tmp, h_tmp, t + dt / 3.0d0)
          case('div')
            k2_u = calc_ua(lon_tmp, lat_tmp,t + dt/3.0d0) + calc_ud(lon_tmp, lat_tmp, h_tmp, t+dt/3.0d0)
            k2_v = calc_va(lon_tmp, lat_tmp, t + dt/3.0d0)
            k2_w = calc_w(lon_tmp, lat_tmp, h_tmp, t + dt/3.0d0)
          case default
            print *, 'no matchig velocity'
          end select
          call uv2xyz(k2_u, k2_v, lon_tmp, lat_tmp, k2_xdot, k2_ydot, k2_zdot)
          k2_x = dt * k2_xdot
          k2_y = dt * k2_ydot
          k2_z = dt * k2_zdot
          k2_h = dt * k2_w

          x_tmp = xg + k2_x * 2.0d0 / 3.0d0
          y_tmp = yg + k2_y * 2.0d0 / 3.0d0
          z_tmp = zg + k2_z * 2.0d0 / 3.0d0
          lon_tmp = modulo(atan2(y_tmp, x_tmp)+pi2, pi2)
          lat_tmp = asin(z_tmp / sqrt(x_tmp**2 + y_tmp**2 + z_tmp**2))
          h_tmp = h + k2_h * 2.0d0 / 3.0d0
          select case(case)
          case('hadley')
            k3_u = hadley_u(lat_tmp)
            k3_v = hadley_v(lat_tmp, h_tmp, t + dt*2.0d0 / 3.0d0)
            k3_w = hadley_w(lat_tmp, h_tmp, t + dt*2.0d0 / 3.0d0)
          case('div')
            k3_u = calc_ua(lon_tmp, lat_tmp,t + dt*2.0d0/3.0d0) + calc_ud(lon_tmp, lat_tmp, h_tmp, t+dt*2.d0/3.0d0)
            k3_v = calc_va(lon_tmp, lat_tmp, t + dt*2.0d0/3.0d0)
            k3_w = calc_w(lon_tmp, lat_tmp, h_tmp, t + dt*2.0d0/3.0d0)
          case default
            print *, 'no matchig velocity'
          end select
          call uv2xyz(k3_u, k3_v, lon_tmp, lat_tmp, k3_xdot, k3_ydot, k3_zdot)
          k3_x = dt * k3_xdot
          k3_y = dt * k3_ydot
          k3_z = dt * k3_zdot
          k3_h = dt * k3_w

          midx = xg + k1_x * 0.25d0 + k3_x * 0.75d0
          midy = yg + k1_y * 0.25d0 + k3_y * 0.75d0
          midz = zg + k1_z * 0.25d0 + k3_z * 0.75d0
          midh = h + k1_h * 0.25d0 + k3_h * 0.75d0

          deplon(i, j, k) = modulo(atan2(midy, midx)+pi2,pi2)
          deplat(i, j, k) = asin(midz / sqrt(midx**2 + midy**2 + midz**2))
          deph(i, j, k) = midh
        end do
      end do
    end do

  end subroutine find_forward_points
end module upstream_forward_module