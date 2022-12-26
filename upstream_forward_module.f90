module upstream_forward_module
  use grid_module, only: longitudes=>lon, latitudes=>lat, height
  use uv_hadley_module, only: calc_u, calc_v, calc_w
  use sphere_module, only: lonlat2xyz, uv2xyz
  use math_module, only: pi2=>math_pi2
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

          k1_u = calc_u(lat)
          k1_v = calc_v(lat, h, t)
          k1_w = calc_w(lat, h, t)
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
          k2_u = calc_u(lon_tmp)
          k2_v = calc_v(t + dt/3.0d0, lon_tmp, lat_tmp)
          k2_w = calc_w(lat_tmp, h_tmp, t + dt / 3.0d0)
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
          k3_u = calc_u(lon_tmp)
          k3_v = calc_v(lat_tmp, h_tmp, t + 2.0d0*dt/3.0d0)
          k3_w = calc_w(lat_tmp, h_tmp, t + 2.0d0*dt/3.0d0)
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