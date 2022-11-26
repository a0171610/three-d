module uv_hadley_module
  use math_module, only: pi=>math_pi
  implicit none

  real(8), parameter, private :: tau = 864000.0d0, K = 5.0d0, u0 = 40.0d0, w0 = 0.15d0
  real(8), parameter, private :: Rd = 287.0d-3, T0 = 300.0d0, ztop = 12000.0d0
contains

subroutine uv_hadley(t, lon, lat, height, gu, gv, gw)
  implicit none

  real(8), intent(in) :: t
  real(8), dimension(:), intent(in) :: lon, lat, height
  real(8), dimension(:, :, :), intent(inout) :: gu, gv, gw

  integer(8) :: i, j, k1, nx, ny, nz

  nx = size(lon)
  ny = size(lat)
  nz = size(height)

  do i = 1, nx
    do j = 1, ny
      do k1 = 1, nz
        gu(i, j, k1) = calc_u(lat(j))
        gv(i, j, k1) = calc_v(lat(j), height(k1), t)
        gw(i, j, k1) = calc_w(lat(j), height(k1), t)
      end do
    end do
  end do

end subroutine uv_hadley

  function calc_u(lat) result(ans)
    use planet_module, only: planet_radius
    implicit none
    real(8), intent(in) :: lat
    real(8) :: ans

    ans = u0 * cos(lat) / planet_radius
  end function calc_u

  function calc_v(lat, h, t) result(ans)
    use planet_module, only: transorm_height_to_pressure
    implicit none
    real(8), intent(in) :: lat, h, t
    real(8) :: ans, rho0
    real(8) :: pres

    pres = transorm_height_to_pressure(h)
    rho0 = calc_rho(1000.0d0)

    ans = -(w0 * pi * rho0) / (K * ztop * calc_rho(pres)) * cos(lat) * sin(K * lat)
    ans = ans * cos(pi * h / ztop) * cos(pi * t / tau)

  end function calc_v

  function calc_w(lat, h, t) result(ans)
    use planet_module, only: transorm_height_to_pressure
    implicit none
    real(8), intent(in) :: lat, h, t
    real(8) :: ans, rho0, tmp
    real(8), parameter :: g = 9.81606d0
    real(8) :: pres

    rho0 = calc_rho(1000.0d0)
    pres = transorm_height_to_pressure(h)

    ans = (w0 * rho0 / (K * calc_rho(pres))) * sin(pi * h / ztop) * cos(pi * t / tau)
    tmp = K * cos(lat) * cos(K * lat) - 2.0d0 * sin(K * lat) * sin(lat)
    
    ans = ans * tmp
    ans = ans
  end function calc_w

  function calc_rho(pres) result(ans)
    implicit none
    real(8), intent(in) :: pres
    real(8) :: ans

    ans = pres / (Rd * T0)
  end function calc_rho


end module uv_hadley_module
