module uv_hadley_module
  use math_module, only: pi=>math_pi
  implicit none

  real(8), parameter, private :: tau = 864000.0d0, K = 5.0d0, u0 = 40.0d0, w0 = 0.15d0
  real(8), parameter, private :: Rd = 287.0d-3, T0 = 300.0d0, ztop = 12000.0d0
contains

subroutine uv_hadley(t, lon, lat, pres, gu, gv, gomega)
  implicit none

  real(8), intent(in) :: t
  real(8), dimension(:), intent(in) :: lon, lat, pres
  real(8), dimension(:, :, :), intent(inout) :: gu, gv, gomega

  integer(8) :: i, j, k1, nx, ny, nz

  nx = size(lon)
  ny = size(lat)
  nz = size(pres)

  do i = 1, nx
    do j = 1, ny
      do k1 = 1, nz
        gu(i, j, k1) = calc_u(lat(j))
        gv(i, j, k1) = calc_v(lat(j), pres(k1), t)
        gomega(i, j, k1) = calc_omega(lat(j), pres(k1), t)
      end do
    end do
  end do

end subroutine uv_hadley

  function calc_u(lat) result(ans)
    implicit none
    real(8), intent(in) :: lat
    real(8) :: ans

    ans = u0 * cos(lat)
  end function calc_u

  function calc_v(lat, pres, t) result(ans)
    use planet_module, only: transorm_pressure_to_height
    implicit none
    real(8), intent(in) :: lat, pres, t
    real(8) :: ans, rho0

    rho0 = calc_rho(1000.0d0)

    ans = -(w0 * pi * rho0) / (K * ztop * calc_rho(pres)) * cos(lat) * sin(K * lat)
    ans = ans * cos(pi * transorm_pressure_to_height(pres) / ztop) * cos(pi * t / tau)

  end function calc_v

  function calc_omega(lat, pres, t) result(ans)
    use planet_module, only: transorm_pressure_to_height
    implicit none
    real(8), intent(in) :: lat, pres, t
    real(8) :: ans, rho0, tmp
    real(8), parameter :: g = 9.81606d0

    rho0 = calc_rho(1000.0d0)

    ans = -(g * w0 * rho0 / K) * sin(pi * transorm_pressure_to_height(pres) / ztop) * cos(pi * t / tau)
    tmp = K * cos(lat) * cos(K * lat) - 2.0d0 * sin(K * lat) * sin(lat)
    
    ans = ans * tmp
  end function calc_omega

  function calc_rho(pres) result(ans)
    implicit none
    real(8), intent(in) :: pres
    real(8) :: ans

    ans = pres / (Rd * T0)
  end function calc_rho


end module uv_hadley_module