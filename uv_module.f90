module uv_module
  use math_module, only: pi=>math_pi, pi2=>math_pi2
  use planet_module, only: d=>day_in_sec
  implicit none
  private

  real(8), parameter :: t1 = 5.0d0, kappa = 1.0d0, pt = 254.944d0, p0 = 1000.0d0, b = 0.2d0, tau = 12.0d0*d
  real(8), parameter :: Rd = 287.0d0, T0 = 300.0d0

  public :: uv_div, calc_w, calc_ua, calc_ud, calc_va

contains

  subroutine uv_div(t, lon, lat, height, gu, gv, gw)
    implicit none

    real(8), intent(in) :: t
    real(8), dimension(:), intent(in) :: lon, lat, height
    real(8), dimension(:, :, :), intent(inout) :: gu, gv, gw
    real(8), allocatable :: gua(:, :, :), gva(:, :, :), gud(:, :, :)

    integer(8) :: i, j, k, nx, ny, nz

    nx = size(lon)
    ny = size(lat)
    nz = size(height)

    allocate(gua(nx, ny, nz), gva(nx, ny, nz), gud(nx, ny, nz))

    do i = 1, nx
      do j = 1, ny
        do k = 1, nz
          gua(i, j, k) = calc_ua(lon(i), lat(j), t)
          gva(i, j, k) = calc_va(lon(i), lat(j), t)
          gud(i, j, k) = calc_ud(lon(i), lat(j), height(k), t)
          gw(i, j, k) = calc_w(lon(i), lat(j), height(k), t)
        end do
      end do
    end do

    gu = gua + gud
    gv = gva

  end subroutine uv_div

  function calc_ua(lon, lat, t) result(ans)
    implicit none
    real(8), intent(in) :: lon, lat, t
    real(8) :: ans
    real(8) :: lambda1

    lambda1 = lon - 2.0d0 * pi * t / tau
    ans = 10.0d0  * sin(lambda1)**2 * sin(2.0d0*lat) * cos(pi*t/tau) + 2.0d0*pi*cos(lat)
    ans = ans / tau
  end function calc_ua

  function calc_ud(lon, lat, h, t) result(ans)
    use planet_module, only: transorm_height_to_pressure
    implicit none
    real(8), intent(in) :: lon, lat, h, t
    real(8) :: ans
    real(8) :: lambda1, tmp, omega1
    real(8) :: p

    p = transorm_height_to_pressure(h)

    lambda1 = lon - 2.0d0 * pi * t / tau
    tmp = exp((pt - p) / (b * pt)) - exp((p - p0) / (b * pt))
    omega1 = 23000.0d0 * pi / tau
    ans = omega1 * cos(lambda1) * cos(lat)**2 * cos(2.0d0*pi*t/tau) * tmp
    ans = ans / (b * pt * 100.0d0)

  end function calc_ud

  function calc_va(lon, lat, t) result(ans)
    implicit none
    real(8), intent(in) :: lon, lat, t
    real(8) :: ans
    real(8) :: lambda1

    lambda1 = lon - 2.0d0 * pi * t / tau

    ans = 10.0d0 * sin(2.0d0*lambda1) * cos(lat) * cos(pi * t / tau)
    ans = ans / tau
  end function calc_va

  function calc_w(lon, lat, z, t) result(ans)
    use planet_module, only: transorm_height_to_pressure
    real(8), intent(in) :: lon, lat, z, t
    real(8) :: ans
    real(8) :: omega1, lambda1
    real(8) :: pres, rho

    pres = transorm_height_to_pressure(z)
    rho = pres / (Rd * T0)

    omega1 = 23000.0d0 * pi / tau
    lambda1 = lon - 2.0d0 * pi * t / tau

    ans = 1.0d0 + exp((pt - p0) / (b*pt)) - exp((pres-p0) / (b*pt)) - exp((pt-pres) / (b*pt))
    ans = omega1 * sin(lambda1) * cos(lat) * cos(2.0d0*pi*t/tau) * ans

    ans = -ans / (9.80616d0 * rho)
    ans = ans / 100.0d0
  end function calc_w

end module uv_module
