module planet_module
  use math_module, only: pi=>math_pi
  implicit none

! Geophysical constants
  !real(8), public :: planet_radius = 6.371d6, day_in_sec = 86400.0d0, angular_velocity
  real(8), public, parameter :: planet_radius = 1.0d0, day_in_sec = 1.0d0
  real(8), public :: angular_velocity
  real(8), private, parameter :: p0 = 1000.0d0, T0 = 300.0d0, Rd = 287.0d0, g = 9.80616d0

  public :: planet_init, transorm_height_to_pressure, transorm_pressure_to_height

contains

  subroutine planet_init()
    implicit none

    angular_velocity= 2.0d0*pi/day_in_sec

  end subroutine planet_init

  function transorm_height_to_pressure(z) result(ans)
    implicit none
    real(8), intent(in) :: z
    real(8) :: ans
    ans = p0 * exp((-g * z) / (Rd * T0))

  end function transorm_height_to_pressure

  function transorm_pressure_to_height(pres) result(ans)
    real(8) :: pres, ans

    ans = (Rd * T0 / g) * log(p0 / pres)
  end function transorm_pressure_to_height

end module planet_module