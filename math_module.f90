module math_module
  implicit none

! Mathematical constants
  real(8), parameter, public :: &
    math_pi = acos(-1.0d0), &
    math_pir = 1.0d0/math_pi, &
    math_pi2 = 2.0d0*math_pi, math_pih = 0.5d0*math_pi, &
    math_deg2rad = math_pi/180.0d0

contains

  function math_atan2(y,x) result(theta)
    implicit none
  
    real(8), intent(in) :: x, y
    real(8) :: theta
  
    if ((x /= 0.0d0) .and. (y /= 0.0d0)) then
          theta = atan2(y,x)
          if (theta < 0) then
            theta = theta + math_pi * 2.0d0
          end if
    else if (x == 0.0d0) then
        if (y > 0.0d0) then
            theta = math_pi * 0.5d0
        else if (y < 0.0d0) then
            theta = math_pi + 0.5d0 * math_pi
        else
            theta = 0.0d0
        end if
    else if (y == 0.0d0) then
        if (x >= 0.0d0) then
            theta = 0.0d0
        else
            theta = math_pi
        end if 
    end if
  
  end function math_atan2

  function sphere_orthodrome(lon1, lat1, lon2, lat2) result(l)
    ! returns the shortest distance between two points on the unit sphere
      implicit none
  
      real(8), intent(in) :: lon1, lat1, lon2, lat2
      real(8) :: l
  
      l = acos(sphere_cosine(sin(lat1),sin(lat2),cos(lat1),cos(lat2),lon1-lon2))
  
    end function sphere_orthodrome

    ! spherical triagle composed of points pa, pb, pc and arcs a, b, c
  function sphere_cosine(cosb, cosc, sinb, sinc, pa) result(cosa)
    implicit none

    real(8), intent(in) :: cosb, cosc, sinb, sinc, pa
    real(8) :: cosa

    cosa = cosb*cosc+sinb*sinc*cos(pa)

  end function sphere_cosine
end module math_module