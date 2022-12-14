module sphere_module
  ! utility for a spherical topology
    implicit none
    private
  
    public :: lonlat2xyz, uv2xyz, xyz2uv, lon2i, lat2j, orthodrome
  
  contains
  
    subroutine lonlat2xyz(lon, lat, x, y, z)
      implicit none
  
      real(8), intent(in) :: lon, lat
      real(8), intent(out) :: x, y, z
  
      x = cos(lon)*cos(lat)
      y = sin(lon)*cos(lat)
      z = sin(lat)
  
    end subroutine lonlat2xyz
  
    subroutine uv2xyz(u, v, lon, lat, xd, yd, zd)
      implicit none
  
      real(8), intent(in) :: u, v, lon, lat
      real(8), intent(out) :: xd, yd, zd
  
      xd = -u*sin(lon) - v*cos(lon)*sin(lat)
      yd =  u*cos(lon) - v*sin(lon)*sin(lat)
      zd =  v*cos(lat)
  
    end subroutine uv2xyz
  
    subroutine xyz2uv(xd, yd, zd, lon, lat, u, v)   ! Ritchie1987 式(49)
      use math_module, only: pih=>math_pih  
      implicit none
  
      real(8), intent(in) :: xd, yd, zd, lon, lat
      real(8), intent(out) :: u, v
  
      if (abs(lat) == pih) then ! lon = 0 is chosen
        u =  yd
        v = -xd ! omitted division by sin(pi/2) = 1
      else
        u = cos(lon)*yd - sin(lon)*xd
        v = zd/cos(lat)
      end if
  
    end subroutine xyz2uv
  
    function lon2i(lon,nx) result(i)
    ! returns the closest longitudinal point i, not exceeding lon
    ! 1 <= return value <= nx
      use math_module, only: pi2=>math_pi2
      implicit none
  
      real(8), intent(in) :: lon
      integer(8), intent(in) :: nx
  
      integer(8) :: i
      real(8) :: dlonr
  
      dlonr = nx/(pi2)
      i = floor(lon*dlonr+1.0d0) ! lon = 2pi/nx*(i-1)=dlon*(i-1)
  
    end function lon2i
  
    function lat2j(lat,ny) result(j)
    ! returns the closest Gaussian point j using approximation
    ! lat varies from NP to SP
    ! 1 <= return value <= ny
      use math_module, only: pi=>math_pi
      implicit none
  
      real(8), intent(in) :: lat
      integer(8) :: ny
  
      integer(8) :: j
  
      j = int(anint(0.5d0*(ny+1-(2.0d0*ny+1.0d0)*lat/pi))) ! lat = (J+1-2j)pi/(2J+1)
  
    end function lat2j
  
    function orthodrome(x1, y1, x2, y2) result(l)
    ! returns the shortest distance between two points on the unit sphere
      implicit none
  
      real(8), intent(in) :: x1, y1, x2, y2
      real(8) :: l, num

      num = cos(x1-x2)*cos(y1)*cos(y2)+sin(y1)*sin(y2)

      if (num > 1.0d0) num = 0.999999999999d0
  
      l = acos(num)
  
    end function orthodrome
  
  end module sphere_module
