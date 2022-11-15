module init_module
  use math_module, only: pi=>math_pi
  implicit none
  private

  integer, parameter, private :: nc = 2
  real(8), dimension(nc), private :: &
    lonc = (/5.0d0*pi/6.0d0, 7.0d0*pi/6.0d0/), &
    latc = (/0.0d0, 0.0d0/)

  public :: init_ghill, init_ghill2, init_cbell2, init_scyli2, init_ccbel2, init_cbell_2

contains

! Ritche 1987
  subroutine init_ghill(lon,lat,gphi)
    use planet_module, only: a=>planet_radius
    use math_module, only: deg2rad=>math_deg2rad
    use sphere_module, only: orthodrome
    implicit none

    real(8), dimension(:), intent(in) :: lon, lat
    real(8), dimension(:,:), intent(inout) :: gphi

    real(8), parameter :: &
      phimax = 100.0d0, xi = 0.0d0, yi = 0.0d0, L = 2.5d6 ! Gaussian hill

    integer(8) :: i, j, nx, ny
    real(8) :: r, loni, lati

    nx = size(lon)
    ny = size(lat)
    loni = xi*deg2rad
    lati = yi*deg2rad
    do j=1, ny
      do i=1, nx
        r = a*orthodrome(lon(i),lat(j),loni,lati)
        gphi(i,j) = phimax * exp(-(pi*r/L)**2)
      end do
    end do

  end subroutine init_ghill

  subroutine init_ghill2(lon,lat,gphi)
    use sphere_module, only: lonlat2xyz
    implicit none

    real(8), dimension(:), intent(in) :: lon, lat
    real(8), dimension(:,:), intent(inout) :: gphi

    real(8), parameter :: hmax = 0.95d0, b0 = 5.0d0

    integer(8) :: i, j, k, nx, ny
!    real(8) :: x0, y0, z0, x, y, z

    nx = size(lon)
    ny = size(lat)
    gphi(:,:) = 0.0d0
    do k=1, nc
!      call lonlat2xyz(lonc(k),latc(k),x0,y0,z0)
      do j=1, ny
        do i=1, nx
!          call lonlat2xyz(lon(i),lat(j),x,y,z)
!          gphi(i,j) = gphi(i,j) + hmax*exp(-b0*((x-x0)**2+(y-y0)**2+(z-z0)**2))
          gphi(i,j) = gphi(i,j) + hmax*exp(-2*b0*(1.0d0-cos(lat(j))*cos(lon(i)-lonc(k))))
        end do
      end do
    end do

  end subroutine init_ghill2

  subroutine init_cbell2(lon,lat,gphi)
    use sphere_module, only: orthodrome
    implicit none

    real(8), dimension(:), intent(in) :: lon, lat
    real(8), dimension(:,:), intent(inout) :: gphi

    real(8), parameter :: &
      hmax = 1.0d0, r0 = 0.5d0, b = 0.1d0, c = 0.9d0

    integer(8) :: i, j, k, nx, ny
    real(8) :: hch, pir0r, r

    nx = size(lon)
    ny = size(lat)
    gphi(:,:) = b
    pir0r = pi/r0
    hch = c*0.5d0*hmax
    do k=1, nc
      do j=1, ny
        do i=1, nx
          r = orthodrome(lon(i),lat(j),lonc(k),latc(k))
          if (r<r0) then
            gphi(i,j) = gphi(i,j) + hch*(1.0d0+cos(pir0r*r))
          end if
        end do
      end do
    end do

  end subroutine init_cbell2

  subroutine init_cbell_2(lon, lat, z, gphi)
    use sphere_module, only: orthodrome
    implicit none

    real(8), dimension(:), intent(in) :: lon, lat, z
    real(8), dimension(:, :, :), intent(inout) :: gphi

    real(8), parameter :: &
      hmax = 1.0d0, r0 = 0.5d0, b = 0.1d0, c = 0.9d0

    integer(8) :: i, j, id, k, nx, ny, nz
    real(8) :: r, dist

    nx = size(lon)
    ny = size(lat)
    nz = size(z)

    do id = 1, 2
      do i = 1, nx
        do j = 1, ny
          do k = 1, nz
            r = orthodrome(lon(i), lat(j), lonc(id), latc(id))
            dist = threed_dist(r, z(k))
            gphi(i, j, k) = gphi(i, j, k) + (1.0d0 + cos(dist * pi)) / 2.0d0
          end do
        end do
      end do
    end do
  end subroutine init_cbell_2

  subroutine init_scyli2(lon,lat,gphi)
    use sphere_module, only: orthodrome
    implicit none

    real(8), dimension(:), intent(in) :: lon, lat
    real(8), dimension(:,:), intent(inout) :: gphi

    real(8), parameter :: &
      hmax = 1.0d0, r0 = 0.5d0, b = 0.1d0, c = 1.0d0

    integer(8) :: i, j, nx, ny
    real(8) :: &
      r1, r2, dlon0, dlon1, dlon2, dlat0, dlat1, dlat2

    nx = size(lon)
    ny = size(lat)
    gphi(:,:) = b
    dlon0 = r0/6.0d0
    dlat0 = 5.0d0/12.0d0*r0
    do j=1, ny
      do i=1, nx
        r1 = orthodrome(lon(i),lat(j),lonc(1),latc(1))
        r2 = orthodrome(lon(i),lat(j),lonc(2),latc(2))
        dlon1 = abs(lon(i)-lonc(1))
        dlon2 = abs(lon(i)-lonc(2))
        dlat1 = lat(j)-latc(1)
        dlat2 = lat(j)-latc(2)
        if (((r1<=r0).and.(dlon1>=dlon0)).or. &
            ((r2<=r0).and.(dlon2>=dlon0)).or. &
            ((r1<=r0).and.(dlon1<dlon0).and.dlat1<-dlat0).or. &
            ((r2<=r0).and.(dlon2<dlon0).and.dlat2>dlat0)) then
          gphi(i,j) = c
        end if
      end do
    end do

  end subroutine init_scyli2

  subroutine init_ccbel2(lon,lat,gphi)
    use sphere_module, only: orthodrome
    implicit none

    real(8), dimension(:), intent(in) :: lon, lat
    real(8), dimension(:,:), intent(inout) :: gphi

    real(8), parameter :: a = -0.8, b = 0.9

    call init_cbell2(lon,lat,gphi)
    gphi(:,:) = a*gphi(:,:)**2 + b

  end subroutine init_ccbel2

  function threed_dist(r, z) result(ans)
    implicit none
    real(8), intent(in) :: r, z
    real(8) :: ans
    real(8), parameter :: Zt = 1000.0d0, zc = 5000.0d0, Rt = 0.5d0
    ans = (r / Rt) ** 2 + ((z - zc)/ Zt) ** 2
    ans = min(ans, 1.0d0)
  end function threed_dist

end module init_module
