module grid_module
  implicit none
  private

  integer(8), parameter, public :: nz = 61
  integer(8), parameter, public ::  ntrunc = 39, nlon = 120, nlat = 60
  !integer(8), parameter, public ::  ntrunc = 79, nlon = 240, nlat = 120
  !integer(8), parameter, public ::  ntrunc = 159, nlon = 480, nlat = 240
  !integer(8), parameter, public ::  ntrunc = 319, nlon = 960, nlat = 480

  complex(8), dimension(:, :, :), allocatable, public :: sphi, sphi_old
  real(8), dimension(:, :, :), allocatable, public :: gphi, gphi_initial, gu, gv, gomega
  real(8), dimension(:), allocatable, public :: lon, lat, coslat, coslatr, wgt, pres, sigma, height, rho
  real(8), public :: Umax, dlat(nlat), dlat4(nlat)
  real(8), public, parameter :: ps = 1000.0d0, Rd = 287.0d-3, T0 = 300.0d0

  public :: grid_init, grid_clean, pole_regrid

contains

  subroutine grid_init()
    use math_module, only: pi2=>math_pi2, pih=>math_pih, math_pi
    use legendre_transform_module, only: &
      legendre_init, legendre_analysis
    use init_module, only: &
      init_ghill, init_ghill2, init_cbell2, init_scyli2, init_ccbel2, init_cbell_2, init_hadley
    use uv_module, only: uv_div
    use uv_hadley_module, only: uv_hadley
    use time_module, only: case
    use planet_module, only: planet_radius, transorm_height_to_pressure
    implicit none


    integer(8) :: i, j, k, m, n
    real(8) :: dlon, eps

    allocate(lon(nlon), lat(nlat), coslat(nlat), coslatr(nlat), wgt(nlat), pres(nz), sigma(nz), height(nz), rho(nz))
    allocate(gphi(nlon,nlat, nz), gphi_initial(nlon, nlat, nz), gu(nlon,nlat, nz), gv(nlon,nlat, nz))
    allocate(sphi(0:ntrunc,0:ntrunc, nz), sphi_old(0:ntrunc,0:ntrunc, nz), gomega(nlon, nlat, nz))
 
    dlon = pi2/nlon
    do i=1, nlon
      lon(i) = dlon * dble(i-1)
    end do
    call legendre_init(nlon,nlat,ntrunc,lat,wgt)
    coslat(:) = cos(lat(:))
    coslatr(:) = 1.0d0 / coslat(:)
    eps = pih-lat(1)
    dlat(1) = pih+eps-lat(2)
    dlat(nlat) = -pih-eps-lat(nlat-1)
    do j = 2, nlat - 1
      dlat(j) = lat(j+1) - lat(j-1)
    end do

    dlat4(1) = -(math_pi - lat(3) + lat(nlat - 1))
    dlat4(2) = -(math_pi - lat(4) + lat(nlat))
    dlat4(nlat-1) = -(math_pi - lat(1) + lat(nlat - 3))
    dlat4(nlat) = -(math_pi - lat(2) + lat(nlat - 2))
    do j = 3, nlat-2
      dlat4(j) = lat(j+2) - lat(j - 2)
    end do

    do i = 1, nz
      height(i) = dble(i - 1) * 200.0d0
      pres(i) = transorm_height_to_pressure(height(i))
      sigma(i) = pres(i) / ps
      rho(i) = pres(i) / (Rd * T0)
    end do

    select case(case)
    case('hadley')
      call init_hadley(lon, lat, height, gphi)
    case('div')
      call init_cbell_2(lon, lat, height, gphi)
    case default
      print *, "No matching initial field"
    stop
  end select

    gphi_initial(:, :, :) = gphi(:, :, :)
    do k = 1, nz
      call legendre_analysis(gphi(:,:,k), sphi(:,:,k))
    enddo
      do m = 0, ntrunc
      do n = m, ntrunc
        do k = 1, nz
          sphi_old(n, m, k) = sphi(n, m, k)
        enddo
      enddo
    enddo

    select case(case)
      case('hadley')
        call uv_hadley(0.0d0, lon, lat, pres, gu, gv, gomega)
      case('div')
        call uv_div(0.0d0, lon, lat, pres, gu, gv, gomega)
      case default
        print *, "No matching initial field"
      stop
    end select

  Umax = real(maxval(gu) * planet_radius)

  end subroutine grid_init

  subroutine grid_clean
    implicit none

    deallocate(lon, lat, coslat, coslatr, gphi, gu, gv, sphi, sphi_old)

  end subroutine grid_clean

  subroutine pole_regrid(x, y)
    implicit none
    integer(8), intent(inout) :: x, y
    if (x < 1) then
        x = x + nlon
    endif
    if (x > nlon) then
        x = x - nlon
    endif
    if ( y < 1 ) then
        y = 1 - y
        if ( x > nlon/2) then
            x = x - nlon/2
        else
            x = x + nlon/2
        endif
    endif
    if ( y > nlat ) then
        y = 2 * nlat + 1 - y
        if ( x > nlon/2) then
            x = x - nlon/2
        else
            x = x + nlon/2
        endif
    endif
end subroutine pole_regrid

end module grid_module
