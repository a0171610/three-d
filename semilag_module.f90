module semilag_module

  use grid_module, only: nlon, nlat, nz, ntrunc, gphi, height, &
    gu, gv, gw, gphi, sphi_old, sphi, latitudes=>lat, lon, coslatr, wgt
  private
  
  real(8), dimension(:, :, :), allocatable, private :: midlon, midlat, deplon, deplat, deph
  real(8), dimension(:, :, :), allocatable, private :: gphi_old, gphi_initial, gphix, gphiy, gphixy
  real(8), dimension(:, :, :), allocatable, private :: gphiz, gphixz, gphiyz, gphixyz

  private :: update
  public :: semilag_init, semilag_timeint, semilag_clean

contains

  subroutine semilag_init()
    use interpolate3d_module, only: interpolate_init
    use legendre_transform_module, only: legendre_synthesis
    implicit none

    integer(8) :: i,j

    allocate(gphi_old(nlon, nlat, nz), gphix(nlon, nlat, nz),gphiy(nlon, nlat, nz),gphixy(nlon, nlat, nz), &
             midlon(nlon,nlat,nz),midlat(nlon,nlat,nz),deplon(nlon,nlat,nz),deplat(nlon,nlat,nz), gphi_initial(nlon, nlat, nz), &
             deph(nlon, nlat, nz))
    allocate(gphiz(nlon, nlat, nz), gphixz(nlon, nlat, nz), gphiyz(nlon, nlat, nz), gphixyz(nlon, nlat, nz))
    call interpolate_init(gphi)

    gphi_old(:, :, :) = gphi(:, :, :)
    gphi_initial(:, :, :) = gphi(:, :, :)

    do i = 1, nlon
      midlon(i, :, :) = lon(i)
    end do
    do j = 1, nlat
      midlat(:, j, :) = latitudes(j)
    end do
    open(11, file="animation.txt")
    do i = 1, nlat
      do j = 1, nz
        write(11,*) latitudes(i), height(j), gphi(nlon/2, i, j)
      end do        
    end do

  end subroutine semilag_init

  subroutine semilag_clean()
    use interpolate3d_module, only: interpolate_clean
    implicit none

    deallocate(gphi_old,gphix,gphiy,gphixy,midlon,midlat,deplon,deplat)
    call interpolate_clean()

  end subroutine semilag_clean

  subroutine semilag_timeint()
    use math_module, only: pi=>math_pi
    use time_module, only: nstep, deltat, hstep
    use legendre_transform_module, only: legendre_synthesis
    implicit none

    integer(8) :: i, j, k

    open(10, file="uv.txt")
    do i = 1, nlon
      do j = 1, nlat
        write(10,*) lon(i) * 180.0d0 / pi, latitudes(j) * 180.0d0 / pi, gu(i, j, 25)
      end do
    end do
    close(10)

    do i = 1, nstep
      call update((i-0.5d0)*deltat, deltat)
      write(*, *) "step=", i, "maxval = ", maxval(gphi), 'minval = ', minval(gphi)
      if ( mod(i, hstep) == 0 ) then
        do j = 1, nlat
          do k = 1, nz
            write(11,*) latitudes(j), height(k), gphi(nlon/2, j, k)
          end do
        end do
      endif
    end do
    close(11)
    open(10, file="log.txt")
    do i = 1, nlat
      do j = 1, nz
        write(10,*) latitudes(i), height(j), gphi(nlon/2, i, j)
      enddo
    enddo
    close(10)
  end subroutine semilag_timeint

  subroutine update(t, dt)
    use math_module, only: &
      pi=>math_pi, pir=>math_pir, pih=>math_pih
    use upstream3d_module, only: find_points
    use time_module, only: case, imethod
    use uv_module, only: uv_div
    use uv_hadley_module, only: uv_hadley
    use interpolate3d_module, only: interpolate_set, interpolate_setd, interpolate_tricubic
    use legendre_transform_module, only: legendre_analysis, legendre_synthesis, &
        legendre_synthesis_dlon, legendre_synthesis_dlat, legendre_synthesis_dlonlat
    implicit none

    real(8), intent(in) :: t, dt

    integer(8) :: i, j, k

    select case(case)
      case('hadley')
        call uv_hadley(t, lon, latitudes, height, gu, gv, gw)
      case('div')
        call uv_div(t, lon, latitudes, height, gu, gv, gw)
      case default
        print *, "No matching initial field"
      stop
    end select

    call find_points(gu, gv, gw, t, 0.5d0*dt, midlon, midlat, deplon, deplat, deph)
    do k = 1, nz
      call legendre_analysis(gphi_old(:, :, k), sphi_old(:, :, k))
    enddo

! calculate spectral derivatives

    select case(imethod)
    case('sph')
      call sph_derivative
    case('fd')
      call fd_derivative
    case default
      write(*,*) "no mathing imethod for ", imethod
    end select

    call vertical_derivative

! set grids
    call interpolate_set(gphi_old)
    call interpolate_setd(gphix, gphiy, gphiz, gphixy, gphixz, gphiyz, gphixyz)
    do j = 1, nlat
      do i = 1, nlon
        do k = 1, nz
          call interpolate_tricubic(deplon(i,j,k), deplat(i,j,k), deph(i, j, k), gphi(i,j,k))
        enddo
      enddo
    end do

! spectral
    gphi_old = gphi

  end subroutine update

  subroutine sph_derivative
    use legendre_transform_module, only: legendre_analysis, legendre_synthesis, &
        legendre_synthesis_dlon, legendre_synthesis_dlat, legendre_synthesis_dlonlat
    implicit none
    integer(8) :: j, k

    do k = 1, nz
      call legendre_synthesis_dlon(sphi_old(:, :, k), gphix(:, :, k))
      call legendre_synthesis_dlat(sphi_old(:, :, k), gphiy(:, :, k))
      call legendre_synthesis_dlonlat(sphi_old(:, :, k), gphixy(:, :, k))
    enddo

    do j = 1, nlat
      gphiy(: ,j, :) = gphiy(:, j, :) * coslatr(j)
      gphixy(:, j, :) = gphixy(:, j, :) * coslatr(j)
    end do
  end subroutine sph_derivative

  subroutine fd_derivative
    use math_module, only: pir=>math_pir, pih=>math_pih
    implicit none
    integer(8) :: i, j, k
    real(8) :: dlonr, eps, gphitmp(nlon)

    do k = 1, nz
      dlonr = 0.25d0*nlon*pir
      gphix(1,:,k) = dlonr * (gphi_old(2,:,k) - gphi_old(nlon,:,k))
      gphix(nlon, :, k) = dlonr * (gphi_old(1,:,k) - gphi_old(nlon-1,:,k))
      do i=2, nlon-1
        gphix(i,:,k) = dlonr*(gphi_old(i+1,:,k) - gphi_old(i-1,:,k))
      end do
    end do

    do k = 1, nz
      eps = pih-latitudes(1)
      gphitmp = cshift(gphi_old(:,1,k),nlon/2)
      gphiy(:,1,k) = (gphitmp-gphi_old(:,2,k))/(pih+eps-latitudes(2))
      gphitmp = cshift(gphix(:,1,k),nlon/2)
      gphixy(:,1,k) = (gphitmp-gphix(:,2,k))/(pih+eps-latitudes(2))
      gphitmp = cshift(gphi_old(:,nlat,k),nlon/2)
      gphiy(:,nlat,k) = (gphitmp-gphi_old(:,nlat-1,k))/(-pih-eps-latitudes(nlat-1))
      gphitmp = cshift(gphix(:,nlat,k),nlon/2)
      gphixy(:,nlat,k) = (gphitmp-gphix(:,nlat-1,k))/(-pih-eps-latitudes(nlat-1))
      do j=2, nlat-1
        gphiy(:,j,k) = (gphi_old(:,j+1,k)-gphi_old(:,j-1,k))/(latitudes(j+1)-latitudes(j-1))
        gphixy(:,j,k) = (gphix(:,j+1,k)-gphix(:,j-1,k))/(latitudes(j+1)-latitudes(j-1))
      end do
    end do
  end subroutine fd_derivative

  subroutine vertical_derivative
    implicit none

    integer(8) :: k

    do k = 2, nz-1
      gphiz(:, :, k) = (gphi_old(:, :, k + 1) - gphi_old(:, :, k - 1)) / (height(k+1) - height(k-1))
      gphixz(:, :, k) = (gphix(:, :, k + 1) - gphix(:, :, k - 1)) / (height(k+1) - height(k-1))
      gphiyz(:, :, k) = (gphiy(:, :, k + 1) - gphiy(:, :, k - 1)) / (height(k+1) - height(k-1))
      gphixyz(:, :, k) = (gphixy(:, :, k + 1) - gphixy(:, :, k - 1)) / (height(k+1) - height(k-1))
    end do

    gphiz(:, :, 1) = (gphi_old(:, :, 2) - gphi_old(:, :, 1)) / (height(2) - height(1))
    gphixz(:, :, 1) = (gphix(:, :, 2) - gphix(:, :, 1)) / (height(2) - height(1))
    gphiyz(:, :, 1) = (gphiy(:, :, 2) - gphiy(:, :, 1)) / (height(2) - height(1))
    gphixyz(:, :, 1) = (gphixy(:, :, 2) - gphixy(:, :, 1)) / (height(2) - height(1))

    gphiz(:, :, nz) = (gphi_old(:, :, nz) - gphi_old(:, :, nz - 1)) / (height(nz) - height(nz-1))
    gphixz(:, :, nz) = (gphix(:, :, nz) - gphix(:, :, nz - 1)) / (height(nz) - height(nz-1))
    gphiyz(:, :, nz) = (gphiy(:, :, nz) - gphiy(:, :, nz - 1)) / (height(nz) - height(nz-1))
    gphixyz(:, :, nz) = (gphixy(:, :, nz) - gphixy(:, :, nz - 1)) / (height(nz) - height(nz-1))
  end subroutine vertical_derivative

end module semilag_module