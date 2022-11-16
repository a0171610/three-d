module nisl_module

  use grid_module, only: nlon, nlat, ntrunc, nz, lon, pres, coslatr, height, rho, &
    gu, gv, gomega, gphi, gphi_initial, sphi_old, sphi, longitudes=>lon, latitudes=>lat, wgt
  private
  
  real(8), dimension(:, :, :), allocatable, private :: &
    gphi_old, dgphi, dgphim, gphim, gphix, gphiy, gphixy, &
    midlon, midlat, deplon, deplat, gum, gvm, deppres, depheight, &
    gphiz, gphixz, gphiyz, gphixyz
  complex(8), dimension(:, :, :), allocatable, private :: sphi1
  real(8), dimension(:, :, :), allocatable, private :: gw

  private :: update
  public :: nisl_init, nisl_timeint, nisl_clean

contains

  subroutine nisl_init()
    use time_module, only: deltat
    use interpolate3d_module, only: interpolate_init
    use interpolate1d_module, only: interpolate1d_init
    use legendre_transform_module, only: legendre_synthesis
    implicit none

    integer(8) :: i, j, k

    allocate(sphi1(0:ntrunc, 0:ntrunc, nz),gphi_old(nlon, nlat, nz))
    allocate(gphim(nlon, nlat, nz),dgphi(nlon, nlat, nz),dgphim(nlon, nlat, nz))
    allocate(midlon(nlon, nlat, nz), midlat(nlon, nlat, nz), deppres(nlon, nlat, nz))
    allocate(deplon(nlon, nlat, nz), deplat(nlon, nlat, nz))
    allocate(gum(nlon, nlat, nz), gvm(nlon, nlat, nz), depheight(nlon, nlat, nz))
    allocate(gphix(nlon, nlat, nz), gphiy(nlon, nlat, nz), gphixy(nlon, nlat, nz), gw(nlon, nlat, nz))
    allocate(gphiz(nlon, nlat, nz), gphixz(nlon, nlat, nz), gphiyz(nlon, nlat, nz), gphixyz(nlon, nlat, nz))

    call interpolate_init(gphi)
    call interpolate1d_init(gphi(1,1,:))

    do k = 1, nz
      call legendre_synthesis(sphi_old(:,:,k), gphi_old(:,:,k))
    enddo
    gphi(:, :, :) = gphi_old(:, :, :)


    open(11, file="animation.txt")
    do i = 1, nlon
      do j = 1, nlat
          write(11,*) longitudes(i), latitudes(j), gphi(i, j, 25)
      end do        
    end do
    call update(0.0d0, deltat)

  end subroutine nisl_init

  subroutine nisl_clean()
    use interpolate3d_module, only: interpolate_clean
    implicit none

    deallocate(sphi1,gphi_old,gphim,dgphi,dgphim,gum,gvm, &
      midlon,midlat,deplon,deplat)
    call interpolate_clean()

  end subroutine nisl_clean

  subroutine nisl_timeint()
    use time_module, only: nstep, deltat, hstep
    use legendre_transform_module, only: legendre_synthesis
    implicit none

    integer(8) :: i, j, k

    do i = 2, nstep
      call update((i-1)*deltat, 2.0d0*deltat)
      write(*, *) 'step = ', i, "maxval = ", maxval(gphi), 'minval = ', minval(gphi)
      if ( mod(i, hstep) == 0 ) then
        do j = 1, nlon
            do k = 1, nlat
              write(11,*) longitudes(j), latitudes(k), gphi(j, k, 25)
            end do
        end do
      endif
    end do
    close(11)
        open(10, file="log.txt")
    do i = 1, nlon
      do j = 1, nlat
        write(10,*) lon(i), latitudes(j), gphi(i, j, 25)
      enddo
    enddo
    close(10)
    
  end subroutine nisl_timeint

  subroutine update(t, dt)
    use uv_module, only: uv_div
    use time_module, only: case
    use uv_hadley_module, only: uv_hadley
    use planet_module, only: transorm_pressure_to_height, transorm_height_to_pressure
    use upstream3d_module, only: find_points
    use legendre_transform_module, only: legendre_analysis, legendre_synthesis, &
        legendre_synthesis_dlon, legendre_synthesis_dlat, legendre_synthesis_dlonlat
    use interpolate3d_module, only: interpolate_set, interpolate_setd, interpolate_tricubic
    use interpolate1d_module, only: interpolate1d_set, interpolate1d_linear

    implicit none

    integer(8) :: i, j, k, m
    real(8), intent(in) :: t, dt
    real(8) :: tmppres, ans
    real(8), allocatable :: zdot(:, :, :), midh(:, :, :)
    integer(8), allocatable :: id(:, :, :)
    real(8), parameter :: g = 9.80616d0

    allocate(zdot(nlon, nlat, nz), id(nlon, nlat, nz), midh(nlon, nlat, nz))

    select case(case)
      case('hadley')
        call uv_hadley(t, lon, latitudes, pres, gu, gv, gomega)
      case('div')
        call uv_div(t, lon, latitudes, pres, gu, gv, gomega)
      case default
        print *, "No matching initial field"
      stop
    end select
    call find_points(gu, gv, gomega, t, 0.5d0*dt, midlon, midlat, deplon, deplat, deppres)

    do k = 1, nz
      gw(:, :, k) = -gomega(:, :, k) / (g * rho(k))
    enddo

    do k = 1, nz
      do i = 1, nlon
        do j = 1, nlat
          depheight(i, j, k) = transorm_pressure_to_height(deppres(i, j, k))
          call check_height(depheight(i, j, k))
          id(i, j, k) = int(anint(depheight(i, j, k) / 200.0d0)) + 1
          zdot(i, j, k) = gw(i, j, k) - (height(k) - height(id(i, j, k))) / (2.0d0 * dt)
          midh(i, j, k) = (height(k) + height(id(i, j, k))) / 2.0d0
          call check_height(midh(i, j, k))
        enddo
      enddo
    enddo

    do k = 1, nz
      call legendre_synthesis(sphi_old(:, :, k), gphi_old(:, :, k))
    enddo

    ! calculate spectral derivatives

    do k = 1, nz
      call legendre_synthesis_dlon(sphi_old(:, :, k), gphix(:, :, k))
      call legendre_synthesis_dlat(sphi_old(:, :, k), gphiy(:, :, k))
      call legendre_synthesis_dlonlat(sphi_old(:, :, k), gphixy(:, :, k))
    enddo

    do j = 1, nlat
      gphiy(: ,j, :) = gphiy(:, j, :) * coslatr(j)
      gphixy(:, j, :) = gphixy(:, j, :) * coslatr(j)
    end do

    gphiz = 0.0d0; gphixz = 0.0d0; gphiyz = 0.0d0; gphixyz = 0.0d0

    ! set grids
    call interpolate_set(gphi_old)
    call interpolate_setd(gphix, gphiy, gphiz, gphixy, gphixz, gphiyz, gphixyz)
    do j = 1, nlat
      do i = 1, nlon
        do k = 1, nz
          tmppres = transorm_height_to_pressure(height(id(i, j, k)))
          call check_pressure(tmppres)
          call interpolate_tricubic(deplon(i,j,k), deplat(i,j,k), tmppres, gphi(i,j,k))
        enddo
      enddo
    end do

    do k = 2, nz-1
      gphiz(:, :, k) = (gphi(:, :, k + 1) - gphi(:, :, k - 1)) / (height(k+1) - height(k-1))
    end do
    gphiz(:, :, 1) = (gphi(:, :, 2) - gphi(:, :, 1)) / (height(2) - height(1))
    gphiz(:, :, nz) = (gphi(:, :, nz) - gphi(:, :, nz - 1)) / (height(nz) - height(nz-1))

    do i = 1, nlon
      do j = 1, nlat
        call interpolate1d_set(gphiz(i, j, :))
        do k = 1, nz
          call check_height(midh(i, j, k))
          call interpolate1d_linear(midh(i, j, k), ans)
          gphi(i, j, k) = gphi(i, j, k) + zdot(i, j, k) * ans * dt
        enddo
      enddo
    enddo

    do k = 1, nz
      call legendre_analysis(gphi(:,:,k), sphi1(:,:,k))
      do m = 0, ntrunc
        sphi_old(m : ntrunc, m, k) = sphi(m : ntrunc, m, k)
        sphi(m : ntrunc, m, k) = sphi1(m : ntrunc, m, k)
      enddo
    enddo

  end subroutine update

  subroutine check_height(h)
    implicit none
    real(8), intent(inout) :: h
    real(8), parameter :: eps = 1.0d-7

    if (h > 12000.0d0 - eps) then
      h = 12000.d0 - eps
    endif
  end subroutine check_height

  subroutine check_pressure(p)
    implicit none

    real(8), intent(inout) :: p
    real(8), parameter :: pt = 254.944d0, eps = 1.0d-6, ps = 1000.0d0

    if (p < pt - eps) then
      p = pt - eps
    endif
    if (p > ps - eps) then
      p = ps - eps
    endif
  end subroutine check_pressure

end module nisl_module
