module nisl_module

  use grid_module, only: nlon, nlat, ntrunc, nz, lon, coslatr, height, dh, &
    gu, gv, gw, gphi, gphi_initial, longitudes=>lon, latitudes=>lat, wgt
  private
  
  real(8), dimension(:, :, :), allocatable, private :: &
    gphi_old, dgphi, dgphim, gphim, gphix, gphiy, gphixy, &
    midlon, midlat, deplon, deplat, gum, gvm, depheight, &
    gphiz, gphi1

  private :: update
  public :: nisl_init, nisl_timeint, nisl_clean

contains

  subroutine nisl_init()
    use time_module, only: deltat
    use interpolate2d_module, only: interpolate2d_init
    use interpolate1d_module, only: interpolate1d_init
    use legendre_transform_module, only: legendre_synthesis
    implicit none

    integer(8) :: i, j

    allocate(gphi_old(nlon, nlat, nz), gphi1(nlon, nlat, nz))
    allocate(gphim(nlon, nlat, nz),dgphi(nlon, nlat, nz),dgphim(nlon, nlat, nz))
    allocate(midlon(nlon, nlat, nz), midlat(nlon, nlat, nz))
    allocate(deplon(nlon, nlat, nz), deplat(nlon, nlat, nz))
    allocate(gum(nlon, nlat, nz), gvm(nlon, nlat, nz), depheight(nlon, nlat, nz))
    allocate(gphix(nlon, nlat, nz), gphiy(nlon, nlat, nz), gphixy(nlon, nlat, nz))
    allocate(gphiz(nlon, nlat, nz))

    call interpolate2d_init(gphi)
    call interpolate1d_init(gphi(1,1,:))

    gphi_old = gphi
    gphi1 = gphi

    open(11, file="animation.txt")
    do i = 1, nlat
      do j = 1, nz
          write(11,*) latitudes(i), height(j), gphi(nlon/2, i, j)
      end do        
    end do
    call update(0.0d0, deltat)

  end subroutine nisl_init

  subroutine nisl_clean()
    use interpolate2d_module, only: interpolate2d_clean
    implicit none

    deallocate(gphi_old,gphim,dgphi,dgphim,gum,gvm, &
      midlon,midlat,deplon,deplat)
    call interpolate2d_clean()

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
    
  end subroutine nisl_timeint

  subroutine update(t, dt)
    use uv_module, only: uv_div
    use time_module, only: case
    use uv_hadley_module, only: uv_hadley, calc_w
    use upstream3d_module, only: find_points
    use legendre_transform_module, only: legendre_analysis, legendre_synthesis, &
        legendre_synthesis_dlon, legendre_synthesis_dlat, legendre_synthesis_dlonlat
    use interpolate2d_module, only: interpolate2d_set, interpolate2d_setd, interpolate2d_bicubic
    use interpolate1d_module, only: interpolate1d_set, interpolate1d_linear

    implicit none

    integer(8) :: i, j, k
    real(8), intent(in) :: t, dt
    real(8) :: ans
    real(8), allocatable :: zdot(:, :, :), midh(:, :, :)
    integer(8), allocatable :: id(:, :, :)

    allocate(zdot(nlon, nlat, nz), id(nlon, nlat, nz), midh(nlon, nlat, nz))

    select case(case)
      case('hadley')
        call uv_hadley(t, lon, latitudes, height, gu, gv, gw)
      case('div')
        call uv_div(t, lon, latitudes, height, gu, gv, gw)
      case default
        print *, "No matching initial field"
      stop
    end select
    call find_points(gu, gv, gw, t, 0.5d0*dt, midlon, midlat, deplon, deplat, depheight)

    do k = 1, nz
      do i = 1, nlon
        do j = 1, nlat
          call check_height(depheight(i, j, k))
          id(i, j, k) = int(anint(depheight(i, j, k) / dh)) + 1
          zdot(i, j, k) = gw(i, j, k) - (height(k) - height(id(i, j, k))) / (2.0d0 * dt)
          midh(i, j, k) = (height(k) + height(id(i, j, k))) / 2.0d0
          call check_height(midh(i, j, k))
        enddo
      enddo
    enddo

    call fd_derivative

    ! set grids
    call interpolate2d_set(gphi_old)
    call interpolate2d_setd(gphix, gphiy, gphixy)
    do j = 1, nlat
      do i = 1, nlon
        do k = 1, nz
          call check_height(height(id(i, j, k)))
          call interpolate2d_bicubic(deplon(i,j,k), deplat(i,j,k), height(id(i, j, k)), gphi(i,j,k))
        enddo
      enddo
    end do

    do k = 2, nz-1
      gphiz(:, :, k) = (gphi1(:, :, k + 1) - gphi1(:, :, k - 1)) / (height(k+1) - height(k-1))
    end do
    gphiz(:, :, 1) = (gphi1(:, :, 2) - gphi1(:, :, 1)) / (height(2) - height(1))
    gphiz(:, :, nz) = (gphi1(:, :, nz) - gphi1(:, :, nz - 1)) / (height(nz) - height(nz-1))

    do i = 1, nlon
      do j = 1, nlat
        call interpolate1d_set(gphiz(i, j, :))
        do k = 1, nz
          call check_height(midh(i, j, k))
          call interpolate1d_linear(midh(i, j, k), ans)
          gphi(i, j, k) = gphi(i, j, k) + (height(k) - height(id(i, j, k))) * ans
        enddo
      enddo
    enddo

    do i = 1, nlon
      do j = 1, nlat
        do k = 1, nz
          gphiz(i, j, k) = gphiz(i, j, k) * calc_w(latitudes(j), height(k), t)
        end do
      end do
    end do

    do i = 1, nlon
      do j = 1, nlat
        call interpolate1d_set(gphiz(i, j, :))
        do k = 1, nz
          call check_height(midh(i, j, k))
          call interpolate1d_linear(midh(i, j, k), ans)
          gphi(i, j, k) = gphi(i, j, k) - ans * dt
        enddo
      enddo
    enddo

    gphi_old(:,:,:) = gphi1(:,:,:)
    gphi1(:,:,:) = gphi(:,:,:)

  end subroutine update

  subroutine check_height(h)
    implicit none
    real(8), intent(inout) :: h
    real(8), parameter :: eps = 1.0d-7

    if (h > 12000.0d0 - eps) then
      h = 12000.d0 - eps
    endif
  end subroutine check_height

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

end module nisl_module
