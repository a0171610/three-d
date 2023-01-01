module nisl_module

  use grid_module, only: nlon, nlat, ntrunc, nz, lon, coslatr, height, dh, &
    gu, gv, gw, gphi, gphi_initial, longitudes=>lon, latitudes=>lat, wgt, sphi_old
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

    implicit none

    integer(8) :: i, j, k
    real(8), intent(in) :: t, dt
    real(8) :: ans1, ans2
    integer(8), allocatable :: id(:, :, :)

    allocate(id(nlon, nlat, nz))

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
        enddo
      enddo
    enddo

    call fd_derivative(gphi_old)
    call mid_points

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

    call fd_derivative(gphiz)
    call interpolate2d_set(gphiz)
    call interpolate2d_setd(gphix, gphiy, gphixy)
    do k = 1, nz
      do i = 1, nlon
        do j = 1, nlat
          call interpolate2d_bicubic(midlon(i,j,k), midlat(i,j,k), height(k), ans1)
          call interpolate2d_bicubic(midlon(i,j,k), midlat(i,j,k), height(id(i,j,k)), ans2)
          gphi(i, j, k) = gphi(i, j, k) + (height(k) - height(id(i, j, k))) * (ans1 + ans2) * 0.5d0
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

    call fd_derivative(gphiz)
    call interpolate2d_set(gphiz)
    call interpolate2d_setd(gphix, gphiy, gphixy)
    do k = 1, nz
      do i = 1, nlon
        do j = 1, nlat
          call interpolate2d_bicubic(midlon(i,j,k), midlat(i,j,k), height(k), ans1)
          call interpolate2d_bicubic(midlon(i,j,k), midlat(i,j,k), height(id(i,j,k)), ans2)
          gphi(i, j, k) = gphi(i, j, k) - (ans1 + ans2) * 0.5d0 * dt
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

  subroutine mid_points
    use sphere_module, only: lonlat2xyz
    use math_module, only: pi2=>math_pi2
    implicit none
    integer(8) :: i, j, k
    real(8) :: xg, yg, zg, xr, yr, zr, xm, ym, zm

    do k = 1, nz
      do i = 1, nlon
        do j = 1, nlat
          call lonlat2xyz(lon(i), latitudes(j), xg, yg, zg)
          call lonlat2xyz(deplon(i, j, k), deplat(i, j, k), xr, yr, zr)
          xm = (xg + xr) * 0.5d0
          ym = (yg + yr) * 0.5d0
          zm = (zg + zr) * 0.5d0
          midlon(i, j, k) = modulo(atan2(ym, xm) + pi2, pi2)
          midlat(i, j, k) = asin(zm / sqrt(xm**2 + ym**2 + zm**2))
        enddo
      enddo
    enddo
  end subroutine mid_points

  subroutine fd_derivative(f)
    use math_module, only: pir=>math_pir, pih=>math_pih
    implicit none
    real(8), intent(in) :: f(nlon, nlat, nz)
    integer(8) :: i, j, k
    real(8) :: dlonr, eps, gphitmp(nlon)

    do k = 1, nz
      dlonr = 0.25d0*nlon*pir
      gphix(1,:,k) = dlonr * (f(2,:,k) - f(nlon,:,k))
      gphix(nlon, :, k) = dlonr * (f(1,:,k) - f(nlon-1,:,k))
      do i=2, nlon-1
        gphix(i,:,k) = dlonr*(f(i+1,:,k) - f(i-1,:,k))
      end do
    end do

    do k = 1, nz
      eps = pih-latitudes(1)
      gphitmp = cshift(f(:,1,k),nlon/2)
      gphiy(:,1,k) = (gphitmp-f(:,2,k))/(pih+eps-latitudes(2))
      gphitmp = cshift(gphix(:,1,k),nlon/2)
      gphixy(:,1,k) = (gphitmp-gphix(:,2,k))/(pih+eps-latitudes(2))
      gphitmp = cshift(f(:,nlat,k),nlon/2)
      gphiy(:,nlat,k) = (gphitmp-f(:,nlat-1,k))/(-pih-eps-latitudes(nlat-1))
      gphitmp = cshift(gphix(:,nlat,k),nlon/2)
      gphixy(:,nlat,k) = (gphitmp-gphix(:,nlat-1,k))/(-pih-eps-latitudes(nlat-1))
      do j=2, nlat-1
        gphiy(:,j,k) = (f(:,j+1,k)-f(:,j-1,k))/(latitudes(j+1)-latitudes(j-1))
        gphixy(:,j,k) = (gphix(:,j+1,k)-gphix(:,j-1,k))/(latitudes(j+1)-latitudes(j-1))
      end do
    end do
  end subroutine fd_derivative

  subroutine sph_derivative(f)
    use legendre_transform_module, only: legendre_analysis, legendre_synthesis, &
        legendre_synthesis_dlon, legendre_synthesis_dlat, legendre_synthesis_dlonlat
    implicit none
    integer(8) :: j, k
    real(8) :: f(nlon, nlat, nz)

    do k = 1, nz
      call legendre_analysis(f(:,:,k), sphi_old(:,:,k))
      call legendre_synthesis_dlon(sphi_old(:, :, k), gphix(:, :, k))
      call legendre_synthesis_dlat(sphi_old(:, :, k), gphiy(:, :, k))
      call legendre_synthesis_dlonlat(sphi_old(:, :, k), gphixy(:, :, k))
    enddo

    do j = 1, nlat
      gphiy(: ,j, :) = gphiy(:, j, :) * coslatr(j)
      gphixy(:, j, :) = gphixy(:, j, :) * coslatr(j)
    end do
  end subroutine sph_derivative

end module nisl_module
