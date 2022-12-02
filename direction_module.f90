module direction_module

  use grid_module, only: nlon, nlat, ntrunc, nz, lon, coslatr, height, dh, &
    gu, gv, gw, gphi, gphi_initial, sphi_old, sphi, longitudes=>lon, latitudes=>lat, wgt
  use time_module, only: case
  private
  
  real(8), dimension(:, :, :), allocatable, private :: &
    gphi_old, dgphi, dgphim, gphim, gphix, gphiy, gphixy, &
    midlon, midlat, deplon, deplat, gum, gvm, depheight, &
    gphiz, gphixz, gphiyz, gphixyz
  complex(8), dimension(:, :, :), allocatable, private :: sphi1

  private :: update
  public :: direction_init, direction_timeint, direction_clean

contains

  subroutine direction_init()
    use time_module, only: deltat
    use interpolate3d_module, only: interpolate_init
    use interpolate1d_module, only: interpolate1d_init
    use legendre_transform_module, only: legendre_synthesis
    implicit none

    integer(8) :: i, j, k

    allocate(sphi1(0:ntrunc, 0:ntrunc, nz),gphi_old(nlon, nlat, nz))
    allocate(gphim(nlon, nlat, nz),dgphi(nlon, nlat, nz),dgphim(nlon, nlat, nz))
    allocate(midlon(nlon, nlat, nz), midlat(nlon, nlat, nz))
    allocate(deplon(nlon, nlat, nz), deplat(nlon, nlat, nz))
    allocate(gum(nlon, nlat, nz), gvm(nlon, nlat, nz), depheight(nlon, nlat, nz))
    allocate(gphix(nlon, nlat, nz), gphiy(nlon, nlat, nz), gphixy(nlon, nlat, nz))
    allocate(gphiz(nlon, nlat, nz), gphixz(nlon, nlat, nz), gphiyz(nlon, nlat, nz), gphixyz(nlon, nlat, nz))

    call interpolate_init(gphi)
    call interpolate1d_init(gphi(1,1,:))

    do k = 1, nz
      call legendre_synthesis(sphi_old(:,:,k), gphi_old(:,:,k))
    enddo
    gphi(:, :, :) = gphi_old(:, :, :)

    open(11, file="animation.txt")
    do i = 1, nlat
      do j = 1, nz
          write(11,*) latitudes(i), height(j), gphi(nlon/2, i, j)
      end do        
    end do
    call update(0.0d0, deltat)

  end subroutine direction_init

  subroutine direction_clean()
    use interpolate3d_module, only: interpolate_clean
    implicit none

    deallocate(sphi1,gphi_old,gphim,dgphi,dgphim,gum,gvm, &
      midlon,midlat,deplon,deplat)
    call interpolate_clean()

  end subroutine direction_clean

  subroutine direction_timeint()
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
      if ( i == nstep/2 ) then
        open(10, file="log.txt")
        do j = 1, nlat
          do k = 1, nz
            write(10,*) latitudes(j), height(k), gphi(nlon/2, j, k)
          enddo
        enddo
        close(10)
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
    
  end subroutine direction_timeint

  subroutine update(t, dt)
    use uv_module, only: uv_div
    use uv_hadley_module, only: uv_hadley, calc_w
    use upstream3d_module, only: find_points
    use legendre_transform_module, only: legendre_analysis, legendre_synthesis, &
        legendre_synthesis_dlon, legendre_synthesis_dlat, legendre_synthesis_dlonlat
    use interpolate3d_module, only: interpolate_set, interpolate_setd, interpolate_tricubic
    use interpolate1d_module, only: interpolate1d_set, interpolate1d_linear

    implicit none

    integer(8) :: i, j, k, m
    real(8), intent(in) :: t, dt
    real(8) :: ans, tmp1, tmp2
    real(8), allocatable :: zdotA(:, :, :), zdotB(:, :, :), midhA(:, :, :), midhB(:, :, :), ratio(:, :, :), gphi1(:, :, :)
    integer(8), allocatable ::  idA(:, :, :), idB(:, :, :)
    real(8), parameter :: g = 9.80616d0

    allocate(zdotA(nlon, nlat, nz), zdotB(nlon, nlat, nz), midhA(nlon, nlat, nz), midhB(nlon, nlat, nz), gphi1(nlon, nlat, nz))
    allocate(idA(nlon, nlat, nz), idB(nlon, nlat, nz), ratio(nlon, nlat, nz))

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
          idA(i, j, k) = int(aint(depheight(i, j, k) / dh)) + 1
          idB(i, j, k) = idA(i, j, k) + 1
          ratio(i, j, k) = calculate_ratio(depheight(i,j,k) - height(idA(i,j,k)), height(idB(i,j,k)) - depheight(i,j,k))

          zdotA(i, j, k) = gw(i, j, k) - (height(k) - height(idA(i, j, k))) / (2.0d0 * dt)
          zdotB(i, j, k) = gw(i, j, k) - (height(k) - height(idB(i, j, k))) / (2.0d0 * dt)
          midhA(i, j, k) = (height(k) + height(idA(i, j, k))) / 2.0d0
          midhB(i, j, k) = (height(k) + height(idB(i, j, k))) / 2.0d0
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
          call check_height(height(idA(i, j, k)))
          call check_height(height(idB(i, j, k)))
          call interpolate_tricubic(deplon(i,j,k), deplat(i,j,k), height(idA(i, j, k)), tmp1)
          call interpolate_tricubic(deplon(i,j,k), deplat(i,j,k), height(idB(i, j, k)), tmp2)
          gphi(i, j, k) = tmp1 * ratio(i, j, k) + tmp2 * (1.0d0 - ratio(i, j, k))
        enddo
      enddo
    end do

    do k = 1, nz
      call legendre_synthesis(sphi(:, :, k), gphi1(:, :, k))
    enddo

    do k = 2, nz-1
      gphiz(:, :, k) = (gphi1(:, :, k + 1) - gphi1(:, :, k - 1)) / (height(k+1) - height(k-1))
    end do
    gphiz(:, :, 1) = (gphi1(:, :, 2) - gphi1(:, :, 1)) / (height(2) - height(1))
    gphiz(:, :, nz) = (gphi1(:, :, nz) - gphi1(:, :, nz - 1)) / (height(nz) - height(nz-1))

    do i = 1, nlon
      do j = 1, nlat
        call interpolate1d_set(gphiz(i, j, :))
        do k = 1, nz
          call check_height(midhA(i, j, k))
          call interpolate1d_linear(midhA(i, j, k), ans)
          gphi(i, j, k) = gphi(i, j, k) + (height(k) - height(idA(i, j, k))) * ans * ratio(i, j, k)

          call check_height(midhB(i, j, k))
          call interpolate1d_linear(midhB(i, j, k), ans)
          gphi(i, j, k) = gphi(i, j, k) + (height(k) - height(idB(i, j, k))) * ans * (1.0d0 - ratio(i, j, k))
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
          call check_height(midhA(i, j, k))
          call interpolate1d_linear(midhA(i, j, k), ans)
          gphi(i, j, k) = gphi(i, j, k) - ans * dt * ratio(i, j, k)

          call check_height(midhB(i, j, k))
          call interpolate1d_linear(midhB(i, j, k), ans)
          gphi(i, j, k) = gphi(i, j, k) - ans * dt * (1.0d0 - ratio(i, j, k))
        enddo
      enddo
    enddo

    do i = 1, nlon
      do j = 1, nlat
        do k = 1, nz
          if ( gphi(i, j, k) < -0.05d0 ) then
            gphi(i, j, k) = -0.05d0
          endif
        end do
      end do
    end do

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

  function calculate_ratio(a, b) result(ans)
    real(8), intent(in) :: a, b
    real(8) :: ans
    !ans = b ** 4 / (a ** 4 + b ** 4)
    ans = b / (a + b)
  end function calculate_ratio

end module direction_module
