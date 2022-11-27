module new_diagram_module

  use grid_module, only: nlon, nlat, ntrunc, nz, lon, coslatr, height,  &
    gu, gv, gw, gphi, gphi_initial, sphi_old, sphi, longitudes=>lon, latitudes=>lat, wgt
  private
  
  integer(8), dimension(:, :), allocatable, private :: pa, qa, pb, qb, pc, qc, pd, qd
  real(8), dimension(:, :, :), allocatable, private :: midlonA, midlatA, midlonB, midlatB, midlonC, midlatC, midlonD, midlatD
  real(8), dimension(:, :, :), allocatable, private :: &
    gphi_old, dgphi, dgphim, gphim, gphix, gphiy, gphixy, &
    deplon, deplat, depheight, &
    gphiz, gphixz, gphiyz, gphixyz
  complex(8), dimension(:, :, :), allocatable, private :: sphi1
  real(8), allocatable, private :: A(:, :, :), B(:, :, :), C(:, :, :), D(:, :, :)
  real(8), dimension(:, :), allocatable, private :: gumA, gvmA, gumB, gvmB, gumC, gvmC, gumD, gvmD
  real(8), dimension(:, :, :), allocatable, private :: dgphimA, dgphimB, dgphimC, dgphimD, midh

  private :: update
  public :: new_diagram_init, new_diagram_timeint, new_diagram_clean

contains

  subroutine new_diagram_init()
    use time_module, only: deltat
    use interpolate3d_module, only: interpolate_init
    use interpolate1d_module, only: interpolate1d_init
    use interpolate_module, only: interpolate2d_init=>interpolate_init
    use legendre_transform_module, only: legendre_synthesis
    implicit none

    integer(8) :: i, j, k

    allocate(sphi1(0:ntrunc, 0:ntrunc, nz),gphi_old(nlon, nlat, nz))
    allocate(gphim(nlon, nlat, nz),dgphi(nlon, nlat, nz),dgphim(nlon, nlat, nz))
    allocate(deplon(nlon, nlat, nz), deplat(nlon, nlat, nz))
    allocate(depheight(nlon, nlat, nz))
    allocate(gphix(nlon, nlat, nz), gphiy(nlon, nlat, nz), gphixy(nlon, nlat, nz))
    allocate(gphiz(nlon, nlat, nz), gphixz(nlon, nlat, nz), gphiyz(nlon, nlat, nz), gphixyz(nlon, nlat, nz))
    allocate(A(nlon, nlat, nz), B(nlon, nlat, nz), C(nlon, nlat, nz), D(nlon, nlat, nz))
    allocate(midlonA(nlon, nlat, nz), midlatA(nlon, nlat, nz), midlonB(nlon, nlat, nz), midlatB(nlon, nlat, nz))
    allocate(midlonC(nlon, nlat, nz), midlatC(nlon, nlat, nz), midlonD(nlon, nlat, nz), midlatD(nlon, nlat, nz))
    allocate(gumA(nlon, nlat), gvmA(nlon, nlat), gumB(nlon, nlat), gvmB(nlon, nlat))
    allocate(gumC(nlon, nlat), gvmC(nlon, nlat), gumD(nlon, nlat), gvmD(nlon, nlat))
    allocate(dgphimA(nlon, nlat, nz), dgphimB(nlon, nlat, nz), dgphimC(nlon, nlat, nz), dgphimD(nlon, nlat, nz))
    allocate(pa(nlon, nlat), qa(nlon, nlat), pb(nlon, nlat), qb(nlon, nlat))
    allocate(pc(nlon, nlat), qc(nlon, nlat), pd(nlon, nlat), qd(nlon, nlat), midh(nlon, nlat, nz))

    call interpolate_init(gphi)
    call interpolate1d_init(gphi(1, 1, :))
    call interpolate2d_init(gphi(:, :, 1))

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

  end subroutine new_diagram_init

  subroutine new_diagram_clean()
    use interpolate3d_module, only: interpolate_clean
    implicit none

    deallocate(sphi1,gphi_old,gphim,dgphi,dgphim,deplon,deplat)
    call interpolate_clean()

  end subroutine new_diagram_clean

  subroutine new_diagram_timeint()
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
    
  end subroutine new_diagram_timeint

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
    use interpolate_module, only: interpolate_dist_ratio

    implicit none

    integer(8) :: i, j, k, m
    real(8), intent(in) :: t, dt
    real(8) :: tmppres, ans
    real(8), allocatable :: zdot(:, :, :)
    integer(8), allocatable :: id(:, :, :)
    real(8), parameter :: g = 9.80616d0
    real(8) :: tmp

    allocate(zdot(nlon, nlat, nz), id(nlon, nlat, nz))

    select case(case)
      case('hadley')
        call uv_hadley(t, lon, latitudes, height, gu, gv, gw)
      case('div')
        call uv_div(t, lon, latitudes, height, gu, gv, gw)
      case default
        print *, "No matching initial field"
      stop
    end select
    call find_points(gu, gv, gw, t, 0.5d0*dt, midlonA, midlatA, deplon, deplat, depheight)
    call set_niuv(dt)

    do k = 1, nz
      do i = 1, nlon
        do j = 1, nlat
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

    !!!! ここからnislと違う !!!!!!!!!!!!!!

    do k = 1, nz
      do i = 1, nlon
        do j = 1, nlat
          call interpolate_dist_ratio(deplon(i, j, k), deplat(i, j, k), A(i, j, k), B(i, j, k), C(i, j, k), D(i, j, k))
        end do
      end do
    end do

    ! dF/dlon
    do k = 1, nz
      call legendre_synthesis_dlon(sphi(:, :, k), dgphi(:, :, k))
    enddo
    call tricubic_interpolate_set(dgphi, gphix, gphiy, gphiz, gphixy, gphixz, gphiyz, gphixyz)
    call interpolate_set(dgphi)
    call interpolate_setd(gphix, gphiy, gphiz, gphixy, gphixz, gphiyz, gphixyz)
    gphim(:, :, :) = 0.0d0
    do j = 1, nlat
      do i = 1, nlon
        do k = 1, nz
          tmp = transorm_height_to_pressure(midh(i, j, k))
          call interpolate_tricubic(midlonA(i, j, k), midlatA(i, j, k), tmp, dgphimA(i, j, k))
          call interpolate_tricubic(midlonB(i, j, k), midlatB(i, j, k), tmp, dgphimB(i, j, k))
          call interpolate_tricubic(midlonC(i, j, k), midlatC(i, j, k), tmp, dgphimC(i, j, k))
          call interpolate_tricubic(midlonD(i, j, k), midlatD(i, j, k), tmp, dgphimD(i, j, k))

          gphim(i, j, k) = gphim(i, j, k) + A(i, j, k) * gumA(i, j) * dgphimA(i, j, k) / cos(latitudes(j))
         ! gphim(i, j, k) = gphim(i, j, k) + B(i, j, k) * gumB(i, j) * dgphimB(i, j, k) / cos(latitudes(j))
         ! gphim(i, j, k) = gphim(i, j, k) + C(i, j, k) * gumC(i, j) * dgphimC(i, j, k) / cos(latitudes(j))
         ! gphim(i, j, k) = gphim(i, j, k) + D(i, j, k) * gumD(i, j) * dgphimD(i, j, k) / cos(latitudes(j))
        enddo
      enddo
    enddo

    ! cos(lat)dF/dlat
    do k = 1, nz
      call legendre_synthesis_dlat(sphi(:, :, k), dgphi(:, :, k))
    enddo
    call tricubic_interpolate_set(dgphi, gphix, gphiy, gphiz, gphixy, gphixz, gphiyz, gphixyz)
    call interpolate_set(dgphi)
    call interpolate_setd(gphix, gphiy, gphiz, gphixy, gphixz, gphiyz, gphixyz)
    do j = 1, nlat
      do i = 1, nlon
        do k = 1, nz
          tmp = transorm_height_to_pressure(midh(i, j, k))
          call interpolate_tricubic(midlonA(i, j, k), midlatA(i, j, k), tmp, dgphimA(i, j, k))
          call interpolate_tricubic(midlonB(i, j, k), midlatB(i, j, k), tmp, dgphimB(i, j, k))
          call interpolate_tricubic(midlonC(i, j, k), midlatC(i, j, k), tmp, dgphimC(i, j, k))
          call interpolate_tricubic(midlonD(i, j, k), midlatD(i, j, k), tmp, dgphimD(i, j, k))

          !gphim(i, j, k) = gphim(i, j, k) + A(i, j, k) * gvmA(i, j) * dgphimA(i, j, k) / cos(latitudes(j))
         !! gphim(i, j, k) = gphim(i, j, k) + B(i, j, k) * gvmB(i, j) * dgphimB(i, j, k) / cos(latitudes(j))
         ! gphim(i, j, k) = gphim(i, j, k) + C(i, j, k) * gvmC(i, j) * dgphimC(i, j, k) / cos(latitudes(j))
         ! gphim(i, j, k) = gphim(i, j, k) + D(i, j, k) * gvmD(i, j) * dgphimD(i, j, k) / cos(latitudes(j))
        enddo
      enddo
    enddo

    gphi = gphi + dt * gphim

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

  subroutine set_niuv(dt)
    use math_module, only: math_pi, pi2=>math_pi2
    use sphere_module, only: xyz2uv, lonlat2xyz
    use interpolate_module, only: find_stencil_
    use planet_module, only: transorm_height_to_pressure
    implicit none

    real(8), intent(in) :: dt

    integer(8) :: i, j, k
    real(8) :: dlonr, pr
    integer(8), dimension(4) :: tmp1, tmp2

    dlonr = 0.5d0 * nlon / math_pi
    call set_velocity_zero()
    do k = 1, nz
      do j = 1, nlat
        do i = 1, nlon
          ! find grid points near departure points
          call find_stencil_(deplon(i, j, 1), deplat(i, j, 1), tmp1, tmp2)
          pa(i, j) = tmp1(1); qa(i, j) = tmp2(1)
          pb(i, j) = tmp1(2); qb(i, j) = tmp2(2)
          pc(i, j) = tmp1(3); qc(i, j) = tmp2(3)
          pd(i, j) = tmp1(4); qd(i, j) = tmp2(4)

          pr = transorm_height_to_pressure(midh(i, j, k))

          call calc_niuv(dt, pa(i, j), qa(i, j), longitudes(i), latitudes(j), midlonA(i, j, k), &
           midlatA(i, j, k), pr, gumA(i, j), gvmA(i, j))
          call calc_niuv(dt, pb(i, j), qb(i, j), longitudes(i), latitudes(j), midlonB(i, j, k), &
            midlatB(i, j, k), pr, gumB(i, j), gvmB(i, j))
          call calc_niuv(dt, pc(i, j), qc(i, j), longitudes(i), latitudes(j), midlonC(i, j, k), &
            midlatC(i, j, k), pr, gumC(i, j), gvmC(i, j))
          call calc_niuv(dt, pd(i, j), qd(i, j), longitudes(i), latitudes(j), midlonD(i, j, k), &
            midlatD(i, j, k), pr, gumD(i, j), gvmD(i, j))
        end do
      end do
    enddo

  end subroutine  set_niuv

  subroutine set_velocity_zero
    implicit none
    gumA(:, :) = 0.0d0
    gumB(:, :) = 0.0d0
    gumC(:, :) = 0.0d0
    gumD(:, :) = 0.0d0
    gvmA(:, :) = 0.0d0
    gvmB(:, :) = 0.0d0
    gvmC(:, :) = 0.0d0
    gvmD(:, :) = 0.0d0
  end subroutine set_velocity_zero

  ! Ritchie1987 式(45)のu^* - Uをgumに詰める(gvmも)
  subroutine calc_niuv(dt, p, q, lon1, lat, midlon, midlat, midp, gum, gvm)
    use grid_module, only: latitudes => lat, longitudes => lon
    use interpolate_module, only: interpolate_bilinearuv
    use math_module, only: math_pi, pi2=>math_pi2
    use sphere_module, only: xyz2uv, lonlat2xyz
    use uv_module, only: calc_ua, calc_ud, calc_va
    implicit none
    real(8), intent(in) :: dt
    integer(8), intent(in) :: p, q
    real(8), intent(in) :: lon1, lat, midp
    real(8), intent(out) :: gum, gvm
    real(8), intent(out) :: midlon, midlat

    real(8) :: xg, yg, zg, xr, yr, zr, xm, ym, zm, xdot, ydot, zdot, lon_grid, lat_grid, u, v, b1


    lon_grid = longitudes(p)
    lat_grid = latitudes(q)
    call lonlat2xyz(lon_grid, lat_grid, xr, yr, zr)
    ! arrival points
    call lonlat2xyz(lon1, lat, xg, yg, zg)

    b1 = 1.0d0 / sqrt( 2.0d0 * (1.0d0 + (xg*xr + yg*yr + zg*zr)) ) ! Ritchie1987 式(44)
    xm = b1 * (xg + xr)
    ym = b1 * (yg + yr)
    zm = b1 * (zg + zr)
    midlon = modulo(atan2(ym, xm) + pi2, pi2)
    midlat = asin(zm)

    xdot = (xg - xr) / dt
    ydot = (yg - yr) / dt
    zdot = (zg - zr) / dt
    call xyz2uv(xdot, ydot, zdot, midlon, midlat, u, v)  !Richie1987式(49)
    gum = u
    gvm = v

    call interpolate_bilinearuv(midlon, midlat, u, v)
    u = calc_ua(midlon, midlat, midp) + calc_va(midlon, midlat, midp)
    v = calc_va(midlon, midlat, midp)

    gum = gum - u
    gvm = gvm - v

  end subroutine calc_niuv

  subroutine tricubic_interpolate_set(f, fx, fy, fz, fxy, fxz, fyz, fxyz)
    use legendre_transform_module, only: legendre_analysis, legendre_synthesis_dlat, legendre_synthesis_dlon, &
      legendre_synthesis_dlonlat
    implicit none
    real(8), intent(in) :: f(:, :, :)
    real(8), intent(out), dimension(:, :, :) :: fx, fy, fz, fxy, fxz, fyz, fxyz

    integer(8) :: j, k

    do k = 1, nz
      call legendre_analysis(f(:, :, k), sphi1(:, :, k))
      call legendre_synthesis_dlon(sphi1(:, :, k), fx(:, :, k))
      call legendre_synthesis_dlat(sphi1(:, :, k), fy(:, :, k))
      call legendre_synthesis_dlonlat(sphi1(:, :, k), fxy(:, :, k))
      do j = 1, nlat
        fy(:, j, k) = fy(:, j, k) / cos(latitudes(j))
        fxy(:, j, k) = fxy(:, j, k) / cos(latitudes(j))
      end do
    end do

    do k = 2, nz-1
      fz(:, :, k) = (f(:, :, k + 1) - f(:, :, k - 1)) / (height(k+1) - height(k-1))
      fxz(:, :, k) = (fx(:, :, k + 1) - fx(:, :, k - 1)) / (height(k+1) - height(k-1))
      fyz(:, :, k) = (fy(:, :, k + 1) - fy(:, :, k - 1)) / (height(k+1) - height(k-1))
      fxyz(:, :, k) = (fxy(:, :, k + 1) - fxy(:, :, k - 1)) / (height(k+1) - height(k-1))
    end do

    fz(:, :, 1) = (f(:, :, 2) - f(:, :, 1)) / (height(2) - height(1))
    fxz(:, :, 1) = (fx(:, :, 2) - fx(:, :, 1)) / (height(2) - height(1))
    fyz(:, :, 1) = (fy(:, :, 2) - fy(:, :, 1)) / (height(2) - height(1))
    fxyz(:, :, 1) = (fxy(:, :, 2) - fxy(:, :, 1)) / (height(2) - height(1))

    fz(:, :, nz) = (f(:, :, nz) - f(:, :, nz - 1)) / (height(nz) - height(nz-1))
    fxz(:, :, nz) = (fx(:, :, nz) - fx(:, :, nz - 1)) / (height(nz) - height(nz-1))
    fyz(:, :, nz) = (fy(:, :, nz) - fy(:, :, nz - 1)) / (height(nz) - height(nz-1))
    fxyz(:, :, nz) = (fxy(:, :, nz) - fxy(:, :, nz - 1)) / (height(nz) - height(nz-1))

  end subroutine tricubic_interpolate_set

end module new_diagram_module
