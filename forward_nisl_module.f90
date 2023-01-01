module forward_nisl_module

  use grid_module, only: nlon, nlat, ntrunc, nz, lon, coslatr, height, dh, &
    gu, gv, gw, gphi, gphi_initial, longitudes=>lon, latitudes=>lat, wgt, sphi_old
  private

  real(8), dimension(:, :, :), allocatable, private :: &
    gphi_old, dgphi, dgphim, gphim, gphix, gphiy, gphixy, &
    midlon, midlat, midh, deplon, deplat, gum, gvm, gwm, depheight, &
    gphiz, gphi1, gphixz, gphiyz, gphixyz
  integer(8), allocatable, private :: p(:, :, :), q(:, :, :), r(:, :, :)

  private :: update
  public :: forward_nisl_init, forward_nisl_timeint, forward_nisl_clean

contains

  subroutine forward_nisl_init()
    use time_module, only: deltat
    use interpolate2d_module, only: interpolate2d_init
    use interpolate3d_module, only: interpolate3d_init=>interpolate_init
    use interpolate_module, only: interpolate_init
    use legendre_transform_module, only: legendre_synthesis
    implicit none

    integer(8) :: i, j

    allocate(gphi_old(nlon, nlat, nz), gphi1(nlon, nlat, nz))
    allocate(gphim(nlon, nlat, nz),dgphi(nlon, nlat, nz),dgphim(nlon, nlat, nz))
    allocate(midlon(nlon, nlat, nz), midlat(nlon, nlat, nz), midh(nlon, nlat, nz))
    allocate(deplon(nlon, nlat, nz), deplat(nlon, nlat, nz))
    allocate(gum(nlon, nlat, nz), gvm(nlon, nlat, nz), depheight(nlon, nlat, nz))
    allocate(gphix(nlon, nlat, nz), gphiy(nlon, nlat, nz), gphixy(nlon, nlat, nz))
    allocate(gphiz(nlon, nlat, nz), gwm(nlon, nlat, nz))
    allocate(p(nlon, nlat, nz), q(nlon, nlat, nz), r(nlon, nlat, nz))
    allocate(gphixz(nlon, nlat, nz), gphiyz(nlon, nlat, nz), gphixyz(nlon, nlat, nz))

    call interpolate3d_init(gphi)
    call interpolate2d_init(gphi)
    call interpolate_init(gphi(:,:,1))

    gphi_old = gphi
    gphi1 = gphi

    open(11, file="animation.txt")
    do i = 1, nlat
      do j = 1, nz
          write(11,*) latitudes(i), height(j), gphi(nlon/2, i, j)
      end do
    end do

  end subroutine forward_nisl_init

  subroutine forward_nisl_clean()
    use interpolate2d_module, only: interpolate2d_clean
    implicit none

    deallocate(gphi_old,gphim,dgphi,dgphim,gum,gvm, gwm, &
      midlon,midlat,deplon,deplat, p, q)
    call interpolate2d_clean()

  end subroutine forward_nisl_clean

  subroutine forward_nisl_timeint()
    use time_module, only: nstep, deltat, hstep
    use legendre_transform_module, only: legendre_synthesis
    implicit none

    integer(8) :: i, j, k

    call update(0.5d0*deltat, deltat)
    write(*, *) 'step = 0 ', "maxval = ", maxval(gphi), 'minval = ', minval(gphi)
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
    !do i = 1, nlat
    !  do j = 1, nz
    !    write(10,*) latitudes(i), height(j), gphi(nlon/2, i, j)
    !  enddo
    !enddo
    do i = 1, nlon
      do j = 1, nlat
        write(10, *) lon(i), latitudes(j), gphi(i, j, 25)
      enddo
    enddo
    close(10)

  end subroutine forward_nisl_timeint

  subroutine update(t, dt)
    use uv_module, only: uv_div
    use upstream_forward_module, only: find_forward_points
    use legendre_transform_module, only: legendre_analysis, legendre_synthesis, &
        legendre_synthesis_dlon, legendre_synthesis_dlat, legendre_synthesis_dlonlat
    use interpolate2d_module, only: interpolate2d_set, interpolate2d_setd, interpolate2d_bicubic
    use spline_interpolate_module, only: interpolate_spline_1index
    use interpolate3d_module, only: interpolate_set, interpolate_setd, interpolate_tricubic

    implicit none

    integer(8) :: i, j, k
    real(8), intent(in) :: t, dt
    real(8) :: ans
    real(8), allocatable :: f(:, :, :)

    allocate(f(nlon, nlat, nz))

    call find_forward_points(t-0.5d0*dt, dt, deplon, deplat, depheight)

    call find_points()

    do i = 1, nlon
      do j = 1, nlat
        do k = 1, nz
          gphi(i, j, k) = gphi_old(p(i,j,k), q(i,j,k), r(i,j,k))
        enddo
      enddo
    enddo

    call mid_points(t, dt)

    call lon_derivative(gphi1, f)
    do j = 1, nlat
      f(:, j, :) = f(:, j, :) / cos(latitudes(j))
    enddo
    call fd_derivative(f)
    call interpolate_set(f)
    call interpolate_setd(gphix, gphiy, gphiz, gphixy, gphixz, gphiyz, gphixyz)
    do i = 1, nlon
      do j = 1, nlat
        do k = 1, nz
          call interpolate_tricubic(midlon(i,j,k), midlat(i,j,k), midh(i,j,k), ans)
          gphi(i, j, k) = gphi(i, j, k) + gum(i, j, k) * ans * dt
        enddo
      enddo
    enddo

    call lat_derivative(gphi1, f)
    call fd_derivative(f)
    call interpolate_set(f)
    call interpolate_setd(gphix, gphiy, gphiz, gphixy, gphixz, gphiyz, gphixyz)
    do i = 1, nlon
      do j = 1, nlat
        do k = 1, nz
          call interpolate_tricubic(midlon(i,j,k), midlat(i,j,k), midh(i,j,k), ans)
          gphi(i, j, k) = gphi(i, j, k) + gvm(i, j, k) * ans * dt
        enddo
      enddo
    enddo

    call vertical_derivative(gphi1, f)
    call fd_derivative(f)
    call interpolate_set(f)
    call interpolate_setd(gphix, gphiy, gphiz, gphixy, gphixz, gphiyz, gphixyz)
    do i = 1, nlon
      do j = 1, nlat
        do k = 1, nz
          call interpolate_tricubic(midlon(i,j,k), midlat(i,j,k), midh(i,j,k), ans)
          gphi(i, j, k) = gphi(i, j, k) + gwm(i, j, k) * ans * dt
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

  subroutine mid_points(t, dt)
    use sphere_module, only: lonlat2xyz, xyz2uv
    use math_module, only: pi2=>math_pi2
    use time_module, only: case
    use uv_hadley_module, only: calc_hadleyu=>calc_u, calc_hadleyv=>calc_v, calc_hadleyw=>calc_w
    use uv_module, only: calc_ua, calc_ud, calc_va, calc_w
    implicit none
    real(8), intent(in) :: t, dt
    integer(8) :: i, j, k
    real(8) :: xg, yg, zg, xr, yr, zr, xm, ym, zm
    real(8) :: xdot, ydot, zdot
    real(8) :: u1, v1, u2, v2, w1, w2

    do k = 1, nz
      do i = 1, nlon
        do j = 1, nlat
          call lonlat2xyz(lon(i), latitudes(j), xg, yg, zg)
          call lonlat2xyz(lon(p(i,j,k)), latitudes(q(i,j,k)), xr, yr, zr)
          xm = (xg + xr) * 0.5d0
          ym = (yg + yr) * 0.5d0
          zm = (zg + zr) * 0.5d0
          midlon(i, j, k) = modulo(atan2(ym, xm) + pi2, pi2)
          midlat(i, j, k) = asin(zm / sqrt(xm**2 + ym**2 + zm**2))
          midh(i, j, k) = (height(k) + height(r(i, j, k))) / 2.0d0
          call check_height(midh(i, j, k))

          xdot = (xg - xr) / dt
          ydot = (yg - yr) / dt
          zdot = (zg - zr) / dt
          call xyz2uv(xdot, ydot, zdot, midlon(i,j,k), midlat(i,j,k), u1, v1)  !Richie1987式(49)
          w1 = (height(k) - height(r(i, j, k))) / dt

          select case(case)
          case('hadley')
            u2 = calc_hadleyu(midlat(i,j,k))
            v2 = calc_hadleyv(midlat(i,j,k), midh(i,j,k), t)
            w2 = calc_hadleyw(midlat(i,j,k), midh(i,j,k), t)
          case('div')
            u2 = calc_ua(midlon(i,j,k), midlat(i,j,k), t) &
              + calc_ud(midlon(i,j,k), midlat(i,j,k), midh(i,j,k), t)
            v2 = calc_va(midlon(i,j,k), midlat(i,j,k), t)
            w2 = calc_w(midlon(i,j,k), midlat(i,j,k), midh(i,j,k), t)
          end select
          gum(i, j, k) = u1 - u2
          gvm(i, j, k) = v1 - v2
          gwm(i, j, k) = w1 - w2
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
    call vertical_derivative(gphix, gphixz)
    call vertical_derivative(gphiy, gphiyz)
    call vertical_derivative(gphixy, gphixyz)
  end subroutine fd_derivative

  subroutine vertical_derivative(f, fz)
    implicit none
    integer(8) :: k
    real(8), intent(in) :: f(nlon, nlat, nz)
    real(8), intent(out) :: fz(nlon, nlat, nz)
    do k = 2, nz-1
      fz(:, :, k) = (f(:, :, k + 1) - f(:, :, k - 1)) / (height(k+1) - height(k-1))
    end do
    fz(:, :, 1) = (f(:, :, 2) - f(:, :, 1)) / (height(2) - height(1))
    fz(:, :, nz) = (f(:, :, nz) - f(:, :, nz - 1)) / (height(nz) - height(nz-1))
  end subroutine vertical_derivative

  subroutine lon_derivative(f, fx)
    use math_module, only: pir=>math_pir
    implicit none
    real(8), intent(in) :: f(nlon, nlat, nz)
    real(8), intent(out) :: fx(nlon, nlat, nz)
    integer(8) :: i, k
    real(8) :: dlonr

    do k = 1, nz
      dlonr = 0.25d0*nlon*pir
      fx(1,:,k) = dlonr * (f(2,:,k) - f(nlon,:,k))
      fx(nlon, :, k) = dlonr * (f(1,:,k) - f(nlon-1,:,k))
      do i=2, nlon-1
        fx(i,:,k) = dlonr*(f(i+1,:,k) - f(i-1,:,k))
      end do
    end do
  end subroutine lon_derivative

  subroutine lat_derivative(f, fy)
    use math_module, only: pir=>math_pir, pih=>math_pih
    implicit none
    real(8), intent(in) :: f(nlon, nlat, nz)
    real(8), intent(out) :: fy(nlon, nlat, nz)
    integer(8) :: j, k
    real(8) :: gphitmp(nlon), eps

    do k = 1, nz
      eps = pih-latitudes(1)
      gphitmp = cshift(f(:,1,k),nlon/2)
      fy(:,1,k) = (gphitmp-f(:,2,k))/(pih+eps-latitudes(2))
      gphitmp = cshift(f(:,nlat,k),nlon/2)
      fy(:,nlat,k) = (gphitmp-f(:,nlat-1,k))/(-pih-eps-latitudes(nlat-1))
      do j=2, nlat-1
        fy(:,j,k) = (f(:,j+1,k)-f(:,j-1,k))/(latitudes(j+1)-latitudes(j-1))
      end do
    end do
  end subroutine lat_derivative

  subroutine find_points
    use math_module, only: pir=>math_pir, pi=>math_pi
    use interpolate_module, only: find_stencil_
    use grid_module, only: pole_regrid
    implicit none
    integer(8) :: i, j, k, i1, j1, k1, it, jt
    real(8) :: dlonr
    integer(8) :: p1, q1, r1
    real(8), allocatable :: use_dist(:,:,:)
    real(8) :: min_dist, tmp_dist

    allocate(use_dist(nlon, nlat, nz))
    dlonr = 0.5d0*nlon*pir
    use_dist(:, :, :) = 1.0d9

    p(:, :, :) = -100
    q(:, :, :) = -100
    r(:, :, :) = -100

    do j = 1, nlat
      do i = 1, nlon
        do k = 1, nz
          p1 = int(anint( deplon(i, j, k) * dlonr + 1.0d0 ))
          if ( p1 > nlon ) then
            p1 = p1 - nlon
          end if
          q1 = int(anint( 0.5d0 * (nlat + 1.0d0 - (2.0d0*dble(nlat)+1.0d0)*deplat(i,j,k) / pi) ))  !latitudesは大きい順で詰められているので注意
          call check_height(depheight(i, j, k))
          r1 = int(anint(depheight(i, j, k) / dh)) + 1
          p(p1, q1, r1) = i
          q(p1, q1, r1) = j
          r(p1, q1, r1) = k
        enddo
      end do
    end do

    do i = 1, nlon
      do j = 1, nlat
        do k = 1, nz
          if (p(i, j, k) /= -100) cycle
          min_dist = 1.0d9
          do i1 = i-2, i+2
            do j1 = j-2, j+2
              do k1 = max(k-1,1), min(k+1, nz)
                it = i1; jt = j1
                call pole_regrid(it, jt)
                tmp_dist = sphere_dist(lon(i), latitudes(j), height(k), &
                  deplon(it,jt,k1), deplat(it,jt,k1), depheight(it,jt,k1))
                  if (min_dist > tmp_dist) then
                    min_dist = tmp_dist
                    p(i, j, k) = it
                    q(i, j, k) = jt
                    r(i, j, k) = k1
                  endif
              enddo
            enddo
          enddo
        enddo
      enddo
    enddo

  end subroutine find_points

  function cartesian_dist(x1, y1, z1, x2, y2, z2) result(l)
    implicit none
    real(8), intent(in) :: x1, y1, z1, x2, y2, z2
    real(8) :: l

    l = sqrt((x1 - x2)**2 + (y1 - y2)**2 + (z1 - z2)**2)
  end function cartesian_dist

  function sphere_dist(lon1, lat1, h1, lon2, lat2, h2) result(l)
    use planet_module, only: planet_radius
    implicit none
    real(8) :: lon1, lat1, h1, lon2, lat2, h2
    real(8) :: l
    real(8) :: x1, y1, z1, x2, y2, z2

    x1 = (planet_radius + h1) * cos(lon1) * cos(lat1)
    y1 = (planet_radius + h1) * sin(lon1) * cos(lat1)
    z1 = (planet_radius + h1) * sin(lat1)
    x2 = (planet_radius + h2) * cos(lon2) * cos(lat2)
    y2 = (planet_radius + h2) * sin(lon2) * cos(lat2)
    z2 = (planet_radius + h2) * sin(lat2)

    l = cartesian_dist(x1, y1, z1, x2, y2, z2)
  end function sphere_dist
end module forward_nisl_module
