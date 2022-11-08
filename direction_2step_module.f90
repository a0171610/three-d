module direction_2step_module

  use grid_module, only: nlon, nlat, ntrunc, dlat, dlat4, &
    gu, gv, gphi, gphi_initial, sphi_old, sphi, longitudes=>lon, latitudes=>lat, wgt
  use time_module, only: velocity, imethod, conserve
  private
  
  integer(8), dimension(:, :), allocatable, private :: pa, qa, pb, qb, pc, qc, pd, qd
  real(8), allocatable, private :: wa(:, :), wb(:, :), wc(:, :), wd(:, :) 
  real(8), dimension(:,:), allocatable, private :: &
    gphi_old, dgphi, gphim, gphix, gphiy, gphixy, deplon, deplat
  real(8), dimension(:, :), allocatable, private :: midlonA, midlatA, midlonB, midlatB, midlonC, midlatC, midlonD, midlatD
  real(8), dimension(:, :), allocatable, private :: gumA, gvmA, gumB, gvmB, gumC, gvmC, gumD, gvmD
  real(8), dimension(:, :), allocatable, private :: dgphimA, dgphimB, dgphimC, dgphimD
  real(8), private :: dlon, right(nlon, nlat)

  private :: update
  public :: direction_2step_init, direction_2step_timeint, direction_2step_clean

contains

  subroutine direction_2step_init()
    use interpolate_module, only: interpolate_init
    use legendre_transform_module, only: legendre_synthesis
    implicit none

    integer(8) :: i, j

    allocate(gphi_old(nlon, nlat))
    allocate(gphim(nlon, nlat),dgphi(nlon, nlat))
    allocate(midlonA(nlon, nlat), midlatA(nlon, nlat), midlonB(nlon, nlat), midlatB(nlon, nlat))
    allocate(midlonC(nlon, nlat), midlatC(nlon, nlat), midlonD(nlon, nlat), midlatD(nlon, nlat))
    allocate(deplon(nlon, nlat), deplat(nlon, nlat), pa(nlon, nlat), qa(nlon, nlat), pb(nlon, nlat), qb(nlon, nlat))
    allocate(pc(nlon, nlat), qc(nlon, nlat), pd(nlon, nlat), qd(nlon, nlat))
    allocate(wa(nlon, nlat), wb(nlon, nlat), wc(nlon, nlat), wd(nlon, nlat))
    allocate(guma(nlon, nlat), gvma(nlon, nlat), gumb(nlon, nlat), gvmb(nlon, nlat))
    allocate(gumc(nlon, nlat), gvmc(nlon, nlat), gumd(nlon, nlat), gvmd(nlon, nlat))
    allocate(dgphimA(nlon, nlat), dgphimB(nlon, nlat), dgphimC(nlon, nlat), dgphimD(nlon, nlat))
    allocate(gphix(nlon, nlat), gphiy(nlon, nlat), gphixy(nlon, nlat))

    call interpolate_init(gphi)

    call legendre_synthesis(sphi_old, gphi_old)
    gphi(:, :) = gphi_old(:, :)
    dlon = longitudes(3) - longitudes(1)

    open(11, file="animation.txt")
    do i = 1, nlon
      do j = 1, nlat
          write(11,*) longitudes(i), latitudes(j), gphi(i, j)
      end do
    end do

  end subroutine direction_2step_init

  subroutine direction_2step_clean()
    use interpolate_module, only: interpolate_clean
    implicit none

    deallocate(gphi_old, gphim, dgphi, deplon, deplat)
    deallocate(midlonA, midlatA, midlonB, midlatB, midlonC, midlatC, midlonD, midlatD)
    deallocate(guma, gvma, gumb, gvmb, gumc, gvmc, gumd, gvmd)
    deallocate(wa, wb, wc, wd, pa, qa, pb, qb, pc, qc, pd, qd)
    deallocate(dgphimA, dgphimB, dgphimC, dgphimD)
    call interpolate_clean()

  end subroutine direction_2step_clean

  subroutine direction_2step_timeint()
    use time_module, only: nstep, deltat, hstep, field
    implicit none

    integer(8) :: i, j, k

    do i = 1, nstep
      call update((i-0.5d0)*deltat, deltat)
      write(*, *) 'step = ', i, "maxval = ", maxval(gphi), 'minval = ', minval(gphi)
      if ( mod(i, hstep) == 0 ) then
        do j = 1, nlon
            do k = 1, nlat
              write(11,*) longitudes(j), latitudes(k), gphi(j, k)
            end do
        end do
      endif
      if (i == nstep / 2 .and. field == "cbell2") then
        open(10, file="log_cbell.txt")
        do j = 1, nlon
          do k = 1, nlat
            write(10,*) gphi(j, k)
          enddo
        enddo
        close(10)
      endif
      if (i == nstep / 2 .and. field == "ccbell2") then
        open(10, file="log_ccbell.txt")
        do j = 1, nlon
          do k = 1, nlat
            write(10,*) wgt(k), gphi(j, k)
          enddo
        enddo
        close(10)
      endif
    end do
    close(11)
    
  end subroutine direction_2step_timeint

  ! dt = leapfrog法の+と-の時刻差, t=中央の時刻
  subroutine update(t, dt)
    use uv_module, only: uv_nodiv, uv_div, uv_sbody
    use upstream_module, only: find_points
    use legendre_transform_module, only: legendre_analysis, legendre_synthesis, &
        legendre_synthesis_dlon, legendre_synthesis_dlat
    use interpolate_module, only: interpolate_set, interpolate_setd, find_stencil_
    use interpolate_module, only: interpolate_bicubic, interpolate_bilinear, interpolate_bilinear_ratio
    use interpolate_module, only: interpolate_dist, interpolate_dist_ratio, interpolate_setuv
    use math_module, only: pih=>math_pih
    implicit none

    integer(8) :: i, j, m
    real(8), intent(in) :: t, dt
    real(8), allocatable :: x(:), b(:)
    allocate(x(nlon * nlat), b(nlon * nlat))

    select case(velocity)
    case("nodiv")
      call uv_nodiv(t,longitudes,latitudes,gu,gv)
    case("div")
      call uv_div(t,longitudes,latitudes,gu,gv)
    end select
    call find_points(gu, gv, 0.5d0*dt, midlonA, midlatA, deplon, deplat)

    do i = 1, nlon
      do j = 1, nlat
        ! AとDの経度番号がi, BとCの経度番号がi+1, AとBの緯度番号がj, CとDの緯度番号がj+1
        call interpolate_dist_ratio(deplon(i, j), deplat(i, j), wa(i, j), wb(i, j), wc(i, j), wd(i, j))
      end do
    end do

    call interpolate_setuv(gu, gv)
    call set_niuv(t, dt)

    call legendre_synthesis(sphi_old, gphi_old)

    ! まずはgphiにbilinear法で求めた値を詰めていく
    call interpolate_set(gphi_old)
    do j = 1, nlat
      do i = 1, nlon
        call interpolate_dist(deplon(i, j), deplat(i, j), gphi(i, j))
      end do
    end do

    ! dF/dlon
    call calculate_diff(gphi_old, gphix, gphiy)
    call calculate_scale_diff(dt)

    do j = 1, nlat
      do i = 1, nlon
        gphi(i, j) = gphi(i, j) + 0.5d0 * dt * wa(i, j) * gumA(i, j) * gphix(pa(i, j), qa(i, j)) / cos(latitudes(j))
        gphi(i, j) = gphi(i, j) + 0.5d0 * dt * wb(i, j) * gumB(i, j) * gphix(pb(i, j), qb(i, j)) / cos(latitudes(j))
        gphi(i, j) = gphi(i, j) + 0.5d0 * dt * wc(i, j) * gumC(i, j) * gphix(pc(i, j), qc(i, j)) / cos(latitudes(j))
        gphi(i, j) = gphi(i, j) + 0.5d0 * dt * wd(i, j) * gumD(i, j) * gphix(pd(i, j), qd(i, j)) / cos(latitudes(j))
      enddo
    enddo

    ! cos(lat)dF/dlat
    do j = 1, nlat
      do i = 1, nlon
        gphi(i, j) = gphi(i, j) + 0.5d0 * dt * wa(i, j) * gvmA(i, j) * gphiy(pa(i, j), qa(i, j))
        gphi(i, j) = gphi(i, j) + 0.5d0 * dt * wb(i, j) * gvmB(i, j) * gphiy(pb(i, j), qb(i, j))
        gphi(i, j) = gphi(i, j) + 0.5d0 * dt * wc(i, j) * gvmC(i, j) * gphiy(pc(i, j), qc(i, j))
        gphi(i, j) = gphi(i, j) + 0.5d0 * dt * wd(i, j) * gvmD(i, j) * gphiy(pd(i, j), qd(i, j))
      enddo
    enddo

    do i = 1, nlon
      do j = 1, nlat
        m = i + (j - 1) * nlon
        b(m) = gphi(i, j) * wgt(j) / right(i, j)
      end do
    end do

    call set_niuv(t+0.5d0*dt, dt)
    call solve_sparse_matrix(dt, b, x)
    do m = 1, nlon*nlat
      call convert_line_to_array(m, i, j)
      gphi(i, j) = x(m)
    end do

    call legendre_analysis(gphi, sphi)
    do m = 0, ntrunc
      sphi_old(m : ntrunc, m) = sphi(m : ntrunc, m)
    enddo

  end subroutine update

  subroutine solve_sparse_matrix(dt, b, x)
    use lsqr_module, only: lsqr_solver_ez
    use sphere_module, only: orthodrome
    use grid_module, only: pole_regrid
    use math_module, only : pih=>math_pih, pir=>math_pir
    implicit none

    integer(8), parameter :: sz = nlat * nlon
    real(8), intent(in) :: dt
    integer(8) :: x1, y1, x2, y2, x3, y3, x4, y4
    real(8), intent(in) :: b(sz)
    real(8), intent(out) :: x(sz)
    type(lsqr_solver_ez) :: solver
    integer :: istop
    integer(8) :: i, j, m, id, row, col
    integer(8), allocatable :: icol(:), irow(:)
    real(8), allocatable :: a(:)
    real(8) :: val
    real(8), allocatable :: gum(:, :), gvm(:, :)
    real(8), allocatable :: scale(:)

    allocate(gum(nlon, nlat), gvm(nlon, nlat), scale(nlon * nlat))

    allocate( icol(sz * 9), irow(sz * 9), a(sz * 9) )

    do i = 1, nlon
      do j = 1, nlat
        if (wa(i, j) > wb(i, j) .and. wa(i, j) > wc(i, j) .and. wa(i, j) > wd(i, j)) then
          wa(i, j) = 1.0d0; wb(i, j) = 0.0d0; wc(i, j) = 0.0d0; wd(i, j) = 0.0d0
        else if (wb(i, j) > wa(i, j) .and. wb(i, j) > wc(i, j) .and. wb(i, j) > wd(i, j)) then
          wa(i, j) = 0.0d0; wb(i, j) = 1.0d0; wc(i, j) = 0.0d0; wd(i, j) = 0.0d0
        else if (wc(i, j) > wa(i, j) .and. wc(i, j) > wb(i, j) .and. wc(i, j) > wd(i, j)) then
          wa(i, j) = 0.0d0; wb(i, j) = 0.0d0; wc(i, j) = 1.0d0; wd(i, j) = 0.0d0
        else
          wa(i, j) = 0.0d0; wb(i, j) = 0.0d0; wc(i, j) = 0.0d0; wd(i, j) = 1.0d0
        endif
      end do
    end do

    gum(:, :) = wa(:, :) * gumA(:, :) + wb(:, :) * gumB(:, :) + wc(:, :) * gumC(:, :) + wd(:, :) * gumD(:, :)
    gvm(:, :) = wa(:, :) * gvmA(:, :) + wb(:, :) * gvmB(:, :) + wc(:, :) * gvmC(:, :) + wd(:, :) * gvmD(:, :)

    scale(:) = 1.0d0

    id = 1
    do i = 1, nlon
      do j = 1, nlat
        row = i + (j-1) * nlon
        x1 = i + 1; y1 = j
        x2 = i - 1; y2 = j
        x3 = i; y3 = j + 1
        x4 = i; y4 = j - 1
        call pole_regrid(x1, y1)
        call pole_regrid(x2, y2)
        call pole_regrid(x3, y3)
        call pole_regrid(x4, y4)

        val = -2.0d0 * gum(i, j) * dt / (3.0d0 * dlon * cos(latitudes(j)))
        col = x1 + (y1 - 1) * nlon
        irow(id) = row
        icol(id) = col
        a(id) = val
        id = id + 1
        scale(col) = scale(col) + val

        val = 2.0d0 * gum(i, j) * dt / (3.0d0 * dlon * cos(latitudes(j)))
        col = x2 + (y2 - 1) * nlon
        irow(id) = row
        icol(id) = col
        a(id) = val
        id = id + 1
        scale(col) = scale(col) + val

        val = -4.0d0 * gvm(i, j) * dt / (3.0d0 * dlat4(j))
        col = x3 + (y3 - 1) * nlon
        irow(id) = row
        icol(id) = col
        a(id) = val
        id = id + 1
        scale(col) = scale(col) + val

        val = 4.0d0 * gvm(i, j) * dt / (3.0d0 * dlat4(j))
        col = x4 + (y4 - 1) * nlon
        irow(id) = row
        icol(id) = col
        a(id) = val
        id = id + 1
        scale(col) = scale(col) + val
      end do
    end do

    do i = 1, nlon
      do j = 1, nlat
        row = i + (j-1) * nlon
        x1 = i + 2; y1 = j
        x2 = i - 2; y2 = j
        x3 = i; y3 = j + 2
        x4 = i; y4 = j - 2
        call pole_regrid(x1, y1)
        call pole_regrid(x2, y2)
        call pole_regrid(x3, y3)
        call pole_regrid(x4, y4)

        val = gum(i, j) * dt / (12.0d0 * dlon * cos(latitudes(j)))
        col = x1 + (y1 - 1) * nlon
        irow(id) = row
        icol(id) = col
        a(id) = val
        id = id + 1
        scale(col) = scale(col) + val

        val = -gum(i, j) * dt / (12.0d0 * dlon * cos(latitudes(j)))
        col = x2 + (y2 - 1) * nlon
        irow(id) = row
        icol(id) = col
        a(id) = val
        id = id + 1
        scale(col) = scale(col) + val

        val = gvm(i, j) * dt / (6.0d0 * dlat4(j))
        col = x3 + (y3 - 1) * nlon
        irow(id) = row
        icol(id) = col
        a(id) = val
        id = id + 1
        scale(col) = scale(col) + val

        val = -gvm(i, j) * dt / (6.0d0 * dlat4(j))
        col = x4 + (y4 - 1) * nlon
        irow(id) = row
        icol(id) = col
        a(id) = val
        id = id + 1
        scale(col) = scale(col) + val
      end do
    end do

    do i = 8*sz+1, 9*sz
      irow(i) = i - 8*sz
      icol(i) = i - 8*sz
      a(i) = 1.0d0
    end do

    do m = 1, nlon*nlat
      call convert_line_to_array(m, i, j)
      scale(m) = scale(m) / wgt(j)
    enddo

    do i = 1, 9*sz
      a(i) = a(i) / scale(icol(i))
    end do

    call solver%initialize(int(sz), int(sz), a, int(irow), int(icol)) ! use defaults for other optional inputs
    call solver%solve(b, 0.0d0, x, istop)       ! solve the linear system
    if (istop /= 1) then
      write(*,*) '収束しませんでした'
    endif
  end subroutine solve_sparse_matrix

  subroutine set_niuv(t, dt)
    use sphere_module, only: xyz2uv, lonlat2xyz
    use interpolate_module, only: find_stencil_
    use upstream_module, only: find_points
    use uv_module, only: uv_div, uv_nodiv
    implicit none

    real(8), intent(in) :: t, dt

    integer(8) :: i, j
    integer(8), dimension(4) :: tmp1, tmp2

    select case(velocity)
    case("nodiv ")
      call uv_nodiv(t,longitudes,latitudes,gu,gv)
    case("div   ")
      call uv_div(t,longitudes,latitudes,gu,gv)
    end select

    call find_points(gu, gv, 0.5d0*dt, midlonA, midlatA, deplon, deplat)

    do j = 1, nlat
      do i = 1, nlon
        ! find grid points near departure points
        call find_stencil_(deplon(i, j), deplat(i, j), tmp1, tmp2)
        pa(i, j) = tmp1(1); qa(i, j) = tmp2(1)
        pb(i, j) = tmp1(2); qb(i, j) = tmp2(2)
        pc(i, j) = tmp1(3); qc(i, j) = tmp2(3)
        pd(i, j) = tmp1(4); qd(i, j) = tmp2(4)

        call calc_niuv(dt, pa(i, j), qa(i, j), longitudes(i), latitudes(j), midlonA(i, j), midlatA(i, j), gumA(i, j), gvmA(i, j))
        call calc_niuv(dt, pb(i, j), qb(i, j), longitudes(i), latitudes(j), midlonB(i, j), midlatB(i, j), gumB(i, j), gvmB(i, j))
        call calc_niuv(dt, pc(i, j), qc(i, j), longitudes(i), latitudes(j), midlonC(i, j), midlatC(i, j), gumC(i, j), gvmC(i, j))
        call calc_niuv(dt, pd(i, j), qd(i, j), longitudes(i), latitudes(j), midlonD(i, j), midlatD(i, j), gumD(i, j), gvmD(i, j))
      end do
    end do
        
  end subroutine  set_niuv

  ! Ritchie1987 式(45)のu^* - Uをgumに詰める(gvmも)
  subroutine calc_niuv(dt, p, q, lon, lat, midlon, midlat, gum, gvm)
    use grid_module, only: latitudes => lat, longitudes => lon
    use interpolate_module, only: interpolate_bilinearuv
    use math_module, only: math_pi, pi2=>math_pi2
    use sphere_module, only: xyz2uv, lonlat2xyz
    implicit none
    real(8), intent(in) :: dt
    integer(8), intent(in) :: p, q
    real(8), intent(in) :: lon, lat
    real(8), intent(out) :: gum, gvm
    real(8), intent(out) :: midlon, midlat

    real(8) :: xg, yg, zg, xr, yr, zr, xm, ym, zm, xdot, ydot, zdot, lon_grid, lat_grid, u, v, b1

    lon_grid = longitudes(p)
    lat_grid = latitudes(q)
    call lonlat2xyz(lon_grid, lat_grid, xr, yr, zr)
    ! arrival points
    call lonlat2xyz(lon, lat, xg, yg, zg)

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
    gum = gum - u
    gvm = gvm - v

  end subroutine calc_niuv

  subroutine convert_line_to_array(m1, i1, j1)
    implicit none
    integer(8), intent(in) :: m1
    integer(8), intent(out) :: i1, j1

    i1 = mod(m1, nlon)
    if (i1 == 0) then
      i1 = nlon
    endif
    j1 = (m1-i1)/nlon + 1
  end subroutine convert_line_to_array

  subroutine calculate_diff(arr, arrx, arry)
    real(8), intent(in) :: arr(nlon, nlat)
    real(8), intent(out) :: arrx(nlon, nlat), arry(nlon, nlat)

    arrx(1,:) = (-arr(3,:) + 8.0d0*arr(2,:) - 8.0d0*arr(nlon,:) + arr(nlon-1,:)) / (6.0d0 * dlon)
    arrx(2,:) = (-arr(4,:) + 8.0d0*arr(3,:) - 8.0d0*arr(1,:) + arr(nlon,:)) / (6.0d0 * dlon)
    arrx(nlon,:) = (-arr(2,:) + 8.0d0*arr(1,:) - 8.0d0*arr(nlon-1,:) + arr(nlon-2, :)) / (6.0d0*dlon)
    arrx(nlon-1,:) = (-arr(1,:) + 8.0d0*arr(nlon,:) - 8.0d0*arr(nlon-2,:) + arr(nlon-3, :)) / (6.0d0*dlon)
    do i=3, nlon-2
      arrx(i,:) = (-arr(i+2,:) + 8.0d0*arr(i+1,:) - 8.0d0*arr(i-1,:) + arr(i-2,:)) / (6.0d0 * dlon)
    end do

    arry(:,1) = (-arr(:,3) + 8.0d0*arr(:,2) - 8.0d0*arr(:,nlat) + arr(:,nlat-1)) / (3.0d0 * dlat4(1))
    arry(:,2) = (-arr(:,4) + 8.0d0*arr(:,3) - 8.0d0*arr(:,1) + arr(:,nlat)) / (3.0d0 * dlat4(2))
    arry(:,nlat) = (-arr(:,2) + 8.0d0*arr(:,1) - 8.0d0*arr(:,nlat-1) + arr(:, nlat-2)) / (3.0d0 * dlat4(nlat))
    arry(:,nlat-1) = (-arr(:,1) + 8.0d0*arr(:,nlat) - 8.0d0*arr(:,nlat-2) + arr(:, nlat-3)) / (3.0d0 * dlat4(nlat-1))
    do j=3, nlat-2
      arry(:,j) = (-arr(:,j+2) + 8.0d0*arr(:,j+1) - 8.0d0*arr(:,j-1) + arr(:,j-2)) / (3.0d0 * dlat4(j))
    end do

  end subroutine calculate_diff

  subroutine calculate_scale_diff(dt)
    use grid_module, only: pole_regrid
    real(8), intent(in) :: dt
    integer(8) :: x, y

    right(:, :) = 1.0d0
    do i = 1, nlon
      do j = 1, nlat
        y = qa(i, j)
        x = pa(i, j) + 2; call pole_regrid(x, y)
        right(x, y) = right(x, y) - dt * wa(i,j) * gumA(i,j) / (12.0d0 * dlon * cos(latitudes(j)))
        x = pa(i, j) + 1; call pole_regrid(x, y)
        right(x, y) = right(x, y) + 8.0d0 * dt * wb(i,j) * gumB(i,j) / (12.0d0 * dlon * cos(latitudes(j)))
        x = pa(i, j) - 1; call pole_regrid(x, y)
        right(x, y) = right(x, y) - 8.0d0 * dt * wc(i,j) * gumC(i,j) / (12.0d0 * dlon * cos(latitudes(j)))
        x = pa(i, j) - 2; call pole_regrid(x, y)
        right(x, y) = right(x, y) + dt * wd(i,j) * gumD(i,j) / (12.0d0 * dlon * cos(latitudes(j)))
      end do
    end do

    do i = 1, nlon
      do j = 1, nlat
        x = pa(i, j); y = qa(i, j) + 2; call pole_regrid(x, y)
        right(x, y) = right(x, y) - dt * wa(i,j) * gvmA(i,j) / (6.0d0 * dlat4(j))
        x = pa(i, j); y = qa(i, j) + 1; call pole_regrid(x, y)
        right(x, y) = right(x, y) + 8.0d0 * dt * wb(i,j) * gvmB(i,j) / (6.0d0 * dlat4(j))
        x = pa(i, j); y = qa(i, j) - 1; call pole_regrid(x, y)
        right(x, y) = right(x, y) - 8.0d0 * dt * wc(i,j) * gvmC(i,j) / (6.0d0 * dlat4(j))
        x = pa(i, j); y = qa(i, j) - 2; call pole_regrid(x, y)
        right(x, y) = right(x, y) + dt * wd(i,j) * gvmD(i,j) / (6.0d0 * dlat4(j))
      end do
    end do
  end subroutine calculate_scale_diff

end module direction_2step_module
