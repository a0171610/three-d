module nisl_module

  use grid_module, only: nlon, nlat, ntrunc, &
    gu, gv, gphi, gphi_initial, sphi_old, sphi, longitudes=>lon, latitudes=>lat, wgt
  use mass_module, only: mass_correct
  use time_module, only: conserve, velocity
  private
  
  integer(8), allocatable, private :: p(:,:), q(:,:)
  real(8), dimension(:,:), allocatable, private :: &
    gphi_old, dgphi, dgphim, gphim, gphix, gphiy, gphixy, &
    midlon, midlat, deplon, deplat, gum, gvm, gmin, gmax, w
  complex(8), dimension(:,:), allocatable, private :: sphi1

  private :: update, bicubic_interpolation_set
  public :: nisl_init, nisl_timeint, nisl_clean

contains

  subroutine nisl_init()
    use time_module, only: deltat
    use interpolate_module, only: interpolate_init
    use legendre_transform_module, only: legendre_synthesis
    implicit none

    integer(8) :: i,j

    allocate(sphi1(0:ntrunc, 0:ntrunc),gphi_old(nlon, nlat))
    allocate(gphim(nlon, nlat),dgphi(nlon, nlat),dgphim(nlon, nlat))
    allocate(midlon(nlon, nlat), midlat(nlon, nlat))
    allocate(deplon(nlon, nlat), deplat(nlon, nlat), p(nlon, nlat), q(nlon, nlat))
    allocate(gum(nlon, nlat), gvm(nlon, nlat))
    allocate(gphix(nlon, nlat), gphiy(nlon, nlat), gphixy(nlon, nlat))

    call interpolate_init(gphi)

    call legendre_synthesis(sphi_old,gphi_old)
    gphi(:, :) = gphi_old(:, :)

    if (conserve) then
      allocate(gmax(nlon,nlat),gmin(nlon,nlat),w(nlon,nlat))
      do j=1, nlat
        w(:,j) = wgt(j)
      end do
      gmin(:,:) = minval(gphi)
      gmax(:,:) = maxval(gphi)
    endif

    open(11, file="animation.txt")
    do i = 1, nlon
      do j = 1, nlat
          write(11,*) longitudes(i), latitudes(j), gphi(i, j)
      end do        
    end do
    call update(0.0d0, deltat)

  end subroutine nisl_init

  subroutine nisl_clean()
    use interpolate_module, only: interpolate_clean
    implicit none

    deallocate(sphi1,gphi_old,gphim,dgphi,dgphim,gum,gvm, &
      midlon,midlat,deplon,deplat,p,q)
    call interpolate_clean()

  end subroutine nisl_clean

  subroutine nisl_timeint()
    use time_module, only: nstep, deltat, hstep, field
    use legendre_transform_module, only: legendre_synthesis
    implicit none

    integer(8) :: i, j, k

    do i = 2, nstep
      call update((i-1)*deltat, 2.0d0*deltat)
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
    
  end subroutine nisl_timeint

  subroutine update(t, dt)
    use uv_module, only: uv_nodiv, uv_div
    use upstream_module, only: find_points
    use legendre_transform_module, only: legendre_analysis, legendre_synthesis, &
        legendre_synthesis_dlon, legendre_synthesis_dlat, legendre_synthesis_dlonlat
    use interpolate_module, only: interpolate_set, interpolate_bilinear, interpolate_setd, interpolate_bicubic
    use interpolate_module, only: interpolate_setuv
    implicit none

    integer(8) :: i, j, m
    real(8), intent(in) :: t, dt

    select case(velocity)
    case("nodiv ")
      call uv_nodiv(t,longitudes,latitudes,gu,gv)
    case("div   ")
      call uv_div(t,longitudes,latitudes,gu,gv)
    end select

    call interpolate_setuv(gu, gv)

    call find_points(gu, gv, 0.5d0*dt, midlon, midlat, deplon, deplat)   ! dtに0.5をかけているのは引数のdtが最初のステップ以外は2.0*deltatを渡しているから
    call set_niuv(dt)

    call legendre_synthesis(sphi_old, gphi_old)

    if (conserve) then
      gmin(:, :) = min(0.0d0, minval(gphi))
      gmax(:, :) = max(0.0d0, maxval(gphi))
    end if

    do j = 1, nlat
      do i = 1, nlon
        gphi(i,j) = gphi_old(p(i, j), q(i, j))
      end do
    end do

    ! dF/dlon
    call legendre_synthesis_dlon(sphi, dgphi)
    call bicubic_interpolation_set(dgphi)
    call interpolate_set(dgphi)
    call interpolate_setd(gphix, gphiy, gphixy)
    do j = 1, nlat
      do i = 1, nlon
        call interpolate_bicubic(midlon(i, j), midlat(i, j), dgphim(i, j))
      enddo
      gphim(:, j) = gum(:, j) * dgphim(:,j) / cos(latitudes(j)) ! gum: -u'
    enddo

    ! cos(lat)dF/dlat
    call legendre_synthesis_dlat(sphi, dgphi)
    call bicubic_interpolation_set(dgphi)
    call interpolate_set(dgphi)
    call interpolate_setd(gphix, gphiy, gphixy)
    do j = 1, nlat
      do i = 1, nlon
        call interpolate_bicubic(midlon(i, j), midlat(i, j), dgphim(i, j))
      enddo
      gphim(:, j) = gphim(:, j) + gvm(:, j) * dgphim(:,j) / cos(latitudes(j)) ! gvm: -v'
    enddo

    open(12, file="res_wind.txt")
    do i = 1, nlon
      do j = 1, nlat
        write(12,*) gum(i, j), gvm(i, j)
      end do
    end do
    close(12)

    gphi = gphi + dt * gphim

    if(conserve) then
      call mass_correct(gphi, gphi_old, gmax, gmin, w)
    endif

! time filter
    call legendre_analysis(gphi, sphi1)
    do m = 0, ntrunc
      sphi_old(m : ntrunc, m) = sphi(m : ntrunc, m)       
      sphi(m : ntrunc, m) = sphi1(m : ntrunc, m)
    enddo

  end subroutine update

  subroutine set_niuv(dt)
    use math_module, only: math_pi, pi2=>math_pi2
    use sphere_module, only: xyz2uv, lonlat2xyz
    use grid_module, only: pole_regrid
    implicit none

    real(8), intent(in) :: dt

    integer(8) :: i, j
    real(8) :: dlonr

    dlonr = 0.5d0 * nlon / math_pi
    gum(:, :) = 0.0d0
    gvm(:, :) = 0.0d0
    do j = 1, nlat
      do i = 1, nlon
        ! find grid points near departure points

        p(i, j) = int(anint( deplon(i, j) * dlonr + 1.0d0 ))
        if ( p(i,j) > nlon ) then
          p(i, j) = p(i, j) - nlon
        end if
        ! lat = (J+1-2j)pi/(2J+1)
        q(i, j) = int(anint( 0.5d0 * (nlat + 1.0d0 - (2.0d0*dble(nlat)+1.0d0)*deplat(i, j) / math_pi) ))  !latitudesは大きい順で詰められているので注意
        call pole_regrid(p(i, j), q(i, j))
        call calc_niuv(dt, p(i, j), q(i, j), longitudes(i), latitudes(j), midlon(i, j), midlat(i, j), gum(i, j), gvm(i, j))
      end do
    end do
        
  end subroutine  set_niuv

  subroutine bicubic_interpolation_set(f)
    use legendre_transform_module, only: legendre_analysis, legendre_synthesis_dlat, legendre_synthesis_dlon, &
      legendre_synthesis_dlonlat
    implicit none
    integer(8) :: j
    real(8), intent(in) :: f(nlon, nlat)

    call legendre_analysis(f, sphi1)
    call legendre_synthesis_dlon(sphi1, gphix)
    call legendre_synthesis_dlat(sphi1, gphiy)
    call legendre_synthesis_dlonlat(sphi1, gphixy)
    do j = 1, nlat
      gphiy(:,j) = gphiy(:,j) / cos(latitudes(j))
      gphixy(:,j) = gphixy(:,j) / cos(latitudes(j))
    end do

  end subroutine bicubic_interpolation_set

  ! Ritchie1987 式(45)のu^* - Uをgumに詰める(gvmも)
  subroutine calc_niuv(dt, p1, q1, lon, lat, midlon1, midlat1, gum1, gvm1)
    use grid_module, only: latitudes => lat, longitudes => lon
    use math_module, only: math_pi, pi2=>math_pi2
    use sphere_module, only: xyz2uv, lonlat2xyz
    use interpolate_module, only: interpolate_bilinearuv
    implicit none
    real(8), intent(in) :: dt
    integer(8), intent(in) :: p1, q1
    real(8), intent(in) :: lon, lat
    real(8), intent(out) :: gum1, gvm1
    real(8), intent(out) :: midlon1, midlat1

    real(8) :: xg, yg, zg, xr, yr, zr, xm, ym, zm, xdot, ydot, zdot, lon_grid, lat_grid, u, v, b


    lon_grid = longitudes(p1)
    lat_grid = latitudes(q1)
    call lonlat2xyz(lon_grid, lat_grid, xr, yr, zr)
    ! arrival points
    call lonlat2xyz(lon, lat, xg, yg, zg)

    b = 1.0d0 / sqrt( 2.0d0 * (1.0d0 + (xg*xr + yg*yr + zg*zr)) ) ! Ritchie1987 式(44)
    xm = b * (xg + xr)
    ym = b * (yg + yr)
    zm = b * (zg + zr)
    midlon1 = modulo(atan2(ym, xm) + pi2, pi2)
    midlat1 = asin(zm)

    xdot = (xg - xr) / dt
    ydot = (yg - yr) / dt
    zdot = (zg - zr) / dt
    call xyz2uv(xdot, ydot, zdot, midlon1, midlat1, u, v)  !Richie1987式(49)
    gum1 = gum1 + u
    gvm1 = gvm1 + v

    call interpolate_bilinearuv(midlon1, midlat1, u, v)
    gum1 = gum1 - u
    gvm1 = gvm1 - v

  end subroutine calc_niuv

end module nisl_module