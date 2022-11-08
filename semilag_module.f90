module semilag_module

  use grid_module, only: nlon, nlat, ntrunc, gphi_, &
    gu, gv, gphi, sphi_old, sphi, latitudes=>lat, lon, coslatr, wgt
  use field_module, only : X, Y
  use time_module, only: velocity
  private
  
  real(8), dimension(:,:), allocatable, private :: &
    gphi_old, gphix, gphiy, gphixy, midlon, midlat, deplon, deplat, gphi_initial

  private :: update
  public :: semilag_init, semilag_timeint, semilag_clean

contains

  subroutine semilag_init()
    use interpolate_module, only: interpolate_init
    use legendre_transform_module, only: legendre_synthesis
    implicit none

    integer(8) :: i,j

    allocate(gphi_old(nlon,nlat), gphix(nlon,nlat),gphiy(nlon,nlat),gphixy(nlon,nlat), &
             midlon(nlon,nlat),midlat(nlon,nlat),deplon(nlon,nlat),deplat(nlon,nlat), gphi_initial(nlon, nlat))
    call interpolate_init(gphi)

    gphi_old(:,:) = gphi(:,:)
    gphi_initial(:, :) = gphi(:, :)

    do i=1, nlon
      midlon(i,:) = lon(i)
    end do
    do j=1, nlat
      midlat(:,j) = latitudes(j)
    end do
    open(11, file="animation.txt")
    do i = 1, nlon
      do j = 1, nlat
          write(11,*) lon(i), latitudes(j), gphi(i, j)
      end do        
    end do

  end subroutine semilag_init

  subroutine semilag_clean()
    use interpolate_module, only: interpolate_clean
    implicit none

    deallocate(gphi_old,gphix,gphiy,gphixy,midlon,midlat,deplon,deplat)
    call interpolate_clean()

  end subroutine semilag_clean

  subroutine semilag_timeint()
    use time_module, only: nstep, deltat, hstep, field
    use legendre_transform_module, only: legendre_synthesis
    implicit none

    integer(8) :: i, j, k

    do i=1, nstep
      call update((i-0.5d0)*deltat, deltat)
      write(*, *) "step=", i, "maxval = ", maxval(gphi), 'minval = ', minval(gphi)
      if ( mod(i, hstep) == 0 ) then
        do j = 1, nlon
            do k = 1, nlat
                write(11,*) lon(j), latitudes(k), gphi(j, k)
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
    open(10, file="log.txt")
    do i = 1, nlon
      do j = 1, nlat
        write(10,*) lon(i), latitudes(j), gphi(i, j)
      enddo
    enddo
    close(10)
    open(12, file="error.txt")
    do i = 1, nlon
        do j = 1, nlat
            write(12,*) lon(i), latitudes(j), gphi_initial(i, j) - gphi(i, j)
        end do
    end do

  end subroutine semilag_timeint

  subroutine update(t, dt)
    use math_module, only: &
      pi=>math_pi, pir=>math_pir, pih=>math_pih
    use upstream_module, only: find_points
    use time_module, only: imethod
    use uv_module, only: uv_sbody, uv_nodiv, uv_div
    use interpolate_module, only: &
      interpolate_set, interpolate_setd, interpolate_setdx, &
      interpolate_bilinear, interpolate_bicubic, interpolate_polin2, &
      interpolate_linpol, interpolate_setdx, interpolate_diff
    use legendre_transform_module, only: legendre_analysis, legendre_synthesis, &
        legendre_synthesis_dlon, legendre_synthesis_dlat, legendre_synthesis_dlonlat
    implicit none

    real(8), intent(in) :: t, dt

    integer(8) :: i, j, m, n
    real(8) :: eps, dlonr
    real(8), dimension(nlon) :: gphitmp

    select case(velocity)
    case("nodiv")
      call uv_nodiv(t,lon,latitudes,gu,gv)
    case("div")
      call uv_div(t,lon,latitudes,gu,gv)
    case("sbody")
      call uv_sbody(lon, latitudes, gu, gv)
    case default
      print *, "No matching model for", velocity
    end select
    call find_points(gu, gv, 0.5d0*dt, midlon, midlat, deplon, deplat)

    call legendre_synthesis(sphi_old,gphi_old)

    if ((imethod=="sph").or.(imethod=="fdy")) then
      call legendre_analysis(gphi_old,sphi_old)
    end if

! calculate spectral derivatives

    if ((imethod=="sph").or.(imethod=="fdy")) then
      call legendre_synthesis_dlon(sphi_old,gphix)
    end if
    if (imethod=="sph") then
      call legendre_synthesis_dlat(sphi_old,gphiy)
      call legendre_synthesis_dlonlat(sphi_old,gphixy)
      do j=1, nlat
        gphiy(:,j) = gphiy(:,j)*coslatr(j)
        gphixy(:,j) = gphixy(:,j)*coslatr(j)
      end do
    end if

! calculate fd derivatives

    if (imethod=="fd") then
! d/dlon
      dlonr = 0.25d0*nlon*pir
      gphix(1,:) = dlonr * (gphi_old(2,:) - gphi_old(nlon,:))
      gphix(nlon,:) = dlonr * (gphi_old(1,:) - gphi_old(nlon-1,:))
      do i=2, nlon-1
        gphix(i,:) = dlonr*(gphi_old(i+1,:) - gphi_old(i-1,:))
      end do
    end if
    if ((imethod=="fd").or.(imethod=="fdy")) then
! d/dphi
      eps = pih-latitudes(1)
      gphitmp = cshift(gphi_old(:,1),nlon/2)
      gphiy(:,1) = (gphitmp-gphi_old(:,2))/(pih+eps-latitudes(2))
      gphitmp = cshift(gphix(:,1),nlon/2)
      gphixy(:,1) = (gphitmp-gphix(:,2))/(pih+eps-latitudes(2))
      gphitmp = cshift(gphi_old(:,nlat),nlon/2)
      gphiy(:,nlat) = (gphitmp-gphi_old(:,nlat-1))/(-pih-eps-latitudes(nlat-1))
      gphitmp = cshift(gphix(:,nlat),nlon/2)
      gphixy(:,nlat) = (gphitmp-gphix(:,nlat-1))/(-pih-eps-latitudes(nlat-1))
      do j=2, nlat-1
        gphiy(:,j) = (gphi_old(:,j+1)-gphi_old(:,j-1))/(latitudes(j+1)-latitudes(j-1))
        gphixy(:,j) = (gphix(:,j+1)-gphix(:,j-1))/(latitudes(j+1)-latitudes(j-1))
      end do 
    end if

! set grids
    call interpolate_set(gphi_old)
    if (imethod=="spcher") then
      call interpolate_setdx(gphix)
    end if
    if ((imethod=="fd    ").or.(imethod=="sph   ").or. &
        (imethod=="fdy   ")) then
      call interpolate_setd(gphix, gphiy, gphixy)
    end if
    do j=1, nlat
      do i=1, nlon
        select case (imethod)
          case ("bilin ") ! monotonicity is guranteed
            call interpolate_bilinear(deplon(i,j), deplat(i,j), gphi(i,j))
          case ("fd    ", "sph   ", "fdy   ")
            call interpolate_bicubic(deplon(i,j), deplat(i,j), gphi(i,j))
          case ("polin2")
            call interpolate_polin2(deplon(i,j), deplat(i,j), gphi(i,j))
          case ("linpol")
            call interpolate_linpol(deplon(i,j), deplat(i,j), gphi(i,j))
        end select
      end do
    end do


! spectral
    call legendre_analysis(gphi,sphi)
    do n=1, ntrunc
      do m=0, n
        sphi_old(n,m) = sphi(n,m)
      end do
    end do

  end subroutine update

end module semilag_module