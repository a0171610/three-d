module semilag_module

  use grid_module, only: nlon, nlat, nz, ntrunc, gphi, pres, rho, height, &
    gu, gv, gomega, gphi, sphi_old, sphi, latitudes=>lat, lon, coslatr, wgt
  private
  
  real(8), dimension(:, :, :), allocatable, private :: midlon, midlat, deplon, deplat, deppres
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
             deppres(nlon, nlat, nz))
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
    do i = 1, nlon
      do j = 1, nlat
        write(11,*) lon(i), latitudes(j), gphi(i, j, 24)
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

    do i = 1, nstep/2
      call update((i-0.5d0)*deltat, deltat)
      write(*, *) "step=", i, "maxval = ", maxval(gphi), 'minval = ', minval(gphi)
      if ( mod(i, hstep) == 0 ) then
        do j = 1, nlon
            do k = 1, nlat
                write(11,*) lon(j), latitudes(k), gphi(j, k, 25)
            end do
        end do
      endif
      if (i == nstep / 2) then
        open(10, file="half.txt")
        do j = 1, nlon
          do k = 1, nlat
            write(10,*) lon(j), latitudes(k), gphi(j, k, 25)
          enddo
        enddo
        close(10)
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
    open(12, file="error.txt")
    do i = 1, nlon
        do j = 1, nlat
      !      write(12,*) lon(i), latitudes(j), gphi_initial(i, j) - gphi(i, j)
        end do
    end do

  end subroutine semilag_timeint

  subroutine update(t, dt)
    use math_module, only: &
      pi=>math_pi, pir=>math_pir, pih=>math_pih
    use upstream3d_module, only: find_points
    use time_module, only: case
    use uv_module, only: uv_div
    use uv_hadley_module, only: uv_hadley
    use interpolate3d_module, only: interpolate_set, interpolate_setd, interpolate_tricubic
    use legendre_transform_module, only: legendre_analysis, legendre_synthesis, &
        legendre_synthesis_dlon, legendre_synthesis_dlat, legendre_synthesis_dlonlat
    implicit none

    real(8), intent(in) :: t, dt

    integer(8) :: i, j, k, m, n

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

    do k = 2, nz-1
      gphiz(:, :, k) = (gphi(:, :, k + 1) - gphi(:, :, k - 1)) / (height(k+1) - height(k-1))
      gphixz(:, :, k) = (gphix(:, :, k + 1) - gphix(:, :, k - 1)) / (height(k+1) - height(k-1))
      gphiyz(:, :, k) = (gphiy(:, :, k + 1) - gphiy(:, :, k - 1)) / (height(k+1) - height(k-1))
      gphixyz(:, :, k) = (gphixy(:, :, k + 1) - gphixy(:, :, k - 1)) / (height(k+1) - height(k-1))
    end do

    gphiz(:, :, 1) = (gphi(:, :, 2) - gphi(:, :, 1)) / (height(2) - height(1))
    gphixz(:, :, 1) = (gphix(:, :, 2) - gphix(:, :, 1)) / (height(2) - height(1))
    gphiyz(:, :, 1) = (gphiy(:, :, 2) - gphiy(:, :, 1)) / (height(2) - height(1))
    gphixyz(:, :, 1) = (gphixy(:, :, 2) - gphixy(:, :, 1)) / (height(2) - height(1))

    gphiz(:, :, nz) = (gphi(:, :, nz) - gphi(:, :, nz - 1)) / (height(nz) - height(nz-1))
    gphixz(:, :, nz) = (gphix(:, :, nz) - gphix(:, :, nz - 1)) / (height(nz) - height(nz-1))
    gphiyz(:, :, nz) = (gphiy(:, :, nz) - gphiy(:, :, nz - 1)) / (height(nz) - height(nz-1))
    gphixyz(:, :, nz) = (gphixy(:, :, nz) - gphixy(:, :, nz - 1)) / (height(nz) - height(nz-1))

! set grids
    call interpolate_set(gphi_old)
    call interpolate_setd(gphix, gphiy, gphiz, gphixy, gphixz, gphiyz, gphixyz)
    do j = 1, nlat
      do i = 1, nlon
        do k = 1, nz
          call interpolate_tricubic(deplon(i,j,k), deplat(i,j,k), deppres(i, j, k), gphi(i,j,k))
        enddo
      enddo
    end do

! spectral
    do k = 1, nz
      call legendre_analysis(gphi(:, :, k), sphi(:, :, k))      
    end do
    do n = 1, ntrunc
      do m = 0, n
        do k = 1, nz
          sphi_old(n, m, k) = sphi(n, m, k)
        end do
      enddo
    end do

  end subroutine update

end module semilag_module