module euler_module
  use planet_module, only: a=>planet_radius
  use grid_module, only: &
    nlon, nlat, ntrunc, &
    gu, gv, gphi, sphi, sphi_old, lon, lat, coslatr
  use time_module, only: nstep, hstep, deltat, velocity
  use legendre_transform_module, only: &
    legendre_analysis, legendre_synthesis, &
    legendre_synthesis_dlon, legendre_synthesis_dlat
  implicit none
  private

  real(8), dimension(:,:), allocatable, private :: dgphi

  complex(8), dimension(:,:), allocatable, private :: sphi1

  private :: update
  public :: eulerian_init, eulerian_timeint, eulerian_clean

contains

  subroutine eulerian_init()
    implicit none

    allocate(dgphi(nlon,nlat), sphi1(0:ntrunc,0:ntrunc))

    call legendre_synthesis(sphi,gphi)
    call update(0.0d0*deltat,0.5d0*deltat)
    call update(0.5d0*deltat,deltat)

  end subroutine eulerian_init

  subroutine eulerian_clean()

    deallocate(dgphi, sphi1)

  end subroutine eulerian_clean

  subroutine eulerian_timeint()
    implicit none

    integer(8) :: i

    do i=2, nstep
      call update((i-1)*deltat,2.0d0*deltat)
      write(*, *) "step=", i, "maxval = ", maxval(gphi)
    end do

  end subroutine eulerian_timeint

  subroutine update(t,dt)
    use uv_module, only: uv_nodiv, uv_div
    implicit none

    real(8), intent(in) :: t, dt

    integer(8) :: j, m, n

    select case(velocity)
      case("nodiv ")
        call uv_nodiv(t,lon,lat,gu,gv)
      case("div   ")
        call uv_div(t,lon,lat,gu,gv)
    end select
    call legendre_synthesis(sphi_old, gphi)
! dF/dlon
    call legendre_synthesis_dlon(sphi, dgphi)
    do j=1, nlat
      dgphi(:,j) =  dgphi(:,j)*coslatr(j)
    end do
    gphi(:, :) = gphi(:, :) - dt*gu(:, :)*dgphi(:, :)
! cos(lat)dF/dlat
    call legendre_synthesis_dlat(sphi, dgphi)
    do j=1, nlat
      dgphi(:,j) =  dgphi(:,j)*coslatr(j)
    end do
    gphi(:, :) = gphi(:, :) - dt*gv(:, :)*dgphi(:, :)
    call legendre_analysis(gphi, sphi1)
! update
    do m=0, ntrunc
      sphi_old(m:ntrunc,m) = sphi(m:ntrunc,m)
    end do
    do m=0, ntrunc
      do n=m, ntrunc
        sphi(n,m) = sphi1(n,m)
      end do
    end do

  end subroutine update

end module euler_module
