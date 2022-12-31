module interpolate2d_module
  ! interpolate in a stencil
    use grid_module, only: latitudes=>lat, wgt, longitudes=>lon, nlon, nlat, dh
    use sphere_module, only: lon2i, lat2j
    implicit none
    private

    integer(8), private :: nx, ny, nz, n=3, nh, nhalo, nx1, nx2, ny1, ny2
    integer(8), dimension(4), private :: is, js
    real(8), private :: u, t, dlon                            ! tとuは線分比
    real(8), dimension(:), allocatable, public :: lon_extend, lat_extend
    real(8), dimension(:,:,:), allocatable, private :: ff, ffx, ffy, ffxy

    public :: interpolate2d_init, interpolate2d_clean, &
              interpolate2d_set, interpolate2d_setd, interpolate2d_bicubic

  contains

    subroutine interpolate2d_init(f)
      use math_module, only: pi2=>math_pi2, pih=>math_pih
      implicit none

      real(8), dimension(:,:,:) :: f

      integer(8) :: i, j

      nh = n/2
      nhalo = nh

      nx = size(f,1)
      ny = size(f,2)
      nz = size(f, 3)
      nx1 = 1 - nhalo
      nx2 = nx + nhalo + 1
      ny1 = 1 - nhalo - 1
      ny2 = ny + nhalo + 1

      allocate(lon_extend(nx1:nx2), lat_extend(ny1:ny2), ff(nx1:nx2,ny1:ny2, nz), &
               ffx(nx1:nx2,ny1:ny2, nz), ffy(nx1:nx2,ny1:ny2, nz), ffxy(nx1:nx2,ny1:ny2, nz))

      dlon = pi2/nx
      do i=nx1, nx2
        lon_extend(i) = dlon*(i-1)
      end do
      lat_extend(1:ny) = latitudes
      do j=1, nhalo+1
        lat_extend(1-j)   = pih + (pih - latitudes(j))
        lat_extend(ny+j)   = -pih + (-pih - latitudes(ny-j+1))
      end do

    end subroutine interpolate2d_init

    subroutine interpolate2d_clean()
      implicit none

      deallocate(lon_extend, lat_extend, ff, ffx, ffy, ffxy)

    end subroutine  interpolate2d_clean

    subroutine interpolate2d_bicubic(lon, lat, h, fi)
      use bicubic_module, only: bcucof, bcuint, bcuintp
      implicit none

      real(8), intent(in) :: lon, lat, h
      real(8), intent(out) :: fi

      real(8), dimension(4) :: z, zx, zy, zxy
      integer(8) :: k, hid
      real(8) :: dlat

      hid = int(anint(h / dh)) + 1
      call find_stencil(lon, lat)
      dlat = lat_extend(js(4)) - lat_extend(js(1))
      do k=1, 4
        z(k) = ff(is(k),js(k), hid)
        zx(k) = ffx(is(k),js(k), hid)
        zy(k) = ffy(is(k),js(k), hid)
        zxy(k) = ffxy(is(k),js(k), hid)
      end do

      call bcucof(z,zx,zy,zxy,dlon,dlat)
      fi = bcuint(t,u)

    end subroutine interpolate2d_bicubic

    subroutine interpolate2d_set(f)
      implicit none

      real(8), dimension(:,:,:), intent(in) :: f

      integer(8) :: i, j

      ff(1:nx, 1:ny, :) = f
      do j=1, nh+1
        ff(1:nx, 1-j, :) = cshift(ff(1:nx, j, :),nx/2)
        ff(1:nx,ny+j, :) = cshift(ff(1:nx, ny-(j-1), :),nx/2)
      end do
      do i=1, nh
        ff(1-i,:,:) = ff(nx-(i-1),:,:)
      end do
      do i=1, nh+1
        ff(nx+i,:,:) = ff(1+(i-1),:,:)
      end do

    end subroutine interpolate2d_set

    subroutine interpolate2d_setd(fx,fy,fxy)
      implicit none

      real(8), dimension(:,:,:), intent(in) :: fx, fy, fxy

      integer(8) :: i, j

      ffx(1:nx,1:ny,:) = fx
      ffy(1:nx,1:ny,:) = fy
      ffxy(1:nx,1:ny,:) = fxy
  ! directions of d/dx and d/dy are reversed beyond poles
      do j=1, nh+1
        ffx(1:nx,1-j,:) = -cshift(ffx(1:nx,j,:),nx/2)
        ffx(1:nx,ny+j,:) = -cshift(ffx(1:nx,ny-(j-1),:),nx/2)
        ffy(1:nx,1-j,:) = -cshift(ffy(1:nx,j,:),nx/2)
        ffy(1:nx,ny+j,:) = -cshift(ffy(1:nx,ny-(j-1),:),nx/2)
        ffxy(1:nx,1-j,:) = cshift(ffxy(1:nx,j,:),nx/2)
        ffxy(1:nx,ny+j,:) = cshift(ffxy(1:nx,ny-(j-1),:),nx/2)
      end do
      do i=1, nh
        ffx(1-i,:,:) = ffx(nx-(i-1),:,:)
        ffy(1-i,:,:) = ffy(nx-(i-1),:,:)
        ffxy(1-i,:,:) = ffxy(nx-(i-1),:,:)
      end do
      do i=1, nh+1
        ffx(nx+i,:,:) = ffx(1+(i-1),:,:)
        ffy(nx+i,:,:) = ffy(1+(i-1),:,:)
        ffxy(nx+i,:,:) = ffxy(1+(i-1),:,:)
      end do

    end subroutine interpolate2d_setd

    subroutine find_stencil(lon, lat)
      implicit none

      real(8), intent(in) :: lon, lat

      integer(8) :: j

      is(1) = lon2i(lon,nx)
      is(2) = is(1) + 1
      t = lon/dlon - is(1) + 1.0d0 ! t = (lon - dlon*(i-1))/dlon
      is(3:4) = is(2:1:-1)

      j = lat2j(lat, ny) 
      if (lat > lat_extend(j)) then 
        j = j - 1
      end if
      js(1:2) = j
      js(3:4) = j + 1
      u = (lat-lat_extend(j))/(lat_extend(j+1)-lat_extend(j))

    end subroutine find_stencil
  end module interpolate2d_module
