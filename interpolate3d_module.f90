! 3次元の鉛直軸はzでやりたい。
module interpolate3d_module
  ! interpolate in a stencil
    use grid_module, only: latitudes=>lat, wgt, longitudes=>lon, nlon, nlat, dh
    use sphere_module, only: lon2i, lat2j
    implicit none
    private
  
    integer(8), private :: nx, ny, nz, n=3, nh, nhalo, nx1, nx2, ny1, ny2
    integer(8), dimension(8), private :: is, js, ks
    real(8), private :: u, t, v, dlon                            ! tとuは線分比
    real(8), dimension(:), allocatable, public :: lon_extend, lat_extend, height
    real(8), dimension(:, :, :), allocatable, private :: ff, ffx, ffy, ffxy, fu, fv
    real(8), dimension(:, :, :), allocatable, private :: ffz, ffxz, ffyz, ffxyz
  
    public :: interpolate_init, interpolate_clean, interpolate_set, interpolate_setd, interpolate_tricubic
  
  contains
  
    subroutine interpolate_init(f)
      use math_module, only: pi2=>math_pi2, pih=>math_pih
      use tricubic_module, only: tricubic_init
      implicit none
  
      real(8), dimension(:, :, :) :: f
  
      integer(8) :: i, j
  
      nh = n/2
      nhalo = nh
    
      nx = size(f, 1)
      ny = size(f, 2)
      nz = size(f, 3)
      nx1 = 1 - nhalo
      nx2 = nx + nhalo + 1
      ny1 = 1 - nhalo - 1
      ny2 = ny + nhalo + 1
  
      allocate(lon_extend(nx1:nx2), lat_extend(ny1:ny2), height(nz), ff(nx1:nx2,ny1:ny2,1:nz), &
               ffx(nx1:nx2,ny1:ny2,1:nz), ffy(nx1:nx2,ny1:ny2,1:nz), ffxy(nx1:nx2,ny1:ny2,1:nz), &
               fu(nx1:nx2,ny1:ny2,1:nz),fv(nx1:nx2,ny1:ny2,1:nz), ffxz(nx1:nx2,ny1:ny2,1:nz), &
               ffz(nx1:nx2,ny1:ny2,1:nz), ffyz(nx1:nx2,ny1:ny2,1:nz), ffxyz(nx1:nx2,ny1:ny2,1:nz))
  
      dlon = pi2/nx
      do i = nx1, nx2
        lon_extend(i) = dlon*(i-1)
      end do
      lat_extend(1:ny) = latitudes
      do j = 1, nhalo+1
        lat_extend(1-j)   = pih + (pih - latitudes(j))
        lat_extend(ny+j)   = -pih + (-pih - latitudes(ny-j+1))
      end do

      do i = 1, nz
        height(i) = dble(i - 1) * dh
      end do

      call tricubic_init()
  
    end subroutine interpolate_init
  
    subroutine interpolate_clean()
      implicit none
  
      deallocate(lon_extend, lat_extend, ff, ffx, ffy, ffxy, fu, fv)
  
    end subroutine  interpolate_clean
  
    subroutine interpolate_tricubic(lon, lat, h, fi)
      use tricubic_module, only: tricof, trint
      use planet_module, only: transorm_pressure_to_height
      implicit none
  
      real(8), intent(in) :: lon, lat, h
      real(8), intent(out) :: fi
  
      real(8), dimension(8) :: z, zx, zy, zz, zxy, zxz, zyz, zxyz
      integer(8) :: k
      real(8) :: dlat
      real(8), parameter :: eps = 1.0d-7

      call find_stencil(lon, lat, h)
      dlat = lat_extend(js(4)) - lat_extend(js(1))
      do k = 1, 8
        z(k) = ff(is(k), js(k), ks(k))
        zx(k) = ffx(is(k), js(k), ks(k))
        zy(k) = ffy(is(k), js(k), ks(k))
        zz(k) = ffz(is(k), js(k), ks(k))
        zxy(k) = ffxy(is(k), js(k), ks(k))
        zxz(k) = ffxz(is(k), js(k), ks(k))
        zyz(k) = ffyz(is(k), js(k), ks(k))
        zxyz(k) = ffxyz(is(k), js(k), ks(k))
      end do
  
      call tricof(z, zx, zy, zz, zxy, zxz, zyz, zxyz, dlon, dlat, dh)
      fi = trint(t, u, v)
     
    end subroutine interpolate_tricubic
  
    subroutine interpolate_set(f)
      implicit none
  
      real(8), dimension(:, :, :), intent(in) :: f
  
      integer(8) :: i, j
  
      ff(1:nx, 1:ny, 1:nz) = f
      do j=1, nh+1
        ff(1:nx, 1-j, 1:nz) = cshift(ff(1:nx, j, 1:nz),nx/2)
        ff(1:nx, ny+j, 1:nz) = cshift(ff(1:nx,ny-(j-1), 1:nz),nx/2)
      end do
      do i=1, nh
        ff(1-i, :, :) = ff(nx-(i-1), :, :)
      end do
      do i=1, nh+1
        ff(nx+i, :, :) = ff(1+(i-1), :, :)
      end do
  
    end subroutine interpolate_set
  
    subroutine interpolate_setd(fx, fy, fz, fxy, fxz, fyz, fxyz)
      implicit none
  
      real(8), dimension(:, :, :), intent(in) :: fx, fy, fz, fxy, fxz, fyz, fxyz
  
      integer(8) :: i, j
  
      ffx(1:nx, 1:ny, 1:nz) = fx
      ffy(1:nx, 1:ny, 1:nz) = fy
      ffz(1:nx, 1:ny, 1:nz) = fz
      ffxy(1:nx, 1:ny, 1:nz) = fxy
      ffxz(1:nx, 1:ny, 1:nz) = fxz
      ffyz(1:nx, 1:ny, 1:nz) = fyz
      ffxyz(1:nx, 1:ny, 1:nz) = fxyz
  ! directions of d/dx and d/dy are reversed beyond poles
      do j=1, nh+1
        ffx(1:nx, 1-j, 1:nz) = -cshift(ffx(1:nx, j, 1:nz),nx/2)
        ffx(1:nx,ny+j, 1:nz) = -cshift(ffx(1:nx, ny-(j-1), 1:nz),nx/2)

        ffy(1:nx, 1-j, 1:nz) = -cshift(ffy(1:nx,j, 1:nz),nx/2)
        ffy(1:nx, ny+j, 1:nz) = -cshift(ffy(1:nx,ny-(j-1), 1:nz),nx/2)

        ffz(1:nx, 1-j, 1:nz) = -cshift(ffz(1:nx, j, 1:nz),nx/2)
        ffz(1:nx,ny+j, 1:nz) = -cshift(ffz(1:nx, ny-(j-1), 1:nz),nx/2)

        ffxy(1:nx, 1-j, 1:nz) = cshift(ffxy(1:nx, j, 1:nz),nx/2)
        ffxy(1:nx,ny+j, 1:nz) = cshift(ffxy(1:nx,ny-(j-1), 1:nz),nx/2)

        ffxz(1:nx, 1-j, 1:nz) = cshift(ffxz(1:nx, j, 1:nz),nx/2)
        ffxz(1:nx,ny+j, 1:nz) = cshift(ffxz(1:nx,ny-(j-1), 1:nz),nx/2)

        ffyz(1:nx, 1-j, 1:nz) = cshift(ffyz(1:nx, j, 1:nz),nx/2)
        ffyz(1:nx,ny+j, 1:nz) = cshift(ffyz(1:nx,ny-(j-1), 1:nz),nx/2)

        ffxyz(1:nx, 1-j, 1:nz) = cshift(ffxyz(1:nx, j, 1:nz),nx/2)
        ffxyz(1:nx,ny+j, 1:nz) = cshift(ffxyz(1:nx,ny-(j-1), 1:nz),nx/2)
      end do
      do i=1, nh
        ffx(1-i, :, :) = ffx(nx-(i-1), :, :)
        ffy(1-i, :, :) = ffy(nx-(i-1), :, :)
        ffxy(1-i, :, :) = ffxy(nx-(i-1), :, :)
        ffz(1-i, :, :) = ffz(nx-(i-1), :, :)
        ffxz(1-i, :, :) = ffxz(nx-(i-1), :, :)
        ffyz(1-i, :, :) = ffyz(nx-(i-1), :, :)
        ffxyz(1-i,:,:) = ffxyz(nx-(i-1), :, :)
      end do
      do i=1, nh+1
        ffx(nx+i, :, :) = ffx(1+(i-1), :, :)
        ffy(nx+i, :, :) = ffy(1+(i-1), :, :)
        ffxy(nx+i, :, :) = ffxy(1+(i-1), :, :)
        ffz(nx+i, :, :) = ffz(1+(i-1), :, :)
        ffxz(nx+i, :, :) = ffxz(1+(i-1), :, :)
        ffyz(nx+i, :, :) = ffyz(1+(i-1), :, :)
        ffxyz(nx+i, :, :) = ffxyz(1+(i-1), :, :)
      end do
  
    end subroutine interpolate_setd
  
    subroutine find_stencil(lon, lat, h)
      implicit none
  
      real(8), intent(in) :: lon, lat, h
   
      integer(8) :: j
  
      is(1) = lon2i(lon, nx); is(2) = is(1) + 1;
      is(3) = is(1); is(5) = is(1); is(7) = is(1)
      is(4) = is(2); is(6) = is(2); is(8) = is(2)
      t = lon/dlon - is(1) + 1.0d0

      j = lat2j(lat, ny)
      if (lat > lat_extend(j)) then 
        j = j - 1
      end if
      js(1:2) = j; js(5:6) = j
      js(3:4) = j + 1; js(7:8) = j + 1
      u = (lat-lat_extend(j))/(lat_extend(j+1)-lat_extend(j))

      ks(1) = int(h/dh) + 1; ks(5) = ks(1) + 1
      ks(2) = ks(1); ks(3) = ks(1); ks(4) = ks(1)
      ks(6) = ks(5); ks(7) = ks(5); ks(8) = ks(5)

      v = (h - (ks(1)-1)*dh) / dh

      if(t > 1.1d0 .or. u > 1.1d0 .or. v > 1.1d0) then
        write(*,*) 'tuv', t, u, v
      endif

      if(t < -0.01d0 .or. u < -0.01d0 .or. v < -0.01d0) then
        write(*,*) 'tuv', t, u, v
      endif
  
    end subroutine find_stencil
end module interpolate3d_module
