module interpolate_module
  ! interpolate in a stencil
    use grid_module, only: latitudes=>lat, wgt, longitudes=>lon, nlon, nlat
    use sphere_module, only: lon2i, lat2j
    implicit none
    private
  
    integer(8), private :: nx, ny, n=3, nh, nhalo, nx1, nx2, ny1, ny2 
    integer(8), dimension(4), private :: is, js
    real(8), private :: u, t, dlon                            ! tとuは線分比
    real(8), dimension(:), allocatable, public :: lon_extend, lat_extend
    real(8), dimension(:,:), allocatable, private :: ff, ffx, ffy, ffxy, fu, fv, ffxl, ffyl
    integer(8), dimension(:, :, :), allocatable, public :: record_i, record_j, record_d
    integer(8), dimension(:, :), allocatable, private :: record_count
  
    public :: interpolate_init, interpolate_clean, &
              interpolate_set, interpolate_setuv, &
              interpolate_setd, interpolate_setdx, &
              interpolate_bilinear, interpolate_bilinearuv, &
              interpolate_polin2, interpolate_polin2uv, &
              interpolate_bicubic, interpolate_linpol, &
              interpolate_diff, interpolate_bilinear_ratio, find_stencil_, &
              interpolate_dist, interpolate_dist_ratio
  
  contains
  
    subroutine interpolate_init(f)
      use math_module, only: pi2=>math_pi2, pih=>math_pih
      implicit none
  
      real(8), dimension(:,:) :: f
  
      integer(8) :: i, j
  
      nh = n/2
      nhalo = nh
    
      nx = size(f,1)
      ny = size(f,2)
      nx1 = 1 - nhalo
      nx2 = nx + nhalo + 1
      ny1 = 1 - nhalo - 1
      ny2 = ny + nhalo + 1

      allocate(lon_extend(nx1:nx2), lat_extend(ny1:ny2), ff(nx1:nx2,ny1:ny2), &
               ffx(nx1:nx2,ny1:ny2), ffy(nx1:nx2,ny1:ny2), ffxy(nx1:nx2,ny1:ny2), &
               ffxl(nx1:nx2,ny1:ny2), ffyl(nx1:nx2,ny1:ny2), &
               fu(nx1:nx2,ny1:ny2),fv(nx1:nx2,ny1:ny2), record_i(nlon, nlat, 60), &
               record_j(nlon, nlat, 60), record_d(nlon, nlat, 60), record_count(nlon, nlat))
  
      dlon = pi2/nx
      do i=nx1, nx2
        lon_extend(i) = dlon*(i-1)
      end do
      lat_extend(1:ny) = latitudes
      do j=1, nhalo+1
        lat_extend(1-j)   = pih + (pih - latitudes(j))
        lat_extend(ny+j)   = -pih + (-pih - latitudes(ny-j+1))
      end do
  
    end subroutine interpolate_init
  
    subroutine interpolate_clean()
      implicit none
  
      deallocate(lon_extend, lat_extend, ff, ffx, ffy, ffxy, fu, fv, ffxl, ffyl)
  
    end subroutine  interpolate_clean
  
    subroutine interpolate_bilinear(lon, lat, fi)
      implicit none
  
      real(8), intent(in) :: lon, lat
      real(8), intent(out) :: fi
  
      real(8), dimension(4) :: fs
      integer(8) :: k 
  
      call find_stencil(lon, lat)
      do k=1, 4
        fs(k) = ff(is(k),js(k))
      end do
  
  ! Due to geographical location weight t and 1-t should be swapped.
  ! However, not so sensitive to t or 1-t at points 3 and 4,
  ! even slightly worse.  Leave it simple.
  
      fi = (1.0d0-u)*((1.0d0-t)*fs(1)+t*fs(2)) + u*(t*fs(3)+(1.0d0-t)*fs(4))
  
    end subroutine interpolate_bilinear

    subroutine interpolate_bilinear_ratio(lon, lat, A, B, C, D)
      implicit none
      real(8), intent(in) :: lon, lat
      real(8), intent(out) :: A, B, C, D

      call find_stencil(lon, lat)
      A = (1.0d0-t) * (1.0d0-u)
      B = t * (1.0d0-u)
      C = t * u
      D = (1.0d0-t) * u        
    end subroutine interpolate_bilinear_ratio

    subroutine interpolate_dist_ratio(lon, lat, A, B, C, D)
      use sphere_module, only: orthodrome
      use grid_module, only: pole_regrid
      implicit none
      real(8), intent(in) :: lon, lat
      real(8), intent(out) :: A, B, C, D
      real(8) :: dist1, dist2, dist3, dist4

      call find_stencil(lon, lat)
      call pole_regrid(is(1), js(1))
      call pole_regrid(is(2), js(2))
      call pole_regrid(is(3), js(3))
      call pole_regrid(is(4), js(4))

      dist1 = 1.0d0/(orthodrome(lon, lat, longitudes(is(1)), latitudes(js(1)))**6)
      dist2 = 1.0d0/(orthodrome(lon, lat, longitudes(is(2)), latitudes(js(2)))**6)
      dist3 = 1.0d0/(orthodrome(lon, lat, longitudes(is(3)), latitudes(js(3)))**6)
      dist4 = 1.0d0/(orthodrome(lon, lat, longitudes(is(4)), latitudes(js(4)))**6)
      A = dist1 / (dist1 + dist2 + dist3 + dist4)
      B = dist2 / (dist1 + dist2 + dist3 + dist4)
      C = dist3 / (dist1 + dist2 + dist3 + dist4)
      D = dist4 / (dist1 + dist2 + dist3 + dist4)

    end subroutine interpolate_dist_ratio

    subroutine interpolate_dist(lon, lat, fi)
      implicit none

      real(8), intent(in) :: lon, lat
      real(8), intent(out) :: fi
      real(8) :: A, B, C, D
      real(8), dimension(4) :: fs
      integer(8) :: k

      call find_stencil(lon, lat)
      do k=1, 4
        fs(k) = ff(is(k),js(k))
      end do

      call interpolate_dist_ratio(lon, lat, A, B, C, D)

      fi = A * fs(1) + B * fs(2) + C * fs(3) + D * fs(4)

    end subroutine interpolate_dist
  
    subroutine interpolate_bilinearuv(lon, lat, fiu, fiv)
      implicit none
  
      real(8), intent(in) :: lon, lat
      real(8), intent(out) :: fiu, fiv
  
      real(8), dimension(4) :: fsu, fsv
      integer(8) :: k
  
      call find_stencil(lon, lat)
      do k=1, 4
        fsu(k) = fu(is(k),js(k))
        fsv(k) = fv(is(k),js(k))
      end do
      fiu = (1.0d0-u)*((1.0d0-t)*fsu(1)+t*fsu(2)) + u*(t*fsu(3)+(1.0d0-t)*fsu(4))
      fiv = (1.0d0-u)*((1.0d0-t)*fsv(1)+t*fsv(2)) + u*(t*fsv(3)+(1.0d0-t)*fsv(4))
  
    end subroutine interpolate_bilinearuv
  
    subroutine interpolate_bicubic(lon, lat, fi)
      use bicubic_module, only: bcucof, bcuint, bcuintp
      implicit none
  
      real(8), intent(in) :: lon, lat
      real(8), intent(out) :: fi
  
      real(8), dimension(4) :: z, zx, zy, zxy
      integer(8) :: k
      real(8) :: dlat
  
      call find_stencil(lon, lat)
      dlat = lat_extend(js(4)) - lat_extend(js(1))
      do k=1, 4
        z(k) = ff(is(k),js(k))
        zx(k) = ffx(is(k),js(k))
        zy(k) = ffy(is(k),js(k))
        zxy(k) = ffxy(is(k),js(k))
      end do
  
      call bcucof(z,zx,zy,zxy,dlon,dlat)
      fi = bcuint(t,u)
  
     
    end subroutine interpolate_bicubic
  
    subroutine interpolate_polin2(lon, lat, fi)
      use polint_module, only : polin2
      implicit none
  
      real(8), intent(in) :: lon, lat
      real(8), intent(out) :: fi
  
      integer(8) :: i0, i1, i2, j0, j1, j2
      real(8) :: dfi
  
      call find_stencil(lon, lat)
      i0 = is(1)
      i1 = i0 - nh
      i2 = i0 + nh + 1
      j0 = js(1)
      j1 = j0 - nh
      j2 = j0 + nh + 1
      call polin2(lon_extend(i1:i2), lat_extend(j1:j2), ff(i1:i2,j1:j2), lon, lat, fi, dfi)
  
  
    end subroutine interpolate_polin2
  
    subroutine interpolate_polin2uv(lon, lat, fiu, fiv, monotonic)
      use polint_module, only : polin2
      implicit none
  
      real(8), intent(in) :: lon, lat
      real(8), intent(out) :: fiu, fiv
      logical, optional, intent(in) :: monotonic
  
      integer(8) :: i0, i1, i2, j0, j1, j2
      real(8) :: dfi
  
      call find_stencil(lon, lat)
      i0 = is(1)
      i1 = i0 - nh
      i2 = i0 + nh + 1
      j0 = js(1)
      j1 = j0 - nh
      j2 = j0 + nh + 1
      call polin2(lon_extend(i1:i2), lat_extend(j1:j2), fu(i1:i2,j1:j2), lon, lat, fiu, dfi)
      call polin2(lon_extend(i1:i2), lat_extend(j1:j2), fv(i1:i2,j1:j2), lon, lat, fiv, dfi)
  
  ! Bermejo and Staniforth 1992
      if (present(monotonic).and.(monotonic)) then
        fiu = min(fiu,maxval(fu(i0:i0+1,j0:j0+1)))
        fiu = max(fiu,minval(fu(i0:i0+1,j0:j0+1)))
        fiv = min(fiv,maxval(fv(i0:i0+1,j0:j0+1)))
        fiv = max(fiv,minval(fv(i0:i0+1,j0:j0+1)))
      end if
  
    end subroutine interpolate_polin2uv
  
    subroutine interpolate_linpol(lon, lat, fi)
      use polint_module, only : polint
      implicit none
  
      real(8), intent(in) :: lon, lat
      real(8), intent(out) :: fi
  
      integer(8) :: i0, i1, i2, j0, j1, j2, j
      real(8) :: dfi
      real(8), dimension(n+1) :: ytmp
  
      call find_stencil(lon, lat)
      i0 = is(1)
      j0 = js(1)
      i1 = i0 - nh
      i2 = i0 + nh + 1
      j1 = j0 - nh
      j2 = j0 + nh + 1
      do j=j1, j0-1
        ytmp(j-j1+1) =  (1.0d0-t)*ff(i0,j)+t*ff(i0+1,j)
      end do
      do j=j0, j0+1
        call polint(lon_extend(i1:i2), ff(i1:i2,j), lon, ytmp(j-j1+1), dfi)
      end do
      do j=j0+2, j2
        ytmp(j-j1+1) =  (1.0d0-t)*ff(i0,j)+t*ff(i0+1,j)
      end do
      call polint(lat_extend(j1:j2), ytmp, lat, fi, dfi)
  
    end subroutine interpolate_linpol
  
    subroutine interpolate_set(f)
      implicit none
  
      real(8), dimension(:,:), intent(in) :: f
  
      integer(8) :: i, j
  
      ff(1:nx,1:ny) = f
      do j=1, nh+1
        ff(1:nx,1-j) = cshift(ff(1:nx,j),nx/2)
        ff(1:nx,ny+j) = cshift(ff(1:nx,ny-(j-1)),nx/2)
      end do
      do i=1, nh
        ff(1-i,:) = ff(nx-(i-1),:)
      end do
      do i=1, nh+1
        ff(nx+i,:) = ff(1+(i-1),:)
      end do
  
    end subroutine interpolate_set
  
    subroutine interpolate_setuv(gu,gv)
      implicit none
  
      real(8), dimension(:,:), intent(in) :: gu, gv
  
      integer(8) :: i, j
  
      do j=1, ny
        fu(1:nx,j) = gu(:,j)
        fv(1:nx,j) = gv(:,j)
      end do
  ! direction of u, v is reversed beyond poles
      do j=1, nh+1
        fu(1:nx,1-j) = -cshift(fu(1:nx,j),nx/2)
        fu(1:nx,ny+j) = -cshift(fu(1:nx,ny-(j-1)),nx/2)
        fv(1:nx,1-j) = -cshift(fv(1:nx,j),nx/2)
        fv(1:nx,ny+j) = -cshift(fv(1:nx,ny-(j-1)),nx/2)
      end do
      do i=1, nh
        fu(1-i,:) = fu(nx-(i-1),:)
        fv(1-i,:) = fv(nx-(i-1),:)
      end do
      do i=1, nh+1
        fu(nx+i,:) = fu(1+(i-1),:)
        fv(nx+i,:) = fv(1+(i-1),:)
      end do
  
    end subroutine interpolate_setuv
  
    subroutine interpolate_setd(fx,fy,fxy)
      implicit none
  
      real(8), dimension(:,:), intent(in) :: fx, fy, fxy
  
      integer(8) :: i, j
  
      ffx(1:nx,1:ny) = fx
      ffy(1:nx,1:ny) = fy
      ffxy(1:nx,1:ny) = fxy
  ! directions of d/dx and d/dy are reversed beyond poles
      do j=1, nh+1
        ffx(1:nx,1-j) = -cshift(ffx(1:nx,j),nx/2)
        ffx(1:nx,ny+j) = -cshift(ffx(1:nx,ny-(j-1)),nx/2)
        ffy(1:nx,1-j) = -cshift(ffy(1:nx,j),nx/2)
        ffy(1:nx,ny+j) = -cshift(ffy(1:nx,ny-(j-1)),nx/2)
        ffxy(1:nx,1-j) = cshift(ffxy(1:nx,j),nx/2)
        ffxy(1:nx,ny+j) = cshift(ffxy(1:nx,ny-(j-1)),nx/2)
      end do
      do i=1, nh
        ffx(1-i,:) = ffx(nx-(i-1),:)
        ffy(1-i,:) = ffy(nx-(i-1),:)
        ffxy(1-i,:) = ffxy(nx-(i-1),:)
      end do
      do i=1, nh+1
        ffx(nx+i,:) = ffx(1+(i-1),:)
        ffy(nx+i,:) = ffy(1+(i-1),:)
        ffxy(nx+i,:) = ffxy(1+(i-1),:)
      end do
  
    end subroutine interpolate_setd
  
    subroutine interpolate_setdx(fx)
      implicit none
  
      real(8), dimension(:,:), intent(in) :: fx
  
      integer(8) :: i, j
  
      ffx(1:nx,1:ny) = fx
  ! direction of d/dx is reversed beyond poles
      do j=1, nh+1
        ffx(1:nx,1-j) = -cshift(ffx(1:nx,j),nx/2)
        ffx(1:nx,ny+j) = -cshift(ffx(1:nx,ny-(j-1)),nx/2)
      end do
      do i=1, nh
        ffx(1-i,:) = fx(nx-(i-1),:)
      end do
      do i=1, nh+1
        ffx(nx+i,:) = fx(1+(i-1),:)
      end do
      ffy(:,:) = 0.d0
      ffxy(:,:) = 0.d0
  
    end subroutine interpolate_setdx
  
    subroutine interpolate_diff()
      implicit none
  
      integer(8) :: i, j
  
      forall(i=1:nx, j=1:ny)
        ffxl(i,j) = ff(i+1,j) - ff(i,j)
        ffyl(i,j) = ff(i,j+1) - ff(i,j)
      end forall
  
    end subroutine interpolate_diff
  
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

    ! find_stencil関数と違って、返り値を極に関して丸めている
    subroutine find_stencil_(lon, lat, is_, js_)
      use grid_module, only: pole_regrid
      implicit none
      real(8), intent(in) :: lon, lat
      integer(8), dimension(:), intent(out) :: is_, js_
   
      integer(8) :: j
  
      is_(1) = lon2i(lon, nx)
      is_(2) = is_(1) + 1
      t = lon/dlon - is_(1) + 1.0d0 ! t = (lon - dlon*(i-1))/dlon
      is_(3:4) = is_(2:1:-1)
  
      j = lat2j(lat, ny)
      if (lat > lat_extend(j)) then
        j = j - 1
      end if
      js_(1 : 2) = j
      js_(3 : 4) = j + 1
      do j = 1, 4
          call pole_regrid(is_(j), js_(j))
      end do
      u = (lat - lat_extend(j)) / (lat_extend(j+1) - lat_extend(j))

    end subroutine find_stencil_

  end module interpolate_module
