module cascade_interpolation_module
  use grid_module, only: latitudes=>lat, lon, nlat, nlon
  implicit none
  real(8), dimension(:, :), allocatable, private :: intersect_lon, intersect
  real(8), dimension(:), allocatable, private ::  intersect_lon_sort, insert_lon_circle, insert_circle, &
    circle_lat_dist, gphi_old_circle, intersect_sort, sphere_dist, sphere_loc

contains

  subroutine cascade_interpolate_init
    allocate(intersect_lon(nlon, nlat), intersect(nlon, nlat), intersect_lon_sort(nlon))
    allocate(insert_lon_circle(0:nlon+1), insert_circle(0:nlon+1), &
      circle_lat_dist(0 : nlat-1), gphi_old_circle(0 : nlat-1), intersect_sort(nlon), sphere_dist(nlat), sphere_loc(nlat))
  end subroutine cascade_interpolate_init
  

  subroutine cascade_interpolate(deplon, deplat, gphi_old, gphi)
    use sphere_module, only: lonlat2xyz, orthodrome
    use math_module, only: pi2=>math_pi2, pih=>math_pih, pi=>math_pi
    use search_module, only: search_bisection
    use sort_module, only: sort_two_array, sort_three_array
    use spline_interpolate_module, only: interpolate_spline
    implicit none

    real(8), intent(in) :: deplon(:, :), deplat(:, :), gphi_old(:, :)
    real(8), intent(out) :: gphi(:, :)

    integer(8) :: i, j, k
    real(8) :: target_lat

    do i = 1, nlon
      do j = 1, nlat
        target_lat = latitudes(j)
        k = search_bisection(deplat(i, :), target_lat)

        if (1 <= k .and. k < nlat) then
          call intersect_points(deplon(i, k), deplat(i, k), deplon(i, k+1), deplat(i,k+1), target_lat, intersect_lon(i, j))
        else
          if (k == 0) then
            call intersect_points(deplon(i, 1), deplat(i, 1), deplon(i, 1), pih, target_lat, intersect_lon(i, j))
          else
            call intersect_points(deplon(i, nlat), deplat(i, nlat), deplon(i, nlat), -pih, target_lat, intersect_lon(i, j))
          endif
        endif
      enddo
    enddo

    do i = 1, nlon
      sphere_dist(:) = 0.0d0
      do j = 2, nlat
        sphere_dist(j) = sphere_dist(j - 1) + orthodrome(deplon(i, j-1), deplat(i, j-1), deplon(i, j), deplat(i, j))
      end do
      do j = 1, nlat
        k = search_bisection(deplat(i, :), latitudes(j))

        if (k == 0) then
          sphere_loc(j) = -orthodrome(intersect_lon(i, j), latitudes(j), deplon(i, 1), deplat(i, 1))
        else
          sphere_loc(j) = sphere_dist(k) + orthodrome(intersect_lon(i, j), latitudes(j), deplon(i, k), deplat(i, k))
        endif
      enddo

      circle_lat_dist(:) = sphere_dist(:)
      gphi_old_circle(:) = gphi_old(i, :)
      call interpolate_spline(nlat-1, circle_lat_dist, gphi_old_circle, sphere_loc, intersect(i, :))
    end do

    do j = 1, nlat
      call sort_two_array(intersect_lon(:, j), intersect(:, j), intersect_lon_sort, intersect_sort)
      insert_lon_circle(1 : nlon) = intersect_lon_sort(:)
      insert_circle(1 : nlon) = intersect_sort(:)
      insert_lon_circle(0) = intersect_lon_sort(nlon) - 2.0d0*pi
      insert_lon_circle(nlon + 1) = intersect_lon_sort(1) + 2.0d0 * pi
      insert_circle(0) = intersect_sort(nlon)
      insert_circle(nlon + 1) = intersect_sort(1)
      call interpolate_spline(nlon+1, insert_lon_circle, insert_circle, lon, gphi(:, j))
    enddo

  end subroutine cascade_interpolate

  subroutine nijihouteishiki(a1, b1, c1, ans1, ans2)
    implicit none
    real(8), intent(in) :: a1, b1, c1
    real(8), intent(out) :: ans1, ans2

    if (b1*b1 - 4.0d0*a1*c1 < 0.0d0) then
      write(*, *) "2次方程式error!"
      return
    endif

    ans1 = (-b1 + sqrt(b1**2 - 4.0d0*a1*c1)) / (2.0d0 * a1)
    ans2 = (-b1 - sqrt(b1**2 - 4.0d0*a1*c1)) / (2.0d0 * a1)
  end subroutine nijihouteishiki

  subroutine intersect_points(lon1, lat1, lon2, lat2, target_lat, lon_ans)
    use sphere_module, only: lonlat2xyz, orthodrome
    use math_module, only: pi2=>math_pi2
    implicit none
    real(8), intent(in) :: lon1, lat1, lon2, lat2, target_lat
    real(8), intent(out) :: lon_ans

    real(8) :: x1, y1, z1, x2, y2, z2
    real(8) :: a2, b2, c2, ansx1, ansx2, ansy1, ansy2, lon_ans1, lon_ans2

    call lonlat2xyz(lon1, lat1, x1, y1, z1)
    call lonlat2xyz(lon2, lat2, x2, y2, z2)

    a2 = y1*z2 - y2*z1
    b2 = z1*x2 - z2*x1
    c2 = x1*y2 - x2*y1
    call nijihouteishiki(a2**2+b2**2, 2.0d0*a2*c2*sin(target_lat), &
      (b2**2+c2**2)*(sin(target_lat)**2) - b2**2, ansx1, ansx2)
    ansy1 = (-c2*sin(target_lat) - a2*ansx1) / b2
    ansy2 = (-c2*sin(target_lat) - a2*ansx2) / b2
    lon_ans1 = modulo(atan2(ansy1, ansx1)+pi2, pi2)
    lon_ans2 = modulo(atan2(ansy2, ansx2)+pi2, pi2)

    if (orthodrome(lon1, lat1, lon_ans1, target_lat) < orthodrome(lon1, lat1, lon_ans2, target_lat))then
      lon_ans = lon_ans1
    else
      lon_ans = lon_ans2
    endif

  end subroutine intersect_points
end module cascade_interpolation_module