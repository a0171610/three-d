module upstream3d_module
  ! finds departure and mid-points
    use grid_module, only: latitudes=>lat, height
    use sphere_module, only: lonlat2xyz, uv2xyz
    use time_module, only: model
    implicit none
    private
  
    integer(8), public :: itermax = 50
    real(8), public :: small = 1.0d-10
  
    public :: find_points
  
  contains
  !! midlon, midlatはdt前の出発点(middle), deplon, deplatは2*dt前, uとvには時刻tのものを入れる(式(31)を参照)

    subroutine find_points(u, v, w, t, dt, midlon, midlat, deplon, deplat, deph)
      use math_module, only: pi2=>math_pi2
      use uv_module, only: calc_w, calc_ua, calc_va, calc_ud
      use uv_hadley_module, only: calc_w_hadley=>calc_w, calc_u, calc_v
      use time_module, only: case
      implicit none
  
      real(8), dimension(:, :, :), intent(in) :: u, v, w
      real(8), intent(in) :: t, dt
      real(8), dimension(:, :, :), intent(inout) :: midlon, midlat
      real(8), dimension(:, :, :), intent(out) :: deplon, deplat, deph
      real(8), parameter :: ht = 12000.0d0, eps = 1.0d-6, hs = 0.0d0

      integer(8) :: nx, ny, nz, i, j, k, step
      real(8) :: un, vn, wn, &     ! normalised velocity
                       bk, &          ! correction factor
                       xd, yd, zd, & ! Cartesian velocity
                       xg, yg, zg, & ! arival point in Cartesian coordinates
                       x0, y0, z0, & ! present point in Cartesian coordinates
                       x1, y1, z1, & ! updated point in Cartesian coordinates
                       lon, lat, h, err, &
                       h0, h1, lon0, lat0

      nx = size(u, 1)
      ny = size(u, 2)
      nz = size(u, 3)

      do j = 1, ny
        do i = 1, nx
          do k = 1, nz
            ! calculate initial values
            un = u(i, j, k)
            vn = v(i, j, k)
            wn = w(i, j, k)
            lon = pi2*dble(i-1)/dble(nx) ! calculate (lon,lat) from (i,j)
            lat = latitudes(j)
            h = height(k)

            call lonlat2xyz(lon, lat, xg, yg, zg) ! transform into Cartesian coordinates
            ! r = g as an initial point for the 1st time step, 最初midlat, midlonには格子点上の値が入っている
            call lonlat2xyz(midlon(i, j, k), midlat(i, j, k), x0, y0, z0)
            h0 = height(k)
            step = 1

            do
              lon0 = modulo(atan2(y0, x0) + pi2, pi2)
              lat0 = asin(z0)
              select case(case)
              case('hadley')
                h1 = h - dt * calc_w_hadley(lat0, h0, t)
              case('div')
                h1 = h - dt * calc_w(lon0, lat0, h0, t)
              case default
                print *, "No matching initial field"
              end select
             ! p1 = p - dt * calc_omega(lon0, lat0, p0, t)
              call uv2xyz(un,vn,lon,lat,xd,yd,zd) ! normalized Cartesian velocity
              ! correction factor
              bk = 1.0d0/sqrt(1.0d0+dt*dt*(xd*xd+yd*yd+zd*zd)-2.0d0*dt*(xd*xg+yd*yg+zd*zg))
              x1 =  bk*(xg - dt*xd) ! calculate new points
              y1 =  bk*(yg - dt*yd)
              z1 =  bk*(zg - dt*zd)
              ! calculate (lon,lat) from (x,y,z)
              lat = asin(z1)
              lon = modulo(atan2(y1,x1)+pi2,pi2)

              select case(case)
                case('hadley')
                  un = calc_u(lat)
                  vn = calc_v(lat, h1, t)
                case('div')
                  un = calc_ua(lon, lat, t) + calc_ud(lon, lat, h1, t)
                  vn = calc_va(lon, lat, t)
                case default
                  print *, "No matching initial field"
                stop
              end select

              err = sqrt((x1-x0)*(x1-x0)+(y1-y0)*(y1-y0)+(z1-z0)*(z1-z0)) ! calculate error

              x0 = x1 ! save as the current point
              y0 = y1
              z0 = z1
              h0 = h1
              step = step + 1
              if ((err<small).or.(step>itermax)) then
                exit
              end if
            end do

            midlon(i, j, k) = lon ! store as the mid-point
            midlat(i, j, k) = lat

            bk = 2.0d0*(x0*xg+y0*yg+z0*zg) ! calculate the departure point
            x1 = bk*x0 - xg
            y1 = bk*y0 - yg
            z1 = bk*z0 - zg

            deplon(i, j, k) = modulo(atan2(y1,x1)+pi2,pi2)
            deplat(i, j, k) = asin(z1)
            select case(case)
            case('hadley')
              deph(i, j, k) = h - 2.0d0 * dt * calc_w_hadley(lat, h1, t)
            case('div')
              deph(i, j, k) = h - 2.0d0 * dt * calc_w(lon, lat, h1, t)
            case default
              print *, "No matching initial field"
            stop
          end select
            if (deph(i, j, k) > ht - eps) then
              deph(i, j, k) = ht - eps
            endif
            if (deph(i, j, k) < hs - eps) then
              deph(i, j, k) = hs - eps
            endif

           ! write(*, *) height(k), deph(i, j, k)

          enddo
        end do
      end do
      
    end subroutine find_points

  end module upstream3d_module
