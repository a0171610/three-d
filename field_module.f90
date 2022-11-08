module field_module
  use grid_module, only : lat, lon, nlat, nlon
  use planet_module, only : a => planet_radius
  use math_module, only : math_atan2
  implicit none
  real(8), allocatable, public :: X(:, :), Y(:, :)
contains

subroutine field_init()
  implicit none
  integer(8) :: i, j
  real(8) :: tmp1, tmp2
  allocate( X(nlon, nlat), Y(nlon, nlat) )

  do i = 1, nlon
    do j = 1, nlat
      call transform_coodinate(lon(i), lat(j), tmp1, tmp2)
      call latlon_to_XY(tmp1, tmp2, X(i,j), Y(i, j))
      X(i, j) = X(i, j) / a
      Y(i, j) = Y(i, j) / a
    enddo
  enddo
end subroutine field_init

subroutine latlon_to_XY(lambda, theta, X1, Y1) ! Ritchie(1986) 式(18)と(19)
  implicit none
  real(8), intent(in):: lambda, theta
  real(8), intent(out):: X1, Y1
  X1 = 2.0d0 * a * cos(theta) * cos(lambda) / (1.0d0 + sin(theta))
  Y1 = 2.0d0 * a * cos(theta) * sin(lambda) / (1.0d0 + sin(theta))
end subroutine latlon_to_XY

subroutine transform_coodinate(lambda, theta, lambda_new, theta_new)
  real(8), intent(in) :: theta, lambda
  real(8), intent(out) :: theta_new, lambda_new
  
  lambda_new = math_atan2(sin(lambda), sqrt(2.0d0) * (cos(lambda) - tan(theta)) / 2.0d0 )

  theta_new = sqrt(2.0d0) * (sin(theta) + cos(theta) * cos(lambda)) / 2.0d0
  theta_new = asin(theta_new)
end subroutine transform_coodinate

end module field_module