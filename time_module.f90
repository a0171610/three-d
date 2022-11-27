module time_module
  implicit none
  private

  integer(8), public, parameter ::  nstep = 200, hstep = 5
  !integer(8), public ::  nstep = 80, hstep = 5

  real(8), public, parameter :: deltat = 86400.0d0/dble(nstep)
  !real(8), public, parameter :: deltat = 86400.0d0/dble(nstep)

  character(len=15), public :: model = "nisl", imethod = "sph"
  character(len=15), public :: case = "hadley"
end module time_module
