!! Ritchie (1987) advection test
program advection

  use planet_module, only: planet_init
  use grid_module, only: grid_init, grid_clean
  use time_module, only: model
  use euler_module, only: eulerian_init, eulerian_timeint, eulerian_clean
  use semilag_module, only: semilag_init, semilag_timeint, semilag_clean
  use nisl_module, only: nisl_init, nisl_timeint, nisl_clean
  use nisl_2step_module, only: nisl_2step_init, nisl_2step_timeint, nisl_2step_clean
  use direction_module, only: direction_init, direction_timeint, direction_clean
  use direction_2step_module, only: direction_2step_init, direction_2step_timeint, direction_2step_clean
  use field_module, only : field_init
  use analysis_module, only: error_log

  implicit none

  call planet_init()
  call grid_init()
  call field_init()
  select case(model)
    case("euler ")
      call eulerian_init()
      call eulerian_timeint()
      call eulerian_clean()
    case("slag  ")
      write(*,*) "hi"
      call semilag_init()
      call semilag_timeint()
      call semilag_clean()
    case("nisl")
      call nisl_init()
      call nisl_timeint()
      call nisl_clean()
    case("nisl_2step")
      call nisl_2step_init()
      call nisl_2step_timeint()
      call nisl_2step_clean()
    case("direction")
      call direction_init()
      call direction_timeint()
      call direction_clean()
    case("direction_2step")
      call direction_2step_init()
      call direction_2step_timeint()
      call direction_2step_clean()
    case default
      print *, "No matching model for", model
  end select
  write(*, *) "model = ", model
  call error_log()
  call grid_clean()

end program advection