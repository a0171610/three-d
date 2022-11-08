module fft_module
  implicit none
  include "fftw3.f"
  private

  integer(8), private :: plan_forward, plan_backward

  public :: fft_init, fft_clean, fft_analysis, fft_synthesis

contains

  subroutine fft_init(g, w)
    implicit none

    real(8), dimension(:,:), intent(in) :: g 
    complex(8), dimension(0:,:), intent(in) :: w 
    integer(8) :: i, j

    i = size(g,1)
    j = size(g,2)

    call dfftw_plan_dft_r2c_1d(plan_forward, i, g(:,1), w(:,1), FFTW_ESTIMATE)
    call dfftw_plan_dft_c2r_1d(plan_backward, i, w(:,1), g(:,1), FFTW_ESTIMATE)

  end subroutine fft_init

  subroutine fft_clean()
    implicit none

    call dfftw_destroy_plan(plan_forward)
    call dfftw_destroy_plan(plan_backward)

  end subroutine fft_clean

  subroutine fft_analysis(g, w)
    implicit none

    real(8), dimension(:,:), intent(in) :: g 
    complex(8), dimension(0:,:), intent(inout) :: w 
    integer(8) :: j

    do j=1, size(g,2)
      call dfftw_execute_dft_r2c(plan_forward, g(:,j), w(:,j))
    end do
    w = w/size(g,1)

  end subroutine fft_analysis

  subroutine fft_synthesis(w,g)
    implicit none

    complex(8), dimension(0:,:), intent(in) :: w 
    real(8), dimension(:,:), intent(inout) :: g 
    integer(8) :: j

    do j=1, size(g,2)
      call dfftw_execute_dft_c2r(plan_backward, w(:,j), g(:,j))
    end do

  end subroutine fft_synthesis

end module fft_module