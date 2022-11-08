module glatwgt_module
  ! calculates Gaussian points and weights.
    use math_module, only: pi=>math_pi
  !  implicit none
    private
  
  ! Source: Based on Swartztrauber (2003)
  ! Author: T. Enomoto
  ! Usage:
  !   Calculate Gaussian points and weights
  !     subroutine f_scgaus(x, w, jMax)
  !       glat(jMax)  cos(Gaussian colatitudes) = sin(Gaussian latitudes)
  !       gwt(jMax)  Gaussian weights
  ! History:
  !   TE 27 Apr 2003  Fixed a bug in setting SH colatitudes
  !   TE 28 Apr 2003  Ported for AFES
  !   TE 24 Apr 2003  Implemented Fourier-Legendre formulation.
  
  ! xacc_min  double minimum accuracy available
    real(8), private, parameter :: xacc_min = 1.0d-15
    integer(8), private :: jMid
    real(8), private, dimension(:), allocatable :: an
  
    private :: legendre_init, legendre_P, legendre_dP, legendre_clean, newton
    public :: glatwgt_calc
  
  contains
  
  ! Public procedures
  
    subroutine glatwgt_calc(glat,gwgt)
      implicit none
  
  ! calculates sin(Gaussian latitudes) between 1 and -1
  ! and Gaussian weights.
  ! NB. Gaussian colatitudes are used during calculation
  
      real(8), dimension(:), intent(inout) :: glat, gwgt
  
      integer(8) :: l, j, jMax
      real(8) :: guess, dpn, s
  
      jMax = size(glat)
      jMid = jMax/2
  
      call legendre_init()
  
      guess = 0.5d0*pi*(1.0d0-1.0d0/(jMax + 1))
      call newton(legendre_P, legendre_dP, guess, glat(jMid))
      guess = 3.0d0*glat(jMid) - pi
      call newton(legendre_P, legendre_dP, guess, glat(jMid-1))
      do l = jMid-2, 1, -1
        guess = 2*glat(l+1) - glat(l+2) 
        call newton(legendre_P, legendre_dP, guess, glat(l))
      end do
      do j = 1, jMid
        call legendre_dP(glat(j), dpn)
        gwgt(j) = (2.0d0*jMax + 1.0d0)/(dpn)**2
      end do
  
    s = sum(gwgt(1:jMid))
  ! *****
  
      gwgt(jMax:jMid+1:-1) = gwgt(1:jMid)
      glat(1:jMid) = cos(glat(1:jMid))
      glat(jMax:jMid+1:-1) = -glat(1:jMid)
  
      call legendre_clean()
    
    end subroutine glatwgt_calc
  
  ! private procedures
  
  subroutine legendre_init()
    implicit none
  
    real(8), parameter :: a0sq = 1.5d0
    integer(8) :: k, l, n
  
    n = jMid*2
    allocate(an(0:jMid))
    an(jMid) = a0sq
    do k=2, n
      an(jMid) = (1.0d0-1.0d0/(4.0d0*k*k))*an(jMid)
    end do
    an(jMid) = sqrt(an(jMid))
  
    do k=1, jMid
      l = 2*k
      an(jMid-k) = (l-1.0d0)*(2.0d0*n-l+2.0d0)/(l*(2.0d0*n-l+1.0d0)) * an(jMid-k+1) 
    end do
    an(0) = 0.5d0*an(0)
  
  end subroutine legendre_init
  
  subroutine legendre_P(theta, pn)
    implicit none
  
    real(8), intent(in) :: theta
    real(8), intent(out) :: pn
    integer(8) :: k, l
  
    pn = 0.0
    do l=0, jMid
      k=l*2 ! k = l*2 + 1 if n odd
      pn = pn + an(l) * cos(k*theta)
    end do
    
  end subroutine legendre_P
  
  subroutine legendre_dP(theta, dpn)
    implicit none
  
    real(8), intent(in) :: theta
    real(8), intent(out) :: dpn
    integer(8) :: k, l
  
    dpn = 0.0
    do l=1, jMid
      k=l*2 ! k = l*2 + 1 if n odd
      dpn = dpn - k * an(l) * sin(k*theta)
    end do
    
  end subroutine legendre_dP
  
  subroutine legendre_clean()
    implicit none
  
    deallocate(an)
  
  end subroutine legendre_clean
  
  subroutine newton(f, df, x0, x, tolerance)
    implicit none
  
  ! finds the root u
  
    interface
      subroutine f(x, f_result)
        real(8), intent(in) :: x
        real(8), intent(out) :: f_result
      end subroutine f
      subroutine df(x, f_result)
        real(8), intent(in) :: x
        real(8), intent(out) :: f_result
      end subroutine df
    end interface
    real(8), intent(in) :: x0
    real(8), intent(out) :: x
    real(8), optional, intent(in) :: tolerance
  
    integer(8), parameter :: newton_max = 500
    real(8) :: xacc = xacc_min
  
    real(8) :: y, dy
    integer(8) :: i
  
    if (present(tolerance)) then
      if (tolerance < xacc_min) then
        print *, "### Error in newton: tolerance too small"
        stop
      else
        xacc = tolerance
      end if
    else
      xacc = xacc_min
    end if
  
    x = x0
    do i = 1, newton_max
      call f(x, y)
      call df(x, dy)
      y = y/dy
      if (abs(y) < xacc) then
        return
      end if
      x = x - y
    end do
    print *, "### Error in newton : Too many refinement."
    
  end subroutine newton
  
  end module glatwgt_module