module alf_module
  ! calculates normalised associated Legendre functions
    implicit none
    private
  
    real(8), public, dimension(:,:,:), allocatable :: pnm
  
    public :: alf_calc, alf_calc_old, alf_clean, alf_enm
  
  contains
  
    function alf_enm(n, m) result(e)
      implicit none
  
      integer(8), intent(in) :: n, m
      real(8) :: e, nn, mm
  
      nn = n*n
      mm = m*m
      e =  sqrt((nn - mm)/(4.0d0*nn - 1.0d0))
  
    end function alf_enm
  
    subroutine alf_calc_old(glat, mmax)
      implicit none
  
      real(8), dimension(:), intent(in) :: glat
      integer(8), intent(in) :: mmax
  
      integer(8) :: n, m, jmax, jmax2, ierr
  
      jmax = size(glat)
      jmax2 = jmax/2
      allocate(pnm(jmax2, -1:mmax, 0:mmax), stat = ierr)
      if (ierr > 0) then
        print *, "Allocation err in alf_calc"
        stop
      end if
      pnm = 0.0d0
  
      pnm(:,0,0) = sqrt(0.5d0)
      pnm(:,1,0) = sqrt(1.5d0)*glat(1:jmax2)
  
      do n=2, mmax
        pnm(:,n,0) = (glat(1:jmax2)*pnm(:,n-1,0_8)-alf_enm(n-1,0_8)*pnm(:,n-2,0_8))/alf_enm(n,0_8)
      end do
      do m=1, mmax
        pnm(:,m,m) = sqrt(0.5d0*(2.0d0*m+1.0d0)/m*(1.0d0-glat(1:jmax2)**2))*pnm(:,m-1,m-1)
        pnm(:,m+1,m) = sqrt(2.0d0*m+3.0d0)*glat(1:jmax2)*pnm(:,m,m)
        do n=m+2, mmax
          pnm(:,n,m) = (glat(1:jmax2)*pnm(:,n-1,m) - alf_enm(n-1,m)*pnm(:,n-2,m)) / alf_enm(n,m)
        end do
      end do
  
    end subroutine alf_calc_old
    
    subroutine alf_calc(lat,mmax)
      implicit none
  
      real(8), dimension(:), intent(in) :: lat
      integer(8), intent(in) :: mmax
  
      integer(8) :: j, m, n, l, k, n2, nmod, jmax, ierr
      real(8) :: sqrt_nn1_rev, theta, pi2
      real(8), dimension(:,:), allocatable :: ank
      
      jmax = size(lat)
      allocate(pnm(jmax/2, -1:mmax, 0:mmax), stat = ierr)
      if (ierr > 0) then
        print *, "Allocation err in alf_calc"
        stop
      end if
      pnm = 0.0d0
      allocate(ank(2:mmax, 0:mmax/2))
      ank = 0.0d0
  
  ! calculate fourier coefficients for Pn
      ank(2,1) = 0.75d0*sqrt(2.5d0) 
      do n=3, mmax
        ank(n,n/2) = &
          sqrt(1.0d0-1.0d0/(4.0d0*n*n)) * ank(n-1,(n-1)/2)
      end do
      do n=2, mmax
        n2  = n/2
        do k=1, n2
          l = 2*k
          ank(n,n2-k) = (l-1.0d0)*(2.0d0*n-l+2.0d0)/&
            (l*(2.0d0*n-l+1.0d0)) * ank(n,n2-k+1)
        end do
        if (n==n2*2) then
          ank(n,0) = 0.5d0*ank(n,0)
        end if
      end do
  
  ! calculate Pnm
      pi2 = 0.5d0 * acos(-1.0d0)
  
      do j=1, jmax/2
  
        theta = pi2 - lat(j)
        pnm(:, 0, 0) = 1.0d0/sqrt(2.0d0)
  
  ! Pmm and Pm,n=m+1
        do m=1, mmax
          pnm(j, m, m) = sqrt(1.0d0 + 0.5d0/m) * sin(theta) * pnm(j, m-1, m-1)
        end do
        do m=1, mmax-1
          pnm(j, m+1, m) = sqrt(2.0d0 * m + 3.0d0) * cos(theta) * pnm(j, m, m)
        end do
  
  ! m = 0
        pnm(j, 1, 0) = sqrt(1.5d0)*cos(theta)
        do n=2, mmax
          nmod = n - n/2*2
          do l=0, n/2
            k = 2*l + nmod ! n even: k=2*l, n odd: k=2*l+1
            pnm(j, n, 0) = pnm(j, n, 0) + ank(n,l)*cos(k*theta)
          end do
        end do
  
  ! m = 1
        do n=3, mmax
          nmod = n - n/2*2
          sqrt_nn1_rev = 1.0d0/sqrt(n*(n+1.0d0))
          do l=0, n/2
            k = 2*l + nmod ! n even: k=2*l, n odd: k=2*l+1
            pnm(j, n, 1) = pnm(j, n, 1) + ank(n,l)*k*sqrt_nn1_rev*sin(k*theta)
          end do
        end do
  
  ! m > 1
        do m=2, mmax-2
          do n=m+2, mmax
            pnm(j, n, m) = ( &
              sqrt((2.0d0*n+1.0d0)/(2.0d0*n-3.0d0))* &
                ( sqrt((n+m-2.0d0)*(n+m-3.0d0)) * pnm(j, n-2, m-2) + &
                   sqrt((n-m)*(n-m-1.0d0)) * pnm(j, n-2, m) ) &
              - sqrt((n-m+1.0d0)*(n-m+2.0d0)) * pnm(j, n, m-2) &
              ) / sqrt((n+m-1.0d0)*(n+m))
          end do
        end do
  
      end do! j
  
      deallocate(ank)
  
    end subroutine alf_calc
  
    subroutine alf_clean()
      implicit none
  
      deallocate(pnm)
  
    end subroutine alf_clean
  
  end module alf_module