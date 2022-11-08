module legendre_transform_module
  use glatwgt_module, only: glatwgt_calc
  use alf_module, only: pnm, alf_calc, alf_enm
  use fft_module, only: fft_init, fft_analysis, fft_synthesis
  implicit none
  private

  integer(8), private :: nlon, nlat, ntrunc
  complex(8), dimension(:,:), allocatable, private :: w
  real(8), dimension(:), allocatable, private :: glat, gwgt

  public :: legendre_init, legendre_analysis, legendre_synthesis, &
            legendre_synthesis_dlon, legendre_synthesis_dlat, legendre_synthesis_dlonlat

contains

  subroutine legendre_init(nx,ny,nt,lat,wgt)
    implicit none

    integer(8), intent(in) :: nx, ny, nt
    real(8), dimension(:), intent(inout) :: lat, wgt

    real(8), dimension(:,:), allocatable :: g

    nlon = nx
    nlat = ny
    ntrunc = nt
    allocate(glat(nlat),gwgt(nlat))
    call glatwgt_calc(glat,gwgt)
    lat(:) = asin(glat(:))
    wgt(:) = gwgt(:)
    call alf_calc(lat, ntrunc+1)

    allocate(g(nlon,nlat), w(0:nlon/2,nlat))
    call fft_init(g, w)
    deallocate(g)

  end subroutine legendre_init

  subroutine legendre_clean
    implicit none

    deallocate(glat,gwgt,w)

  end subroutine legendre_clean

  subroutine legendre_analysis(g,s)
    implicit none

    real(8), dimension(:,:), intent(in) :: g
    complex(8), dimension(0:,0:), intent(inout) :: s
    integer(8) :: n, m, nlat2, flg 
    complex(8) :: tmp

    call fft_analysis(g, w)

    nlat2 = nlat/2

    do m=0, ntrunc
      flg = 1
      do n=m, ntrunc
        tmp = (0.0d0, 0.0d0)
        s(n, m) = sum( gwgt(1:nlat2)  * pnm(1:nlat2, n, m) * &
          (w(m, 1:nlat2) + flg * w(m, nlat:nlat2+1:-1)) )
        flg = -flg
      end do
      flg = -flg
    end do

  end subroutine legendre_analysis

  subroutine legendre_synthesis(s,g)
    implicit none

    complex(8), dimension(0:,0:), intent(in) :: s
    real(8), dimension(:,:), intent(inout) :: g

    integer(8) :: n, m, j, jr, nlat2, flg
    real(8) :: nh, sh
    complex(8) :: ntmp, stmp

    nlat2 = nlat/2

    w = (0.0d0, 0.0d0)
    do j=1, nlat2
      jr = nlat - j + 1
      do m=0, ntrunc
        ntmp = (0.0d0, 0.0d0)
        stmp = ntmp
        flg = 1
        do n=m, ntrunc
          nh = pnm(j, n, m)
          sh = nh * flg
          ntmp = ntmp + nh * s(n, m)
          stmp = stmp + sh * s(n, m)
          flg = -flg
        end do
        w(m, j) = ntmp
        w(m, jr) = stmp
        flg = -flg
      end do
    end do

    call fft_synthesis(w,g)

  end subroutine legendre_synthesis

  subroutine legendre_synthesis_dlon(s,g)
    implicit none

    complex(8), dimension(0:,0:), intent(in) :: s
    real(8), dimension(:,:), intent(inout) :: g

    integer(8) :: n, m, j, jr, nlat2, flg
    real(8) :: nh, sh
    complex(8) :: ntmp, stmp, im

    nlat2 = nlat/2

    w = (0.0d0, 0.0d0)
    do j=1, nlat2
      jr = nlat - j + 1
      do m=0, ntrunc
        im = dcmplx(0.0d0, m)
        ntmp = (0.0d0, 0.0d0)
        stmp = ntmp
        flg = 1
        do n=m, ntrunc
          nh = pnm(j, n, m)
          sh = nh * flg
          ntmp = ntmp + im * nh * s(n, m)
          stmp = stmp + im * sh * s(n, m)
          flg = -flg
        end do
        w(m, j) = ntmp
        w(m, jr) = stmp
        flg = -flg
      end do
    end do

    call fft_synthesis(w,g)

  end subroutine legendre_synthesis_dlon

  subroutine legendre_synthesis_dlat(s,g)
    implicit none

    complex(8), dimension(0:,0:), intent(in) :: s
    real(8), dimension(:,:), intent(inout) :: g

    integer(8) :: n, m, j, jr, nlat2, flg
    real(8) :: nh, sh
    complex(8) :: ntmp, stmp

    nlat2 = nlat/2

    w = (0.0d0, 0.0d0)
    do j=1, nlat2
      jr = nlat - j + 1
      do m=0, ntrunc
        ntmp = (0.0d0, 0.0d0)
        stmp = ntmp
        flg = -1
        do n=m, ntrunc
          nh = (n+1)*alf_enm(n,m)*pnm(j,n-1,m)-n*alf_enm(n+1,m)*pnm(j,n+1,m)
          sh = nh * flg
          ntmp = ntmp + nh * s(n, m)
          stmp = stmp + sh * s(n, m)
          flg = -flg
        end do
        w(m, j) = ntmp
        w(m, jr) = stmp
        flg = -flg
      end do
    end do

    call fft_synthesis(w,g)

  end subroutine legendre_synthesis_dlat

  subroutine legendre_synthesis_dlonlat(s,g)
    implicit none

    complex(8), dimension(0:,0:), intent(in) :: s
    real(8), dimension(:,:), intent(inout) :: g

    integer(8) :: n, m, j, jr, nlat2, flg
    real(8) :: nh, sh
    complex(8) :: ntmp, stmp, im

    nlat2 = nlat/2

    w = (0.0d0, 0.0d0)
    do j=1, nlat2
      jr = nlat - j + 1
      do m=0, ntrunc
        im = dcmplx(0.0d0, m)
        ntmp = (0.0d0, 0.0d0)
        stmp = ntmp
        flg = -1
        do n=m, ntrunc
          nh = (n+1)*alf_enm(n,m)*pnm(j,n-1,m)-n*alf_enm(n+1,m)*pnm(j,n+1,m)
          sh = nh * flg
          ntmp = ntmp + im * nh * s(n, m)
          stmp = stmp + im * sh * s(n, m)
          flg = -flg
        end do
        w(m, j) = ntmp
        w(m, jr) = stmp
        flg = -flg
      end do
    end do

    call fft_synthesis(w,g)

  end subroutine legendre_synthesis_dlonlat

end module legendre_transform_module