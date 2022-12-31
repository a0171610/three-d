module spline_interpolate_module

  private

  public :: interpolate_spline, interpolate_spline_1index

  contains

    subroutine interpolate_spline_1index(n, x, f, ansx, ansf)
      implicit none
      integer(8), intent(in) :: n
      real(8), intent(in) :: x(1:n), f(1:n), ansx(:)
      real(8), intent(out) :: ansf(:)

      real(8) :: x0(0:n-1), f0(0:n-1)
      x0(:) = x(:)
      f0(:) = f(:)

      call interpolate_spline(n-1, x0, f0, ansx, ansf)
    end subroutine interpolate_spline_1index
  
    subroutine interpolate_spline(n, x, f, ansx, ansf)
      implicit none
      integer(8), intent(in) :: n
      real(8), intent(in) :: x(0 : n), f(0 : n), ansx(:)
      real(8), intent(out) :: ansf(:)
      integer(8) :: i, sz
  
      sz = size(ansx)
      call cubic_spline(0.0d0, 0.0d0, x, f, n, 1_8)
      do i = 1, sz
        call cubic_spline(ansx(i), ansf(i), x, f, n, 0_8)
      enddo
    end subroutine interpolate_spline
  
    subroutine cubic_spline(x, y, xp, yp, n, mode)
      implicit none
      real(8) x,y,xp(0:n),yp(0:n)
      real(8), allocatable, save :: x1(:),y1(:),a(:),b(:),c(:)
      real(8), allocatable :: ah(:),bh(:),ch(:),dh(:)
      real(8) h1,h2,x0
      integer(8) n,mode,i,i1,i2
      integer(8), save :: nmax
      if (mode == 1) then
        if (allocated(x1)) deallocate ( x1,y1,a,b,c )
        allocate ( x1(0:n), y1(0:n-1), a(0:n-1), b(0:n-1), c(0:n) )
        allocate ( ah(0:n), bh(0:n), ch(0:n), dh(0:n) )
        h1 = xp(1)-xp(0)
        bh(0) = 2*h1
        ch(0) = h1
        dh(0) = 3*(yp(1)-yp(0))
        do i = 1, n-1
          h1 = xp(i)-xp(i-1)
          h2 = xp(i+1)-xp(i)
          ah(i) = h2
          bh(i) = 2*(xp(i+1)-xp(i-1))
          ch(i) = h1
          dh(i) = 3*((yp(i)-yp(i-1))*h2/h1+(yp(i+1)-yp(i))*h1/h2)
        enddo
        h1 = xp(n)-xp(n-1)
        ah(n) = h1
        bh(n) = 2*h1
        dh(n) = 3*(yp(n)-yp(n-1))
        call tridiagonal_matrix(ah,bh,ch,dh,n+1,c)
        do i = 0, n-1
          h1 = xp(i+1)-xp(i);     h2 = h1*h1
          x1(i) = xp(i)
          y1(i) = yp(i)
          b(i) = (3*(yp(i+1)-yp(i))-(c(i+1)+2*c(i))*h1)/h2
          a(i) = (c(i+1)-c(i)-2*b(i)*h1)/(3*h2)
        enddo
        x1(n) = xp(n)
        nmax  = n
        deallocate ( ah, bh, ch, dh )
        return
      endif
      if (x <= x1(1)) then
        i1 = 0
      else if (x >= x1(nmax-1)) then
        i1 = nmax-1
        else
          i1 = 1;  i2 = nmax-1
        do while (i2-i1 > 1)
          i = (i1+i2)/2
          if (x < x1(i)) then
            i2 = i
          else
            i1 = i
          endif
        enddo
      endif
      x0 = x - x1(i1)
      y = ((a(i1)*x0 + b(i1))*x0 + c(i1))*x0 + y1(i1)
    end subroutine cubic_spline
  
    subroutine tridiagonal_matrix(a,b,c,d,n,x)
      implicit none
      real(8) a(n),b(n),c(n),d(n),x(n)
      real(8) G(n),H(n),den
      integer(8) i,n
      G(1) = -c(1)/b(1)
      H(1) = d(1)/b(1)
      do i = 2, n
        den  = 1/(b(i) + a(i)*G(i-1))
        G(i) = -c(i)*den
        H(i) = (d(i) - a(i)*H(i-1))*den
      enddo
      x(n) = H(n)
      do i = n-1, 1, -1
        x(i) = G(i)*x(i+1) + H(i)
      enddo
    end subroutine tridiagonal_matrix
  
  end module spline_interpolate_module
