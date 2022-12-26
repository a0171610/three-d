module sort_module
  use stdlib_sorting, only: sort, sort_index
  implicit none

contains

  subroutine sort_two_array(a, b, a_sort, b_sort)
    real(8), intent(in) :: a(:), b(:)
    real(8), intent(out) :: a_sort(:), b_sort(:)
    real(8), allocatable :: a_tmp(:)
    integer(8), allocatable :: index(:)
    integer(8) :: sz, i

    sz = size(a)
    allocate(a_tmp(sz), index(sz))

    a_tmp(:) = a(:)
    call sort_index(a_tmp, index)

    do i = 1, sz
      a_sort(i) = a(index(i))
      b_sort(i) = b(index(i))
    enddo
  end subroutine sort_two_array

  subroutine sort_three_array(a, b, c, a_sort, b_sort, c_sort)
    real(8), intent(in) :: a(:), b(:), c(:)
    real(8), intent(out) :: a_sort(:), b_sort(:), c_sort(:)
    real(8), allocatable :: a_tmp(:)
    integer(8), allocatable :: index(:)
    integer(8) :: sz, i

    sz = size(a)
    allocate(a_tmp(sz), index(sz))

    a_tmp(:) = a(:)
    call sort_index(a_tmp, index)

    do i = 1, sz
      a_sort(i) = a(index(i))
      b_sort(i) = b(index(i))
      c_sort(i) = c(index(i))
    enddo
  end subroutine sort_three_array
end module sort_module
