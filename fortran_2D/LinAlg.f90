module LinAlg
  implicit none

contains

  ! -------------------------------------------------------------------
  ! Solve A * X = B for X (n x m) using Gaussian elimination with
  ! partial pivoting.
  ! -------------------------------------------------------------------
  subroutine solve_linear_system(n, m, A_in, B_in, X_out)
    integer, intent(in)  :: n, m
    real(8), intent(in)  :: A_in(n, n)
    real(8), intent(in)  :: B_in(n, m)
    real(8), intent(out) :: X_out(n, m)

    real(8) :: Aug(n, n + m)
    real(8) :: tmp_row(n + m), factor, max_val
    integer :: i, j, k, pivot_row

    Aug(:, 1:n)     = A_in
    Aug(:, n+1:n+m) = B_in

    ! Forward elimination with partial pivoting
    do k = 1, n
      pivot_row = k
      max_val   = abs(Aug(k, k))
      do i = k + 1, n
        if (abs(Aug(i, k)) > max_val) then
          max_val   = abs(Aug(i, k))
          pivot_row = i
        end if
      end do

      if (pivot_row /= k) then
        tmp_row          = Aug(k, :)
        Aug(k, :)        = Aug(pivot_row, :)
        Aug(pivot_row, :) = tmp_row
      end if

      do i = k + 1, n
        factor = Aug(i, k) / Aug(k, k)
        Aug(i, k:n+m) = Aug(i, k:n+m) - factor * Aug(k, k:n+m)
      end do
    end do

    ! Back substitution
    do k = n, 1, -1
      do j = 1, m
        X_out(k, j) = Aug(k, n + j)
        do i = k + 1, n
          X_out(k, j) = X_out(k, j) - Aug(k, i) * X_out(i, j)
        end do
        X_out(k, j) = X_out(k, j) / Aug(k, k)
      end do
    end do
  end subroutine solve_linear_system

  ! -------------------------------------------------------------------
  ! Invert an (n x n) matrix: A_inv = inv(A)
  ! -------------------------------------------------------------------
  subroutine invert_matrix(n, A_in, A_inv)
    integer, intent(in)  :: n
    real(8), intent(in)  :: A_in(n, n)
    real(8), intent(out) :: A_inv(n, n)

    real(8) :: I_mat(n, n)
    integer :: i

    I_mat = 0.0_8
    do i = 1, n
      I_mat(i, i) = 1.0_8
    end do

    call solve_linear_system(n, n, A_in, I_mat, A_inv)
  end subroutine invert_matrix

  ! -------------------------------------------------------------------
  ! Kronecker product of two (n x n) matrices: C = kron(A, B)
  ! Matches MATLAB kron(A, B) of size (n^2 x n^2)
  ! -------------------------------------------------------------------
  subroutine kron_square(n, A, B, C)
    integer, intent(in)  :: n
    real(8), intent(in)  :: A(n, n), B(n, n)
    real(8), intent(out) :: C(n*n, n*n)

    integer :: i, j, p, q, row_idx, col_idx

    do j = 1, n
      do i = 1, n
        do q = 1, n
          do p = 1, n
            row_idx = (i - 1) * n + p
            col_idx = (j - 1) * n + q
            C(row_idx, col_idx) = A(i, j) * B(p, q)
          end do
        end do
      end do
    end do
  end subroutine kron_square

end module LinAlg
