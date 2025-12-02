program example_tridiagonal_determinant
    use stdlib_linalg_constants, only: dp
    use stdlib_specialmatrices, only: tridiagonal_dp_type, tridiagonal, determinant, dense
    use stdlib_linalg, only: det
    implicit none

    integer, parameter :: n = 4
    type(Tridiagonal_dp_type) :: A
    real(dp) :: dl(n-1), dv(n), du(n-1)
    real(dp) :: d

   call random_number(dl)
   call random_number(dv)
   call random_number(du)

    A = Tridiagonal(dl, dv, du)
    print*, dense(A)

    d = determinant(A)
    print*, "det(A):"
    print*, d

    print*, "det(A) using linalg:"
    print*, det(dense(A))

end program example_tridiagonal_determinant