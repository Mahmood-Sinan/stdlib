program example_tridiagonal_solve
    use stdlib_linalg_constants, only: dp
    use stdlib_specialmatrices, only: tridiagonal_dp_type, tridiagonal, dense, solve
    implicit none

    integer, parameter :: n = 4, nrhs = 2
    type(Tridiagonal_dp_type) :: A
    real(dp) :: dl(n-1), dv(n), du(n-1)
    real(dp) :: b(n, nrhs), x(n, nrhs)

    call random_number(du)
    call random_number(dv)
    call random_number(dl)
    call random_number(b)

    A = Tridiagonal(dl, dv, du)
    
    print *, "b:"
    print *, b
    call solve(A, b, x)

    print *, "Solution x:"
    print *, x
    
    print*, "Ax: "
    print *, matmul(dense(A), x)

end program example_tridiagonal_solve
