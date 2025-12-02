program example_tridiagonal_solve
    use stdlib_linalg_constants, only: dp
    use stdlib_specialmatrices, only: tridiagonal_dp_type, tridiagonal, dense, solve
    implicit none

    integer, parameter :: n = 4
    type(Tridiagonal_dp_type) :: A
    real(dp) :: dl(n-1), dv(n), du(n-1)
    real(dp) :: b(n), x(n)

    integer :: i

    ! -------------------------------------------
    ! Set up a simple tridiagonal system:
    !   4*x(i) + 1*x(i+1) + 1*x(i-1) = RHS
    !
    ! Example RHS chosen so that x(i)=1 is solution.
    ! -------------------------------------------
    call random_number(du)
    call random_number(dv)
    call random_number(dl)
    call random_number(b)

    ! -------------------------------------------
    ! Construct the Tridiagonal matrix type
    ! -------------------------------------------
    A = Tridiagonal(dl, dv, du)
    print*, "b: "
    print *, b
    call solve(A, b, x)

    print *, "Solution x:"
    do i = 1, n
        print *, x(i)
    end do
    print*, "Ax: "
    print *, matmul(dense(A), x)

end program example_tridiagonal_solve
