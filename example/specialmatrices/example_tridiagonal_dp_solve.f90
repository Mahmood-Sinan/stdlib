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
    du = [2.0_dp, 5.0_dp, 8.0_dp]
    dv = [1.0_dp, 4.0_dp, 7.0_dp, 10.0_dp]
    dl = [3.0_dp, 6.0_dp, 9.0_dp]

    b = [11.0_dp, 12.0_dp, 13.0_dp, 14.0_dp]

    ! -------------------------------------------
    ! Construct the Tridiagonal matrix type
    ! -------------------------------------------
    A = Tridiagonal(dl, dv, du)

    call solve(A, b, x)

    print *, "Solution x:"
    do i = 1, n
        print *, x(i)
    end do

end program example_tridiagonal_solve
