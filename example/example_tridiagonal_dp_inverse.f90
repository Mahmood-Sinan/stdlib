program example_tridiagonal_inverse
    use stdlib_linalg_constants, only: dp
    use stdlib_specialmatrices, only: tridiagonal_dp_type, tridiagonal, dense, inverse
    implicit none

    integer, parameter :: n = 2
    type(Tridiagonal_dp_type) :: A
    real(dp) :: dl(n-1), dv(n), du(n-1)
    real(dp) :: B(n, n), check(n, n)

    integer :: i, j
    dv = [0.0_dp, 1.0_dp]
    du = [1.0_dp]
    dl = [0.0_dp]

    A = Tridiagonal(dl, dv, du)
    B = inverse(A)

    print *, "Inverse of A:"
    do i = 1, n
        do j = 1, n
        print *, B(i,j)
        end do
        print *
    end do

    check = matmul(B, dense(A))
    print *, "Check: "
    do i = 1, n
        do j = 1, n
            if (i == j) then
                if (abs(check(i,j) - 1.0_dp) < 1.0e-12_dp) then
                    print *, 1.0_dp
                else
                    print *, check(i,j)
                end if
            else
                if (abs(check(i,j)) < 1.0e-12_dp) then
                    print *, 0.0_dp
                else
                    print *, check(i,j)
                end if
            end if
        end do
        print *
    end do

end program example_tridiagonal_inverse
