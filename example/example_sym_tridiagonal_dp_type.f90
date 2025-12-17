program example_sym_tridiagonal_matrix
    use stdlib_linalg_constants, only: dp
    use stdlib_specialmatrices
    implicit none

    integer, parameter :: n = 4
    type(sym_tridiagonal_cdp_type) :: A
    complex(dp) :: du(n- 1), dv(n)

    ! Generate random tridiagonal elements.
    dv = [(2.0_dp, 0.0_dp), (2.0_dp, 0.0_dp), (2.0_dp, 0.0_dp), (2.0_dp, 0.0_dp)]
    du = [(-1.0_dp, 1.0_dp), (-1.0_dp, -2.0_dp), (0.5_dp, -0.5_dp)]

    ! Create the corresponding symmetrical tridiagonal matrix.
    A = sym_tridiagonal(du, dv)
    
end program example_sym_tridiagonal_matrix
