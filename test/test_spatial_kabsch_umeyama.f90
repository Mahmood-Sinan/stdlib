module test_kabsch_umeyama
    use testdrive, only : new_unittest, unittest_type, error_type, check, skip_test
    use stdlib_kinds
    use stdlib_math, only: all_close, is_close
    use stdlib_linalg, only: svd, det, eye
    use stdlib_spatial
    use stdlib_intrinsics, only: stdlib_dot_product_kahan
    use stdlib_constants
    implicit none

contains


    !> Collect all exported unit tests
    subroutine collect_suite(testsuite)
        !> Collection of tests
        type(unittest_type), allocatable, intent(out) :: testsuite(:)

        testsuite = [ &
            new_unittest('kabsch_umeyama_real', test_kabsch_umeyama_real), &
            new_unittest('kabsch_umeyama_complex', test_kabsch_umeyama_complex) &
        ]
    end subroutine

    subroutine generate_random_set_real_sp(P, Q, d, N, scale)
        integer, intent(in) :: d, N
        real(sp), intent(out) :: P(d, N), Q(d, N)
        logical, intent(in) :: scale

        ! Internal variables.
        real(sp) :: R(d, d)
        real(sp) :: U(d,d), Vt(d,d)
        real(sp) :: S(d)
        real(sp) :: t(d)
        real(sp) :: c
        integer :: i, j

        c = one_sp
        ! Random reference points Q
        call random_number(Q)

        ! Random proper rotation matrix R constructed via SVD: R = U * V^T
        call random_number(R)
        call svd(R, S, U, Vt)
        do i = 1,d
            do j = 1,d
                R(i,j) = stdlib_dot_product_kahan(U(i,:), Vt(:, j))
            end do
        end do
        if(det(R) < zero_sp) R(:, d) = -R(:, d)

        ! Random scale and translation
        if(scale) call random_number(c)
        call random_number(t)

        ! Construct P = c*R*Q + t
        do j = 1, N
            do i = 1, d
                P(i,j) = c * &
                    stdlib_dot_product_kahan(R(i,:), Q(:,j)) + &
                    t(i)
            end do
        end do
    end subroutine
    subroutine generate_random_set_real_dp(P, Q, d, N, scale)
        integer, intent(in) :: d, N
        real(dp), intent(out) :: P(d, N), Q(d, N)
        logical, intent(in) :: scale

        ! Internal variables.
        real(dp) :: R(d, d)
        real(dp) :: U(d,d), Vt(d,d)
        real(dp) :: S(d)
        real(dp) :: t(d)
        real(dp) :: c
        integer :: i, j

        c = one_dp
        ! Random reference points Q
        call random_number(Q)

        ! Random proper rotation matrix R constructed via SVD: R = U * V^T
        call random_number(R)
        call svd(R, S, U, Vt)
        do i = 1,d
            do j = 1,d
                R(i,j) = stdlib_dot_product_kahan(U(i,:), Vt(:, j))
            end do
        end do
        if(det(R) < zero_dp) R(:, d) = -R(:, d)

        ! Random scale and translation
        if(scale) call random_number(c)
        call random_number(t)

        ! Construct P = c*R*Q + t
        do j = 1, N
            do i = 1, d
                P(i,j) = c * &
                    stdlib_dot_product_kahan(R(i,:), Q(:,j)) + &
                    t(i)
            end do
        end do
    end subroutine
    subroutine generate_random_set_complex_csp(P, Q, d, N, scale)
        integer, intent(in) :: d, N
        complex(sp), intent(out) :: P(d, N), Q(d, N)
        logical, intent(in) :: scale

        ! Internal variables.
        complex(sp) :: R(d,d), t(d)
        complex(sp) :: U(d,d), Vt(d,d)
        real(sp) :: S(d)
        complex(sp) :: c
        real(sp) :: tempdn(d,N,2)
        real(sp) :: tempdd(d,d,2)
        real(sp) :: r1, r2
        integer :: i, j

        c = one_csp
        ! Random reference points Q
        call random_number(tempdn)
        Q = cmplx(tempdn(:,:,1), tempdn(:,:,2), kind=sp)

        ! Random proper rotation matrix R constructed via SVD: R = U * V^H
        call random_number(tempdd)
        R = cmplx(tempdd(:,:,1), tempdd(:,:,2), kind=sp)
        call svd(R, S, U, Vt)
        do i = 1,d
            do j = 1,d
                R(i,j) = stdlib_dot_product_kahan(U(i,:), Vt(:, j))
            end do
        end do

        ! Random scale and translation
        if(scale) then
            call random_number(r1)
            call random_number(r2)
            c = cmplx(r1, r2, kind=sp)
        end if
        call random_number(tempdd)
        t = cmplx(tempdd(:,1,1), tempdd(:,1,2), kind=sp)

        ! Construct P = c*R*Q + t
        do j = 1, N
            do i = 1, d
                P(i,j) = c * &
                    stdlib_dot_product_kahan(R(i,:), Q(:,j)) + &
                    t(i)
            end do
        end do
    end subroutine
    subroutine generate_random_set_complex_cdp(P, Q, d, N, scale)
        integer, intent(in) :: d, N
        complex(dp), intent(out) :: P(d, N), Q(d, N)
        logical, intent(in) :: scale

        ! Internal variables.
        complex(dp) :: R(d,d), t(d)
        complex(dp) :: U(d,d), Vt(d,d)
        real(dp) :: S(d)
        complex(dp) :: c
        real(dp) :: tempdn(d,N,2)
        real(dp) :: tempdd(d,d,2)
        real(dp) :: r1, r2
        integer :: i, j

        c = one_cdp
        ! Random reference points Q
        call random_number(tempdn)
        Q = cmplx(tempdn(:,:,1), tempdn(:,:,2), kind=dp)

        ! Random proper rotation matrix R constructed via SVD: R = U * V^H
        call random_number(tempdd)
        R = cmplx(tempdd(:,:,1), tempdd(:,:,2), kind=dp)
        call svd(R, S, U, Vt)
        do i = 1,d
            do j = 1,d
                R(i,j) = stdlib_dot_product_kahan(U(i,:), Vt(:, j))
            end do
        end do

        ! Random scale and translation
        if(scale) then
            call random_number(r1)
            call random_number(r2)
            c = cmplx(r1, r2, kind=dp)
        end if
        call random_number(tempdd)
        t = cmplx(tempdd(:,1,1), tempdd(:,1,2), kind=dp)

        ! Construct P = c*R*Q + t
        do j = 1, N
            do i = 1, d
                P(i,j) = c * &
                    stdlib_dot_product_kahan(R(i,:), Q(:,j)) + &
                    t(i)
            end do
        end do
    end subroutine

    subroutine test_kabsch_umeyama_real(error)
        type(error_type), allocatable, intent(out) :: error
        block
            integer, parameter :: d = 3, N = 8
            real(sp) :: P_original(d, N), P_recovered(d, N), Q_original(d, N)
            real(sp) :: R_recovered(d, d)
            real(sp) :: t_recovered(d)
            real(sp) :: c_recovered
            real(sp) :: rmsd
            real(sp) :: Identity(d,d)
            real(sp) :: zero_t(d)
            real(sp) :: W(N)

            integer :: i, j
            real(sp), parameter :: sptol = 100 * epsilon(one_sp)

            R_recovered = zero_sp
            t_recovered = zero_sp
            c_recovered = zero_sp
            rmsd = zero_sp
            Identity = eye(d,d)
            zero_t = zero_sp

            ! Random set of transformation
            call generate_random_set_real_sp(P_original, Q_original, d, N, .true.)
            ! Call Kabsch–Umeyama
            call kabsch_umeyama(P_original, Q_original, R_recovered, t_recovered, c_recovered, rmsd)
            ! Apply recovered transform
            do j = 1, N
                do i = 1, d
                    P_recovered(i,j) = c_recovered * &
                        stdlib_dot_product_kahan(R_recovered(i,:), Q_original(:,j)) + &
                        t_recovered(i)
                end do
            end do
            call check(error, all_close(P_recovered, P_original, abs_tol = sptol), .true.)
            if (allocated(error)) return
            call check(error, is_close(rmsd, zero_sp, abs_tol = sptol), .true.)
            if(allocated(error)) return

            ! Identity transformation (Q -> Q)
            call random_number(Q_original)
            ! Call Kabsch–Umeyama
            call kabsch_umeyama(Q_original, Q_original, R_recovered, t_recovered, c_recovered, rmsd)
            call check(error, all_close(R_recovered, Identity, abs_tol=sptol), .true.)
            if(allocated(error)) return
            call check(error, all_close(t_recovered, zero_t, abs_tol = sptol), .true.)
            if (allocated(error)) return
            call check(error, is_close(rmsd, zero_sp, abs_tol = sptol), .true.)
            if(allocated(error)) return

            ! Random set of transformation but with scale disabled.
            call generate_random_set_real_sp(P_original, Q_original, d, N, .false.)
            ! Call Kabsch–Umeyama
            call kabsch_umeyama(P_original, Q_original, R_recovered, t_recovered, c_recovered, rmsd, scale=.false.)
            ! Apply recovered transform
            do j = 1, N
                do i = 1, d
                    P_recovered(i,j) = c_recovered * &
                        stdlib_dot_product_kahan(R_recovered(i,:), Q_original(:,j)) + &
                        t_recovered(i)
                end do
            end do
            call check(error, all_close(P_recovered, P_original, abs_tol = sptol), .true.)
            if (allocated(error)) return
            call check(error, is_close(c_recovered, one_sp, abs_tol = sptol), .true.)
            if(allocated(error)) return
            call check(error, is_close(rmsd, zero_sp, abs_tol = sptol), .true.)
            if(allocated(error)) return

            ! Random set of transformation but with weights. Since P = c*R*Q+t, weights should not matter
            call random_number(W)
            call generate_random_set_real_sp(P_original, Q_original, d, N, .true.)
            ! Call Kabsch–Umeyama
            call kabsch_umeyama(P_original, Q_original, R_recovered, t_recovered, c_recovered, rmsd, W=W)
            ! Apply recovered transform
            do j = 1, N
                do i = 1, d
                    P_recovered(i,j) = c_recovered * &
                        stdlib_dot_product_kahan(R_recovered(i,:), Q_original(:,j)) + &
                        t_recovered(i)
                end do
            end do
            call check(error, is_close(det(R_recovered), one_sp, abs_tol=sptol), .true.)
            if(allocated(error)) return
            call check(error, all_close(P_recovered, P_original, abs_tol = sptol), .true.)
            if (allocated(error)) return
            call check(error, is_close(rmsd, zero_sp, abs_tol = sptol), .true.)
            if(allocated(error)) return
        end block
        block
            integer, parameter :: d = 3, N = 8
            real(dp) :: P_original(d, N), P_recovered(d, N), Q_original(d, N)
            real(dp) :: R_recovered(d, d)
            real(dp) :: t_recovered(d)
            real(dp) :: c_recovered
            real(dp) :: rmsd
            real(dp) :: Identity(d,d)
            real(dp) :: zero_t(d)
            real(dp) :: W(N)

            integer :: i, j
            real(dp), parameter :: dptol = 100 * epsilon(one_dp)

            R_recovered = zero_dp
            t_recovered = zero_dp
            c_recovered = zero_dp
            rmsd = zero_dp
            Identity = eye(d,d)
            zero_t = zero_dp

            ! Random set of transformation
            call generate_random_set_real_dp(P_original, Q_original, d, N, .true.)
            ! Call Kabsch–Umeyama
            call kabsch_umeyama(P_original, Q_original, R_recovered, t_recovered, c_recovered, rmsd)
            ! Apply recovered transform
            do j = 1, N
                do i = 1, d
                    P_recovered(i,j) = c_recovered * &
                        stdlib_dot_product_kahan(R_recovered(i,:), Q_original(:,j)) + &
                        t_recovered(i)
                end do
            end do
            call check(error, all_close(P_recovered, P_original, abs_tol = dptol), .true.)
            if (allocated(error)) return
            call check(error, is_close(rmsd, zero_dp, abs_tol = dptol), .true.)
            if(allocated(error)) return

            ! Identity transformation (Q -> Q)
            call random_number(Q_original)
            ! Call Kabsch–Umeyama
            call kabsch_umeyama(Q_original, Q_original, R_recovered, t_recovered, c_recovered, rmsd)
            call check(error, all_close(R_recovered, Identity, abs_tol=dptol), .true.)
            if(allocated(error)) return
            call check(error, all_close(t_recovered, zero_t, abs_tol = dptol), .true.)
            if (allocated(error)) return
            call check(error, is_close(rmsd, zero_dp, abs_tol = dptol), .true.)
            if(allocated(error)) return

            ! Random set of transformation but with scale disabled.
            call generate_random_set_real_dp(P_original, Q_original, d, N, .false.)
            ! Call Kabsch–Umeyama
            call kabsch_umeyama(P_original, Q_original, R_recovered, t_recovered, c_recovered, rmsd, scale=.false.)
            ! Apply recovered transform
            do j = 1, N
                do i = 1, d
                    P_recovered(i,j) = c_recovered * &
                        stdlib_dot_product_kahan(R_recovered(i,:), Q_original(:,j)) + &
                        t_recovered(i)
                end do
            end do
            call check(error, all_close(P_recovered, P_original, abs_tol = dptol), .true.)
            if (allocated(error)) return
            call check(error, is_close(c_recovered, one_dp, abs_tol = dptol), .true.)
            if(allocated(error)) return
            call check(error, is_close(rmsd, zero_dp, abs_tol = dptol), .true.)
            if(allocated(error)) return

            ! Random set of transformation but with weights. Since P = c*R*Q+t, weights should not matter
            call random_number(W)
            call generate_random_set_real_dp(P_original, Q_original, d, N, .true.)
            ! Call Kabsch–Umeyama
            call kabsch_umeyama(P_original, Q_original, R_recovered, t_recovered, c_recovered, rmsd, W=W)
            ! Apply recovered transform
            do j = 1, N
                do i = 1, d
                    P_recovered(i,j) = c_recovered * &
                        stdlib_dot_product_kahan(R_recovered(i,:), Q_original(:,j)) + &
                        t_recovered(i)
                end do
            end do
            call check(error, is_close(det(R_recovered), one_dp, abs_tol=dptol), .true.)
            if(allocated(error)) return
            call check(error, all_close(P_recovered, P_original, abs_tol = dptol), .true.)
            if (allocated(error)) return
            call check(error, is_close(rmsd, zero_dp, abs_tol = dptol), .true.)
            if(allocated(error)) return
        end block
    end subroutine

    subroutine test_kabsch_umeyama_complex(error)
        !> Error handling
        type(error_type), allocatable, intent(out) :: error
        block
            integer, parameter :: d = 3, N = 8
            complex(sp) :: P_original(d, N), Q_original(d, N), P_recovered(d, N)
            complex(sp) :: R_recovered(d, d)
            complex(sp) :: t_recovered(d)
            complex(sp) :: c_recovered
            real(sp) :: rmsd
            complex(sp) :: Identity(d,d)
            complex(sp) :: zero_t(d)
            real(sp) :: W(N)

            integer :: i, j
            real(sp), parameter :: sptol = 100 * epsilon(one_sp)

            R_recovered = zero_csp
            t_recovered = zero_csp
            c_recovered = zero_csp
            rmsd = zero_sp
            Identity = eye(d,d)
            zero_t = zero_csp

            ! Random set of transformation
            call generate_random_set_complex_csp(P_original, Q_original, d, N, .true.)
            ! Call Kabsch–Umeyama
            call kabsch_umeyama(P_original, Q_original, R_recovered, t_recovered, c_recovered, rmsd)
            ! Apply recovered transform
            do j = 1, N
                do i = 1, d
                    P_recovered(i,j) = c_recovered * &
                        stdlib_dot_product_kahan(conjg(R_recovered(i,:)), Q_original(:,j)) + &
                        t_recovered(i)
                end do
            end do
            call check(error, is_close(abs(det(R_recovered)), one_sp, abs_tol=sptol), .true.)
            if(allocated(error)) return
            call check(error, all_close(P_recovered, P_original, abs_tol = sptol), .true.)
            if (allocated(error)) return
            call check(error, is_close(rmsd, zero_sp, abs_tol = sptol), .true.)
            if(allocated(error)) return

            ! Identity transformation (Q -> Q)
            call generate_random_set_complex_csp(P_original, Q_original, d, N, .true.)
            ! Call Kabsch–Umeyama
            call kabsch_umeyama(Q_original, Q_original, R_recovered, t_recovered, c_recovered, rmsd)
            call check(error, all_close(R_recovered, Identity, abs_tol=sptol), .true.)
            if(allocated(error)) return
            call check(error, all_close(t_recovered, zero_t, abs_tol = sptol), .true.)
            if (allocated(error)) return
            call check(error, is_close(rmsd, zero_sp, abs_tol = sptol), .true.)
            if(allocated(error)) return

            ! Random set of transformation but with scale disabled.
            call generate_random_set_complex_csp(P_original, Q_original, d, N, .false.)
            ! Call Kabsch–Umeyama
            call kabsch_umeyama(P_original, Q_original, R_recovered, t_recovered, c_recovered, rmsd,scale=.false.)
            ! Apply recovered transform
            do j = 1, N
                do i = 1, d
                    P_recovered(i,j) = c_recovered * &
                        stdlib_dot_product_kahan(conjg(R_recovered(i,:)), Q_original(:,j)) + &
                        t_recovered(i)
                end do
            end do
            call check(error, is_close(abs(det(R_recovered)), one_sp, abs_tol=sptol), .true.)
            if(allocated(error)) return
            call check(error, all_close(P_recovered, P_original, abs_tol = sptol), .true.)
            if (allocated(error)) return
            call check(error, is_close(c_recovered, one_csp, abs_tol = sptol), .true.)
            if(allocated(error)) return
            call check(error, is_close(rmsd, zero_sp, abs_tol = sptol), .true.)
            if(allocated(error)) return

            ! Random set of transformation but with weights. Since P = c*R*Q+t, weights should not matter
            call random_number(W)
            call generate_random_set_complex_csp(P_original, Q_original, d, N, .true.)
            ! Call Kabsch–Umeyama
            call kabsch_umeyama(P_original, Q_original, R_recovered, t_recovered, c_recovered, rmsd, W=W)
            ! Apply recovered transform
            do j = 1, N
                do i = 1, d
                    P_recovered(i,j) = c_recovered * &
                        stdlib_dot_product_kahan(conjg(R_recovered(i,:)), Q_original(:,j)) + &
                        t_recovered(i)
                end do
            end do
            call check(error, is_close(abs(det(R_recovered)), one_sp, abs_tol=sptol), .true.)
            if(allocated(error)) return
            call check(error, all_close(P_recovered, P_original, abs_tol = sptol), .true.)
            if (allocated(error)) return
            call check(error, is_close(rmsd, zero_sp, abs_tol = sptol), .true.)
            if(allocated(error)) return
        end block
        block
            integer, parameter :: d = 3, N = 8
            complex(dp) :: P_original(d, N), Q_original(d, N), P_recovered(d, N)
            complex(dp) :: R_recovered(d, d)
            complex(dp) :: t_recovered(d)
            complex(dp) :: c_recovered
            real(dp) :: rmsd
            complex(dp) :: Identity(d,d)
            complex(dp) :: zero_t(d)
            real(dp) :: W(N)

            integer :: i, j
            real(dp), parameter :: dptol = 100 * epsilon(one_dp)

            R_recovered = zero_cdp
            t_recovered = zero_cdp
            c_recovered = zero_cdp
            rmsd = zero_dp
            Identity = eye(d,d)
            zero_t = zero_cdp

            ! Random set of transformation
            call generate_random_set_complex_cdp(P_original, Q_original, d, N, .true.)
            ! Call Kabsch–Umeyama
            call kabsch_umeyama(P_original, Q_original, R_recovered, t_recovered, c_recovered, rmsd)
            ! Apply recovered transform
            do j = 1, N
                do i = 1, d
                    P_recovered(i,j) = c_recovered * &
                        stdlib_dot_product_kahan(conjg(R_recovered(i,:)), Q_original(:,j)) + &
                        t_recovered(i)
                end do
            end do
            call check(error, is_close(abs(det(R_recovered)), one_dp, abs_tol=dptol), .true.)
            if(allocated(error)) return
            call check(error, all_close(P_recovered, P_original, abs_tol = dptol), .true.)
            if (allocated(error)) return
            call check(error, is_close(rmsd, zero_dp, abs_tol = dptol), .true.)
            if(allocated(error)) return

            ! Identity transformation (Q -> Q)
            call generate_random_set_complex_cdp(P_original, Q_original, d, N, .true.)
            ! Call Kabsch–Umeyama
            call kabsch_umeyama(Q_original, Q_original, R_recovered, t_recovered, c_recovered, rmsd)
            call check(error, all_close(R_recovered, Identity, abs_tol=dptol), .true.)
            if(allocated(error)) return
            call check(error, all_close(t_recovered, zero_t, abs_tol = dptol), .true.)
            if (allocated(error)) return
            call check(error, is_close(rmsd, zero_dp, abs_tol = dptol), .true.)
            if(allocated(error)) return

            ! Random set of transformation but with scale disabled.
            call generate_random_set_complex_cdp(P_original, Q_original, d, N, .false.)
            ! Call Kabsch–Umeyama
            call kabsch_umeyama(P_original, Q_original, R_recovered, t_recovered, c_recovered, rmsd,scale=.false.)
            ! Apply recovered transform
            do j = 1, N
                do i = 1, d
                    P_recovered(i,j) = c_recovered * &
                        stdlib_dot_product_kahan(conjg(R_recovered(i,:)), Q_original(:,j)) + &
                        t_recovered(i)
                end do
            end do
            call check(error, is_close(abs(det(R_recovered)), one_dp, abs_tol=dptol), .true.)
            if(allocated(error)) return
            call check(error, all_close(P_recovered, P_original, abs_tol = dptol), .true.)
            if (allocated(error)) return
            call check(error, is_close(c_recovered, one_cdp, abs_tol = dptol), .true.)
            if(allocated(error)) return
            call check(error, is_close(rmsd, zero_dp, abs_tol = dptol), .true.)
            if(allocated(error)) return

            ! Random set of transformation but with weights. Since P = c*R*Q+t, weights should not matter
            call random_number(W)
            call generate_random_set_complex_cdp(P_original, Q_original, d, N, .true.)
            ! Call Kabsch–Umeyama
            call kabsch_umeyama(P_original, Q_original, R_recovered, t_recovered, c_recovered, rmsd, W=W)
            ! Apply recovered transform
            do j = 1, N
                do i = 1, d
                    P_recovered(i,j) = c_recovered * &
                        stdlib_dot_product_kahan(conjg(R_recovered(i,:)), Q_original(:,j)) + &
                        t_recovered(i)
                end do
            end do
            call check(error, is_close(abs(det(R_recovered)), one_dp, abs_tol=dptol), .true.)
            if(allocated(error)) return
            call check(error, all_close(P_recovered, P_original, abs_tol = dptol), .true.)
            if (allocated(error)) return
            call check(error, is_close(rmsd, zero_dp, abs_tol = dptol), .true.)
            if(allocated(error)) return
        end block
    end subroutine

end module


program tester
    use, intrinsic :: iso_fortran_env, only : error_unit
    use testdrive, only : run_testsuite, new_testsuite, testsuite_type
    use test_kabsch_umeyama, only : collect_suite
    implicit none
    integer :: stat, is
    type(testsuite_type), allocatable :: testsuites(:)
    character(len=*), parameter :: fmt = '("#", *(1x, a))'

    stat = 0

    testsuites = [ &
        new_testsuite("kabsch_umeyama", collect_suite) &
        ]

    do is = 1, size(testsuites)
        write(error_unit, fmt) "Testing:", testsuites(is)%name
        call run_testsuite(testsuites(is)%collect, error_unit, stat)
    end do

    if (stat > 0) then
        write(error_unit, '(i0, 1x, a)') stat, "test(s) failed!"
        error stop
    end if
end program
