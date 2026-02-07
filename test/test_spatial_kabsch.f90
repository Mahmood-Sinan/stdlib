module test_kabsch
    use testdrive, only : new_unittest, unittest_type, error_type, check, skip_test
    use stdlib_kinds
    use stdlib_math, only: all_close, is_close
    use stdlib_linalg, only: svd
    use stdlib_spatial
    implicit none

contains


    !> Collect all exported unit tests
    subroutine collect_suite(testsuite)
        !> Collection of tests
        type(unittest_type), allocatable, intent(out) :: testsuite(:)

        testsuite = [ &
            new_unittest('kabsch_real', test_kabsch_real), &
            new_unittest('kabsch_complex', test_kabsch_complex) &
        ]
    end subroutine

    subroutine test_kabsch_real(error)
        type(error_type), allocatable, intent(out) :: error
        block
            integer, parameter :: d = 3, N = 4
            real(sp) :: P_original(d, N), Q_original(d, N), P_recovered(d, N)
            real(sp) :: R_recovered(d, d), R_original(d, d)
            real(sp) :: t_recovered(d), t_original(d)
            real(sp) :: U(d,d), Vt(d,d), S(d)
            real(sp) :: c_recovered, c_original
            real(sp) :: rmsd

            integer :: i, j
            real(sp) :: r1

            call random_seed()

            ! Random reference points Q
            call random_number(Q_original)

            ! Random proper rotation matrix R_original
            ! Constructed via SVD: R = U * V^T
            call random_number(R_original)
            call svd(R_original, S, U, Vt)
            R_original = matmul(U, Vt)

            ! Random scale and translation
            call random_number(c_original)
            call random_number(t_original)

            ! Construct P = c*R*Q + t
            do j = 1, N
                P_original(:,j) = c_original * matmul(R_original, Q_original(:,j)) + t_original
            end do

            ! Call Kabsch–Umeyama
            call kabsch(P_original, Q_original, R_recovered, t_recovered, c_recovered, rmsd)

            ! Apply recovered transform
            do j = 1, N
                P_recovered(:,j) = c_recovered * matmul(R_recovered, Q_original(:,j)) + t_recovered
            end do

            call check(error, all_close(P_recovered, P_original), .true.)
            if (allocated(error)) return
            call check(error, all_close(R_recovered, R_original), .true.)
            if (allocated(error)) return
            call check(error, is_close(c_recovered, c_original), .true.)
            if (allocated(error)) return
            call check(error, all_close(t_recovered, t_original), .true.)
            if (allocated(error)) return
        end block
        block
            integer, parameter :: d = 3, N = 4
            real(dp) :: P_original(d, N), Q_original(d, N), P_recovered(d, N)
            real(dp) :: R_recovered(d, d), R_original(d, d)
            real(dp) :: t_recovered(d), t_original(d)
            real(dp) :: U(d,d), Vt(d,d), S(d)
            real(dp) :: c_recovered, c_original
            real(dp) :: rmsd

            integer :: i, j
            real(dp) :: r1

            call random_seed()

            ! Random reference points Q
            call random_number(Q_original)

            ! Random proper rotation matrix R_original
            ! Constructed via SVD: R = U * V^T
            call random_number(R_original)
            call svd(R_original, S, U, Vt)
            R_original = matmul(U, Vt)

            ! Random scale and translation
            call random_number(c_original)
            call random_number(t_original)

            ! Construct P = c*R*Q + t
            do j = 1, N
                P_original(:,j) = c_original * matmul(R_original, Q_original(:,j)) + t_original
            end do

            ! Call Kabsch–Umeyama
            call kabsch(P_original, Q_original, R_recovered, t_recovered, c_recovered, rmsd)

            ! Apply recovered transform
            do j = 1, N
                P_recovered(:,j) = c_recovered * matmul(R_recovered, Q_original(:,j)) + t_recovered
            end do

            call check(error, all_close(P_recovered, P_original), .true.)
            if (allocated(error)) return
            call check(error, all_close(R_recovered, R_original), .true.)
            if (allocated(error)) return
            call check(error, is_close(c_recovered, c_original), .true.)
            if (allocated(error)) return
            call check(error, all_close(t_recovered, t_original), .true.)
            if (allocated(error)) return
        end block
    end subroutine

    subroutine test_kabsch_complex(error)
        !> Error handling
        type(error_type), allocatable, intent(out) :: error
        block
            integer, parameter :: d = 3, N = 4
            complex(sp) :: P_original(d, N), Q_original(d, N), P_recovered(d, N)
            complex(sp) :: R_recovered(d, d), R_original(d, d)
            complex(sp) :: t_recovered(d), t_original(d)
            complex(sp) :: c_recovered, c_original
            complex(sp) :: U(d,d), Vt(d,d)
            real(sp) :: S(d)
            real(sp) :: rmsd

            integer :: i, j
            real(sp) :: r1, r2

            call random_seed()

            ! Random complex reference points Q
            do j = 1, N
                do i = 1, d
                    call random_number(r1)
                    call random_number(r2)
                    Q_original(i,j) = cmplx(r1, r2, kind=sp)
                end do
            end do

            ! Random complex unitary matrix R_original
            ! Constructed via SVD: R = U * V^H
            do i = 1, d
                do j = 1, d
                    call random_number(r1)
                    call random_number(r2)
                    R_original(i,j) = cmplx(r1, r2, kind=sp)
                end do
            end do
            call svd(R_original, S, U, Vt)
            R_original = matmul(U, Vt)

            ! Random scale and translation
            call random_number(r1)
            call random_number(r2)
            c_original = cmplx(r1,r2, kind=sp)

            do i = 1, d
                call random_number(r1)
                call random_number(r2)
                t_original(i) = cmplx(r1, r2, kind=sp)
            end do

            ! Construct P = c*R*Q + t
            do j = 1, N
                P_original(:,j) = c_original * matmul(R_original, Q_original(:,j)) + t_original
            end do

            ! Call Kabsch–Umeyama
            call kabsch(P_original, Q_original, R_recovered, t_recovered, c_recovered, rmsd)

            ! Apply recovered transform
            do j = 1, N
                P_recovered(:,j) = c_recovered * matmul(R_recovered, Q_original(:,j)) + t_recovered
            end do
            call check(error, all_close(P_recovered, P_original), .true.)
            if (allocated(error)) return
            call check(error, all_close(c_recovered * R_recovered, c_original * R_original), .true.)
            if (allocated(error)) return
            call check(error, all_close(t_recovered, t_original), .true.)
            if (allocated(error)) return
        end block
        block
            integer, parameter :: d = 3, N = 4
            complex(dp) :: P_original(d, N), Q_original(d, N), P_recovered(d, N)
            complex(dp) :: R_recovered(d, d), R_original(d, d)
            complex(dp) :: t_recovered(d), t_original(d)
            complex(dp) :: c_recovered, c_original
            complex(dp) :: U(d,d), Vt(d,d)
            real(dp) :: S(d)
            real(dp) :: rmsd

            integer :: i, j
            real(dp) :: r1, r2

            call random_seed()

            ! Random complex reference points Q
            do j = 1, N
                do i = 1, d
                    call random_number(r1)
                    call random_number(r2)
                    Q_original(i,j) = cmplx(r1, r2, kind=dp)
                end do
            end do

            ! Random complex unitary matrix R_original
            ! Constructed via SVD: R = U * V^H
            do i = 1, d
                do j = 1, d
                    call random_number(r1)
                    call random_number(r2)
                    R_original(i,j) = cmplx(r1, r2, kind=dp)
                end do
            end do
            call svd(R_original, S, U, Vt)
            R_original = matmul(U, Vt)

            ! Random scale and translation
            call random_number(r1)
            call random_number(r2)
            c_original = cmplx(r1,r2, kind=dp)

            do i = 1, d
                call random_number(r1)
                call random_number(r2)
                t_original(i) = cmplx(r1, r2, kind=dp)
            end do

            ! Construct P = c*R*Q + t
            do j = 1, N
                P_original(:,j) = c_original * matmul(R_original, Q_original(:,j)) + t_original
            end do

            ! Call Kabsch–Umeyama
            call kabsch(P_original, Q_original, R_recovered, t_recovered, c_recovered, rmsd)

            ! Apply recovered transform
            do j = 1, N
                P_recovered(:,j) = c_recovered * matmul(R_recovered, Q_original(:,j)) + t_recovered
            end do
            call check(error, all_close(P_recovered, P_original), .true.)
            if (allocated(error)) return
            call check(error, all_close(c_recovered * R_recovered, c_original * R_original), .true.)
            if (allocated(error)) return
            call check(error, all_close(t_recovered, t_original), .true.)
            if (allocated(error)) return
        end block
    end subroutine

end module


program tester
    use, intrinsic :: iso_fortran_env, only : error_unit
    use testdrive, only : run_testsuite, new_testsuite, testsuite_type
    use test_kabsch, only : collect_suite
    implicit none
    integer :: stat, is
    type(testsuite_type), allocatable :: testsuites(:)
    character(len=*), parameter :: fmt = '("#", *(1x, a))'

    stat = 0

    testsuites = [ &
        new_testsuite("kabsch", collect_suite) &
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
