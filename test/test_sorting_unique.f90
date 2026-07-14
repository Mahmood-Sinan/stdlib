module test_sorting_unique
    use testdrive, only : new_unittest, unittest_type, error_type, check, skip_test
    use stdlib_kinds
    use stdlib_sorting, only: unique
    use stdlib_math, only: all_close
    implicit none
contains

    !> Collect all exported unit tests
    subroutine collect_suite(testsuite)
        !> Collection of tests
        type(unittest_type), allocatable, intent(out) :: testsuite(:)

        testsuite = [ &
            new_unittest('int', test_int), &
            new_unittest('real', test_real) &
        ]
    end subroutine

    subroutine test_int(error)
        !> Error handling
        type(error_type), allocatable, intent(out) :: error
        block
            integer(int8), allocatable :: A(:)
            integer(int8), allocatable :: output(:)
            integer(int8), allocatable :: expected_sorted(:)
            integer(int8), allocatable :: expected_stable(:)

            ! Initialize matrix.
            allocate(A(0))
            output = unique(A, .true.)
            allocate(expected_sorted(0))
            call check(error, all(output==expected_sorted), .true.,&
                "Sorted(int8): empty array")
            if(allocated(error)) return
            output = unique(A, .false.)
            allocate(expected_stable(0))
            call check(error, all(output==expected_stable), .true.,&
                "Stable(int8): empty array")
            if(allocated(error)) return

            A = [1_int8]
            output = unique(A, .true.)
            expected_sorted = [1_int8]
            call check(error, all(output==expected_sorted), .true.,&
                "Sorted(int8): single element")
            if(allocated(error)) return
            output = unique(A, .false.)
            expected_stable = [1_int8]
            call check(error, all(output==expected_stable), .true.,&
                "Stable(int8): single element")
            if(allocated(error)) return

            A = [4_int8, 4_int8, 4_int8, 4_int8]
            output = unique(A, .true.)
            expected_sorted = [4_int8]
            call check(error, all(output==expected_sorted), .true.,&
                "Sorted(int8): all duplicates")
            if(allocated(error)) return
            output = unique(A, .false.)
            expected_stable = [4_int8]
            call check(error, all(output==expected_stable), .true.,&
                "Stable(int8): all duplicates")
            if(allocated(error)) return

            A = [1_int8, 2_int8, 3_int8, 4_int8]
            output = unique(A, .true.)
            expected_sorted = [1_int8, 2_int8, 3_int8, 4_int8]
            call check(error, all(output==expected_sorted), .true.,&
                "Sorted(int8): already unique")
            if(allocated(error)) return
            output = unique(A, .false.)
            expected_stable = [1_int8, 2_int8, 3_int8, 4_int8]
            call check(error, all(output==expected_stable), .true.,&
                "Stable(int8): already unique")
            if(allocated(error)) return

            A = [5_int8, 2_int8, 3_int8, 5_int8, 2_int8, 1_int8]
            output = unique(A, .true.)
            expected_sorted = [1_int8, 2_int8, 3_int8, 5_int8]
            call check(error, all(output==expected_sorted), .true.,&
                "Sorted(int8): mixed duplicates")
            if(allocated(error)) return
            output = unique(A, .false.)
            expected_stable = [5_int8, 2_int8, 3_int8, 1_int8]
            call check(error, all(output==expected_stable), .true.,&
                "Stable(int8): mixed duplicates")
            if(allocated(error)) return

            A = [-2_int8, 5_int8, -2_int8, 1_int8, 5_int8]
            output = unique(A, .true.)
            expected_sorted = [-2_int8, 1_int8, 5_int8]
            call check(error, all(output==expected_sorted), .true.,&
                "Sorted(int8): negatives and duplicates")
            if(allocated(error)) return
            output = unique(A, .false.)
            expected_stable = [-2_int8, 5_int8, 1_int8]
            call check(error, all(output==expected_stable), .true.,&
                "Stable(int8): negatives and duplicates")
            if(allocated(error)) return
        end block
        block
            integer(int16), allocatable :: A(:)
            integer(int16), allocatable :: output(:)
            integer(int16), allocatable :: expected_sorted(:)
            integer(int16), allocatable :: expected_stable(:)

            ! Initialize matrix.
            allocate(A(0))
            output = unique(A, .true.)
            allocate(expected_sorted(0))
            call check(error, all(output==expected_sorted), .true.,&
                "Sorted(int16): empty array")
            if(allocated(error)) return
            output = unique(A, .false.)
            allocate(expected_stable(0))
            call check(error, all(output==expected_stable), .true.,&
                "Stable(int16): empty array")
            if(allocated(error)) return

            A = [1_int16]
            output = unique(A, .true.)
            expected_sorted = [1_int16]
            call check(error, all(output==expected_sorted), .true.,&
                "Sorted(int16): single element")
            if(allocated(error)) return
            output = unique(A, .false.)
            expected_stable = [1_int16]
            call check(error, all(output==expected_stable), .true.,&
                "Stable(int16): single element")
            if(allocated(error)) return

            A = [4_int16, 4_int16, 4_int16, 4_int16]
            output = unique(A, .true.)
            expected_sorted = [4_int16]
            call check(error, all(output==expected_sorted), .true.,&
                "Sorted(int16): all duplicates")
            if(allocated(error)) return
            output = unique(A, .false.)
            expected_stable = [4_int16]
            call check(error, all(output==expected_stable), .true.,&
                "Stable(int16): all duplicates")
            if(allocated(error)) return

            A = [1_int16, 2_int16, 3_int16, 4_int16]
            output = unique(A, .true.)
            expected_sorted = [1_int16, 2_int16, 3_int16, 4_int16]
            call check(error, all(output==expected_sorted), .true.,&
                "Sorted(int16): already unique")
            if(allocated(error)) return
            output = unique(A, .false.)
            expected_stable = [1_int16, 2_int16, 3_int16, 4_int16]
            call check(error, all(output==expected_stable), .true.,&
                "Stable(int16): already unique")
            if(allocated(error)) return

            A = [5_int16, 2_int16, 3_int16, 5_int16, 2_int16, 1_int16]
            output = unique(A, .true.)
            expected_sorted = [1_int16, 2_int16, 3_int16, 5_int16]
            call check(error, all(output==expected_sorted), .true.,&
                "Sorted(int16): mixed duplicates")
            if(allocated(error)) return
            output = unique(A, .false.)
            expected_stable = [5_int16, 2_int16, 3_int16, 1_int16]
            call check(error, all(output==expected_stable), .true.,&
                "Stable(int16): mixed duplicates")
            if(allocated(error)) return

            A = [-2_int16, 5_int16, -2_int16, 1_int16, 5_int16]
            output = unique(A, .true.)
            expected_sorted = [-2_int16, 1_int16, 5_int16]
            call check(error, all(output==expected_sorted), .true.,&
                "Sorted(int16): negatives and duplicates")
            if(allocated(error)) return
            output = unique(A, .false.)
            expected_stable = [-2_int16, 5_int16, 1_int16]
            call check(error, all(output==expected_stable), .true.,&
                "Stable(int16): negatives and duplicates")
            if(allocated(error)) return
        end block
        block
            integer(int32), allocatable :: A(:)
            integer(int32), allocatable :: output(:)
            integer(int32), allocatable :: expected_sorted(:)
            integer(int32), allocatable :: expected_stable(:)

            ! Initialize matrix.
            allocate(A(0))
            output = unique(A, .true.)
            allocate(expected_sorted(0))
            call check(error, all(output==expected_sorted), .true.,&
                "Sorted(int32): empty array")
            if(allocated(error)) return
            output = unique(A, .false.)
            allocate(expected_stable(0))
            call check(error, all(output==expected_stable), .true.,&
                "Stable(int32): empty array")
            if(allocated(error)) return

            A = [1_int32]
            output = unique(A, .true.)
            expected_sorted = [1_int32]
            call check(error, all(output==expected_sorted), .true.,&
                "Sorted(int32): single element")
            if(allocated(error)) return
            output = unique(A, .false.)
            expected_stable = [1_int32]
            call check(error, all(output==expected_stable), .true.,&
                "Stable(int32): single element")
            if(allocated(error)) return

            A = [4_int32, 4_int32, 4_int32, 4_int32]
            output = unique(A, .true.)
            expected_sorted = [4_int32]
            call check(error, all(output==expected_sorted), .true.,&
                "Sorted(int32): all duplicates")
            if(allocated(error)) return
            output = unique(A, .false.)
            expected_stable = [4_int32]
            call check(error, all(output==expected_stable), .true.,&
                "Stable(int32): all duplicates")
            if(allocated(error)) return

            A = [1_int32, 2_int32, 3_int32, 4_int32]
            output = unique(A, .true.)
            expected_sorted = [1_int32, 2_int32, 3_int32, 4_int32]
            call check(error, all(output==expected_sorted), .true.,&
                "Sorted(int32): already unique")
            if(allocated(error)) return
            output = unique(A, .false.)
            expected_stable = [1_int32, 2_int32, 3_int32, 4_int32]
            call check(error, all(output==expected_stable), .true.,&
                "Stable(int32): already unique")
            if(allocated(error)) return

            A = [5_int32, 2_int32, 3_int32, 5_int32, 2_int32, 1_int32]
            output = unique(A, .true.)
            expected_sorted = [1_int32, 2_int32, 3_int32, 5_int32]
            call check(error, all(output==expected_sorted), .true.,&
                "Sorted(int32): mixed duplicates")
            if(allocated(error)) return
            output = unique(A, .false.)
            expected_stable = [5_int32, 2_int32, 3_int32, 1_int32]
            call check(error, all(output==expected_stable), .true.,&
                "Stable(int32): mixed duplicates")
            if(allocated(error)) return

            A = [-2_int32, 5_int32, -2_int32, 1_int32, 5_int32]
            output = unique(A, .true.)
            expected_sorted = [-2_int32, 1_int32, 5_int32]
            call check(error, all(output==expected_sorted), .true.,&
                "Sorted(int32): negatives and duplicates")
            if(allocated(error)) return
            output = unique(A, .false.)
            expected_stable = [-2_int32, 5_int32, 1_int32]
            call check(error, all(output==expected_stable), .true.,&
                "Stable(int32): negatives and duplicates")
            if(allocated(error)) return
        end block
        block
            integer(int64), allocatable :: A(:)
            integer(int64), allocatable :: output(:)
            integer(int64), allocatable :: expected_sorted(:)
            integer(int64), allocatable :: expected_stable(:)

            ! Initialize matrix.
            allocate(A(0))
            output = unique(A, .true.)
            allocate(expected_sorted(0))
            call check(error, all(output==expected_sorted), .true.,&
                "Sorted(int64): empty array")
            if(allocated(error)) return
            output = unique(A, .false.)
            allocate(expected_stable(0))
            call check(error, all(output==expected_stable), .true.,&
                "Stable(int64): empty array")
            if(allocated(error)) return

            A = [1_int64]
            output = unique(A, .true.)
            expected_sorted = [1_int64]
            call check(error, all(output==expected_sorted), .true.,&
                "Sorted(int64): single element")
            if(allocated(error)) return
            output = unique(A, .false.)
            expected_stable = [1_int64]
            call check(error, all(output==expected_stable), .true.,&
                "Stable(int64): single element")
            if(allocated(error)) return

            A = [4_int64, 4_int64, 4_int64, 4_int64]
            output = unique(A, .true.)
            expected_sorted = [4_int64]
            call check(error, all(output==expected_sorted), .true.,&
                "Sorted(int64): all duplicates")
            if(allocated(error)) return
            output = unique(A, .false.)
            expected_stable = [4_int64]
            call check(error, all(output==expected_stable), .true.,&
                "Stable(int64): all duplicates")
            if(allocated(error)) return

            A = [1_int64, 2_int64, 3_int64, 4_int64]
            output = unique(A, .true.)
            expected_sorted = [1_int64, 2_int64, 3_int64, 4_int64]
            call check(error, all(output==expected_sorted), .true.,&
                "Sorted(int64): already unique")
            if(allocated(error)) return
            output = unique(A, .false.)
            expected_stable = [1_int64, 2_int64, 3_int64, 4_int64]
            call check(error, all(output==expected_stable), .true.,&
                "Stable(int64): already unique")
            if(allocated(error)) return

            A = [5_int64, 2_int64, 3_int64, 5_int64, 2_int64, 1_int64]
            output = unique(A, .true.)
            expected_sorted = [1_int64, 2_int64, 3_int64, 5_int64]
            call check(error, all(output==expected_sorted), .true.,&
                "Sorted(int64): mixed duplicates")
            if(allocated(error)) return
            output = unique(A, .false.)
            expected_stable = [5_int64, 2_int64, 3_int64, 1_int64]
            call check(error, all(output==expected_stable), .true.,&
                "Stable(int64): mixed duplicates")
            if(allocated(error)) return

            A = [-2_int64, 5_int64, -2_int64, 1_int64, 5_int64]
            output = unique(A, .true.)
            expected_sorted = [-2_int64, 1_int64, 5_int64]
            call check(error, all(output==expected_sorted), .true.,&
                "Sorted(int64): negatives and duplicates")
            if(allocated(error)) return
            output = unique(A, .false.)
            expected_stable = [-2_int64, 5_int64, 1_int64]
            call check(error, all(output==expected_stable), .true.,&
                "Stable(int64): negatives and duplicates")
            if(allocated(error)) return
        end block
    end subroutine

    subroutine test_real(error)
        !> Error handling
        type(error_type), allocatable, intent(out) :: error
        block
            real(sp), allocatable :: A(:)
            real(sp), allocatable :: output(:)
            real(sp), allocatable :: expected_sorted(:)
            real(sp), allocatable :: expected_stable(:)
            integer(int8), allocatable :: bytes(:)

            ! Initialize matrix.
            allocate(A(0))
            output = unique(A, .true.)
            allocate(expected_sorted(0))
            call check(error, all_close(output, expected_sorted), .true.,&
                "Sorted(sp): empty array")
            if(allocated(error)) return
            output = unique(A, .false.)
            allocate(expected_stable(0))
            call check(error, all_close(output, expected_stable), .true.,&
                "Stable(sp): empty array")
            if(allocated(error)) return

            A = [1.0_sp]
            output = unique(A, .true.)
            expected_sorted = [1.0_sp]
            call check(error, all_close(output, expected_sorted), .true.,&
                "Sorted(sp): single element")
            if(allocated(error)) return
            output = unique(A, .false.)
            expected_stable = [1.0_sp]
            call check(error, all_close(output, expected_stable), .true.,&
                "Stable(sp): single element")
            if(allocated(error)) return

            A = [4.0_sp, 4.0_sp, 4.0_sp, 4.0_sp]
            output = unique(A, .true.)
            expected_sorted = [4.0_sp]
            call check(error, all_close(output, expected_sorted), .true.,&
                "Sorted(sp): all duplicates")
            if(allocated(error)) return
            output = unique(A, .false.)
            expected_stable = [4.0_sp]
            call check(error, all_close(output, expected_stable), .true.,&
                "Stable(sp): all duplicates")
            if(allocated(error)) return

            A = [1.0_sp, 2.0_sp, 3.0_sp, 4.0_sp]
            output = unique(A, .true.)
            expected_sorted = [1.0_sp, 2.0_sp, 3.0_sp, 4.0_sp]
            call check(error, all_close(output, expected_sorted), .true.,&
                "Sorted(sp): already unique")
            if(allocated(error)) return
            output = unique(A, .false.)
            expected_stable = [1.0_sp, 2.0_sp, 3.0_sp, 4.0_sp]
            call check(error, all_close(output, expected_stable), .true.,&
                "Stable(sp): already unique")
            if(allocated(error)) return

            A = [5.0_sp, 2.0_sp, 3.0_sp, 5.0_sp, 2.0_sp, 1.0_sp]
            output = unique(A, .true.)
            expected_sorted = [1.0_sp, 2.0_sp, 3.0_sp, 5.0_sp]
            call check(error, all_close(output, expected_sorted), .true.,&
                "Sorted(sp): mixed duplicates")
            if(allocated(error)) return
            output = unique(A, .false.)
            expected_stable = [5.0_sp, 2.0_sp, 3.0_sp, 1.0_sp]
            call check(error, all_close(output, expected_stable), .true.,&
                "Stable(sp): mixed duplicates")
            if(allocated(error)) return

            A = [-2.0_sp, 5.0_sp, -2.0_sp, 1.0_sp, 5.0_sp]
            output = unique(A, .true.)
            expected_sorted = [-2.0_sp, 1.0_sp, 5.0_sp]
            call check(error, all_close(output, expected_sorted), .true.,&
                "Sorted(sp): negatives and duplicates")
            if(allocated(error)) return
            output = unique(A, .false.)
            expected_stable = [-2.0_sp, 5.0_sp, 1.0_sp]
            call check(error, all_close(output, expected_stable), .true.,&
                "Stable(sp): negatives and duplicates")
            if(allocated(error)) return
        end block
        block
            real(dp), allocatable :: A(:)
            real(dp), allocatable :: output(:)
            real(dp), allocatable :: expected_sorted(:)
            real(dp), allocatable :: expected_stable(:)
            integer(int8), allocatable :: bytes(:)

            ! Initialize matrix.
            allocate(A(0))
            output = unique(A, .true.)
            allocate(expected_sorted(0))
            call check(error, all_close(output, expected_sorted), .true.,&
                "Sorted(dp): empty array")
            if(allocated(error)) return
            output = unique(A, .false.)
            allocate(expected_stable(0))
            call check(error, all_close(output, expected_stable), .true.,&
                "Stable(dp): empty array")
            if(allocated(error)) return

            A = [1.0_dp]
            output = unique(A, .true.)
            expected_sorted = [1.0_dp]
            call check(error, all_close(output, expected_sorted), .true.,&
                "Sorted(dp): single element")
            if(allocated(error)) return
            output = unique(A, .false.)
            expected_stable = [1.0_dp]
            call check(error, all_close(output, expected_stable), .true.,&
                "Stable(dp): single element")
            if(allocated(error)) return

            A = [4.0_dp, 4.0_dp, 4.0_dp, 4.0_dp]
            output = unique(A, .true.)
            expected_sorted = [4.0_dp]
            call check(error, all_close(output, expected_sorted), .true.,&
                "Sorted(dp): all duplicates")
            if(allocated(error)) return
            output = unique(A, .false.)
            expected_stable = [4.0_dp]
            call check(error, all_close(output, expected_stable), .true.,&
                "Stable(dp): all duplicates")
            if(allocated(error)) return

            A = [1.0_dp, 2.0_dp, 3.0_dp, 4.0_dp]
            output = unique(A, .true.)
            expected_sorted = [1.0_dp, 2.0_dp, 3.0_dp, 4.0_dp]
            call check(error, all_close(output, expected_sorted), .true.,&
                "Sorted(dp): already unique")
            if(allocated(error)) return
            output = unique(A, .false.)
            expected_stable = [1.0_dp, 2.0_dp, 3.0_dp, 4.0_dp]
            call check(error, all_close(output, expected_stable), .true.,&
                "Stable(dp): already unique")
            if(allocated(error)) return

            A = [5.0_dp, 2.0_dp, 3.0_dp, 5.0_dp, 2.0_dp, 1.0_dp]
            output = unique(A, .true.)
            expected_sorted = [1.0_dp, 2.0_dp, 3.0_dp, 5.0_dp]
            call check(error, all_close(output, expected_sorted), .true.,&
                "Sorted(dp): mixed duplicates")
            if(allocated(error)) return
            output = unique(A, .false.)
            expected_stable = [5.0_dp, 2.0_dp, 3.0_dp, 1.0_dp]
            call check(error, all_close(output, expected_stable), .true.,&
                "Stable(dp): mixed duplicates")
            if(allocated(error)) return

            A = [-2.0_dp, 5.0_dp, -2.0_dp, 1.0_dp, 5.0_dp]
            output = unique(A, .true.)
            expected_sorted = [-2.0_dp, 1.0_dp, 5.0_dp]
            call check(error, all_close(output, expected_sorted), .true.,&
                "Sorted(dp): negatives and duplicates")
            if(allocated(error)) return
            output = unique(A, .false.)
            expected_stable = [-2.0_dp, 5.0_dp, 1.0_dp]
            call check(error, all_close(output, expected_stable), .true.,&
                "Stable(dp): negatives and duplicates")
            if(allocated(error)) return
        end block
    end subroutine

end module

program tester
    use, intrinsic :: iso_fortran_env, only : error_unit
    use testdrive, only : run_testsuite, new_testsuite, testsuite_type
    use test_sorting_unique, only : collect_suite
    implicit none
    integer :: stat, is
    type(testsuite_type), allocatable :: testsuites(:)
    character(len=*), parameter :: fmt = '("#", *(1x, a))'

    stat = 0

    testsuites = [ &
        new_testsuite("sorting_unique", collect_suite) &
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