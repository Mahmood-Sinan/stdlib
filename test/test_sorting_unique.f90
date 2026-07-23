#include "macros.inc"
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
            integer(int8), allocatable :: expected_unsorted(:)

            ! Initialize matrix.
            allocate(A(0))
            output = unique(A, .true.)
            allocate(expected_sorted(0))
            call check(error, all(output==expected_sorted), .true.,&
                "Sorted(int8): empty array")
            if(allocated(error)) return
#if STDLIB_HASHMAPS
            output = unique(A, .false.)
            allocate(expected_unsorted(0))
            call check(error, all(output==expected_unsorted), .true.,&
                "Unsorted(int8): empty array")
            if(allocated(error)) return
#endif

            A = [1_int8]
            output = unique(A, .true.)
            expected_sorted = [1_int8]
            call check(error, all(output==expected_sorted), .true.,&
                "Sorted(int8): single element")
            if(allocated(error)) return
#if STDLIB_HASHMAPS
            output = unique(A, .false.)
            expected_unsorted = [1_int8]
            call check(error, all(output==expected_unsorted), .true.,&
                "Unsorted(int8): single element")
            if(allocated(error)) return
#endif

            A = [4_int8, 4_int8, 4_int8, 4_int8]
            output = unique(A, .true.)
            expected_sorted = [4_int8]
            call check(error, all(output==expected_sorted), .true.,&
                "Sorted(int8): all duplicates")
            if(allocated(error)) return
#if STDLIB_HASHMAPS
            output = unique(A, .false.)
            expected_unsorted = [4_int8]
            call check(error, all(output==expected_unsorted), .true.,&
                "Unsorted(int8): all duplicates")
            if(allocated(error)) return
#endif

            A = [1_int8, 2_int8, 3_int8, 4_int8]
            output = unique(A, .true.)
            expected_sorted = [1_int8, 2_int8, 3_int8, 4_int8]
            call check(error, all(output==expected_sorted), .true.,&
                "Sorted(int8): already unique")
            if(allocated(error)) return
#if STDLIB_HASHMAPS
            output = unique(A, .false.)
            expected_unsorted = [1_int8, 2_int8, 3_int8, 4_int8]
            call check(error, all(output==expected_unsorted), .true.,&
                "Unsorted(int8): already unique")
            if(allocated(error)) return
#endif

            A = [5_int8, 2_int8, 3_int8, 5_int8, 2_int8, 1_int8]
            output = unique(A, .true.)
            expected_sorted = [1_int8, 2_int8, 3_int8, 5_int8]
            call check(error, all(output==expected_sorted), .true.,&
                "Sorted(int8): mixed duplicates")
            if(allocated(error)) return
#if STDLIB_HASHMAPS
            output = unique(A, .false.)
            expected_unsorted = [5_int8, 2_int8, 3_int8, 1_int8]
            call check(error, all(output==expected_unsorted), .true.,&
                "Unsorted(int8): mixed duplicates")
            if(allocated(error)) return
#endif

            A = [-2_int8, 5_int8, -2_int8, 1_int8, 5_int8]
            output = unique(A, .true.)
            expected_sorted = [-2_int8, 1_int8, 5_int8]
            call check(error, all(output==expected_sorted), .true.,&
                "Sorted(int8): negatives and duplicates")
            if(allocated(error)) return
#if STDLIB_HASHMAPS
            output = unique(A, .false.)
            expected_unsorted = [-2_int8, 5_int8, 1_int8]
            call check(error, all(output==expected_unsorted), .true.,&
                "Unsorted(int8): negatives and duplicates")
            if(allocated(error)) return
#endif
        end block
        block
            integer(int16), allocatable :: A(:)
            integer(int16), allocatable :: output(:)
            integer(int16), allocatable :: expected_sorted(:)
            integer(int16), allocatable :: expected_unsorted(:)

            ! Initialize matrix.
            allocate(A(0))
            output = unique(A, .true.)
            allocate(expected_sorted(0))
            call check(error, all(output==expected_sorted), .true.,&
                "Sorted(int16): empty array")
            if(allocated(error)) return
#if STDLIB_HASHMAPS
            output = unique(A, .false.)
            allocate(expected_unsorted(0))
            call check(error, all(output==expected_unsorted), .true.,&
                "Unsorted(int16): empty array")
            if(allocated(error)) return
#endif

            A = [1_int16]
            output = unique(A, .true.)
            expected_sorted = [1_int16]
            call check(error, all(output==expected_sorted), .true.,&
                "Sorted(int16): single element")
            if(allocated(error)) return
#if STDLIB_HASHMAPS
            output = unique(A, .false.)
            expected_unsorted = [1_int16]
            call check(error, all(output==expected_unsorted), .true.,&
                "Unsorted(int16): single element")
            if(allocated(error)) return
#endif

            A = [4_int16, 4_int16, 4_int16, 4_int16]
            output = unique(A, .true.)
            expected_sorted = [4_int16]
            call check(error, all(output==expected_sorted), .true.,&
                "Sorted(int16): all duplicates")
            if(allocated(error)) return
#if STDLIB_HASHMAPS
            output = unique(A, .false.)
            expected_unsorted = [4_int16]
            call check(error, all(output==expected_unsorted), .true.,&
                "Unsorted(int16): all duplicates")
            if(allocated(error)) return
#endif

            A = [1_int16, 2_int16, 3_int16, 4_int16]
            output = unique(A, .true.)
            expected_sorted = [1_int16, 2_int16, 3_int16, 4_int16]
            call check(error, all(output==expected_sorted), .true.,&
                "Sorted(int16): already unique")
            if(allocated(error)) return
#if STDLIB_HASHMAPS
            output = unique(A, .false.)
            expected_unsorted = [1_int16, 2_int16, 3_int16, 4_int16]
            call check(error, all(output==expected_unsorted), .true.,&
                "Unsorted(int16): already unique")
            if(allocated(error)) return
#endif

            A = [5_int16, 2_int16, 3_int16, 5_int16, 2_int16, 1_int16]
            output = unique(A, .true.)
            expected_sorted = [1_int16, 2_int16, 3_int16, 5_int16]
            call check(error, all(output==expected_sorted), .true.,&
                "Sorted(int16): mixed duplicates")
            if(allocated(error)) return
#if STDLIB_HASHMAPS
            output = unique(A, .false.)
            expected_unsorted = [5_int16, 2_int16, 3_int16, 1_int16]
            call check(error, all(output==expected_unsorted), .true.,&
                "Unsorted(int16): mixed duplicates")
            if(allocated(error)) return
#endif

            A = [-2_int16, 5_int16, -2_int16, 1_int16, 5_int16]
            output = unique(A, .true.)
            expected_sorted = [-2_int16, 1_int16, 5_int16]
            call check(error, all(output==expected_sorted), .true.,&
                "Sorted(int16): negatives and duplicates")
            if(allocated(error)) return
#if STDLIB_HASHMAPS
            output = unique(A, .false.)
            expected_unsorted = [-2_int16, 5_int16, 1_int16]
            call check(error, all(output==expected_unsorted), .true.,&
                "Unsorted(int16): negatives and duplicates")
            if(allocated(error)) return
#endif
        end block
        block
            integer(int32), allocatable :: A(:)
            integer(int32), allocatable :: output(:)
            integer(int32), allocatable :: expected_sorted(:)
            integer(int32), allocatable :: expected_unsorted(:)

            ! Initialize matrix.
            allocate(A(0))
            output = unique(A, .true.)
            allocate(expected_sorted(0))
            call check(error, all(output==expected_sorted), .true.,&
                "Sorted(int32): empty array")
            if(allocated(error)) return
#if STDLIB_HASHMAPS
            output = unique(A, .false.)
            allocate(expected_unsorted(0))
            call check(error, all(output==expected_unsorted), .true.,&
                "Unsorted(int32): empty array")
            if(allocated(error)) return
#endif

            A = [1_int32]
            output = unique(A, .true.)
            expected_sorted = [1_int32]
            call check(error, all(output==expected_sorted), .true.,&
                "Sorted(int32): single element")
            if(allocated(error)) return
#if STDLIB_HASHMAPS
            output = unique(A, .false.)
            expected_unsorted = [1_int32]
            call check(error, all(output==expected_unsorted), .true.,&
                "Unsorted(int32): single element")
            if(allocated(error)) return
#endif

            A = [4_int32, 4_int32, 4_int32, 4_int32]
            output = unique(A, .true.)
            expected_sorted = [4_int32]
            call check(error, all(output==expected_sorted), .true.,&
                "Sorted(int32): all duplicates")
            if(allocated(error)) return
#if STDLIB_HASHMAPS
            output = unique(A, .false.)
            expected_unsorted = [4_int32]
            call check(error, all(output==expected_unsorted), .true.,&
                "Unsorted(int32): all duplicates")
            if(allocated(error)) return
#endif

            A = [1_int32, 2_int32, 3_int32, 4_int32]
            output = unique(A, .true.)
            expected_sorted = [1_int32, 2_int32, 3_int32, 4_int32]
            call check(error, all(output==expected_sorted), .true.,&
                "Sorted(int32): already unique")
            if(allocated(error)) return
#if STDLIB_HASHMAPS
            output = unique(A, .false.)
            expected_unsorted = [1_int32, 2_int32, 3_int32, 4_int32]
            call check(error, all(output==expected_unsorted), .true.,&
                "Unsorted(int32): already unique")
            if(allocated(error)) return
#endif

            A = [5_int32, 2_int32, 3_int32, 5_int32, 2_int32, 1_int32]
            output = unique(A, .true.)
            expected_sorted = [1_int32, 2_int32, 3_int32, 5_int32]
            call check(error, all(output==expected_sorted), .true.,&
                "Sorted(int32): mixed duplicates")
            if(allocated(error)) return
#if STDLIB_HASHMAPS
            output = unique(A, .false.)
            expected_unsorted = [5_int32, 2_int32, 3_int32, 1_int32]
            call check(error, all(output==expected_unsorted), .true.,&
                "Unsorted(int32): mixed duplicates")
            if(allocated(error)) return
#endif

            A = [-2_int32, 5_int32, -2_int32, 1_int32, 5_int32]
            output = unique(A, .true.)
            expected_sorted = [-2_int32, 1_int32, 5_int32]
            call check(error, all(output==expected_sorted), .true.,&
                "Sorted(int32): negatives and duplicates")
            if(allocated(error)) return
#if STDLIB_HASHMAPS
            output = unique(A, .false.)
            expected_unsorted = [-2_int32, 5_int32, 1_int32]
            call check(error, all(output==expected_unsorted), .true.,&
                "Unsorted(int32): negatives and duplicates")
            if(allocated(error)) return
#endif
        end block
        block
            integer(int64), allocatable :: A(:)
            integer(int64), allocatable :: output(:)
            integer(int64), allocatable :: expected_sorted(:)
            integer(int64), allocatable :: expected_unsorted(:)

            ! Initialize matrix.
            allocate(A(0))
            output = unique(A, .true.)
            allocate(expected_sorted(0))
            call check(error, all(output==expected_sorted), .true.,&
                "Sorted(int64): empty array")
            if(allocated(error)) return
#if STDLIB_HASHMAPS
            output = unique(A, .false.)
            allocate(expected_unsorted(0))
            call check(error, all(output==expected_unsorted), .true.,&
                "Unsorted(int64): empty array")
            if(allocated(error)) return
#endif

            A = [1_int64]
            output = unique(A, .true.)
            expected_sorted = [1_int64]
            call check(error, all(output==expected_sorted), .true.,&
                "Sorted(int64): single element")
            if(allocated(error)) return
#if STDLIB_HASHMAPS
            output = unique(A, .false.)
            expected_unsorted = [1_int64]
            call check(error, all(output==expected_unsorted), .true.,&
                "Unsorted(int64): single element")
            if(allocated(error)) return
#endif

            A = [4_int64, 4_int64, 4_int64, 4_int64]
            output = unique(A, .true.)
            expected_sorted = [4_int64]
            call check(error, all(output==expected_sorted), .true.,&
                "Sorted(int64): all duplicates")
            if(allocated(error)) return
#if STDLIB_HASHMAPS
            output = unique(A, .false.)
            expected_unsorted = [4_int64]
            call check(error, all(output==expected_unsorted), .true.,&
                "Unsorted(int64): all duplicates")
            if(allocated(error)) return
#endif

            A = [1_int64, 2_int64, 3_int64, 4_int64]
            output = unique(A, .true.)
            expected_sorted = [1_int64, 2_int64, 3_int64, 4_int64]
            call check(error, all(output==expected_sorted), .true.,&
                "Sorted(int64): already unique")
            if(allocated(error)) return
#if STDLIB_HASHMAPS
            output = unique(A, .false.)
            expected_unsorted = [1_int64, 2_int64, 3_int64, 4_int64]
            call check(error, all(output==expected_unsorted), .true.,&
                "Unsorted(int64): already unique")
            if(allocated(error)) return
#endif

            A = [5_int64, 2_int64, 3_int64, 5_int64, 2_int64, 1_int64]
            output = unique(A, .true.)
            expected_sorted = [1_int64, 2_int64, 3_int64, 5_int64]
            call check(error, all(output==expected_sorted), .true.,&
                "Sorted(int64): mixed duplicates")
            if(allocated(error)) return
#if STDLIB_HASHMAPS
            output = unique(A, .false.)
            expected_unsorted = [5_int64, 2_int64, 3_int64, 1_int64]
            call check(error, all(output==expected_unsorted), .true.,&
                "Unsorted(int64): mixed duplicates")
            if(allocated(error)) return
#endif

            A = [-2_int64, 5_int64, -2_int64, 1_int64, 5_int64]
            output = unique(A, .true.)
            expected_sorted = [-2_int64, 1_int64, 5_int64]
            call check(error, all(output==expected_sorted), .true.,&
                "Sorted(int64): negatives and duplicates")
            if(allocated(error)) return
#if STDLIB_HASHMAPS
            output = unique(A, .false.)
            expected_unsorted = [-2_int64, 5_int64, 1_int64]
            call check(error, all(output==expected_unsorted), .true.,&
                "Unsorted(int64): negatives and duplicates")
            if(allocated(error)) return
#endif
        end block
    end subroutine

    subroutine test_real(error)
        !> Error handling
        type(error_type), allocatable, intent(out) :: error
        block
            real(sp), allocatable :: A(:)
            real(sp), allocatable :: output(:)
            real(sp), allocatable :: expected_sorted(:)
            real(sp), allocatable :: expected_unsorted(:)
            integer(int8), allocatable :: bytes(:)

            ! Initialize matrix.
            allocate(A(0))
            output = unique(A, .true.)
            allocate(expected_sorted(0))
            call check(error, all_close(output, expected_sorted), .true.,&
                "Sorted(sp): empty array")
            if(allocated(error)) return
#if STDLIB_HASHMAPS
            output = unique(A, .false.)
            allocate(expected_unsorted(0))
            call check(error, all_close(output, expected_unsorted), .true.,&
                "Unsorted(sp): empty array")
            if(allocated(error)) return
#endif

            A = [1.0_sp]
            output = unique(A, .true.)
            expected_sorted = [1.0_sp]
            call check(error, all_close(output, expected_sorted), .true.,&
                "Sorted(sp): single element")
            if(allocated(error)) return
#if STDLIB_HASHMAPS
            output = unique(A, .false.)
            expected_unsorted = [1.0_sp]
            call check(error, all_close(output, expected_unsorted), .true.,&
                "Unsorted(sp): single element")
            if(allocated(error)) return
#endif

            A = [4.0_sp, 4.0_sp, 4.0_sp, 4.0_sp]
            output = unique(A, .true.)
            expected_sorted = [4.0_sp]
            call check(error, all_close(output, expected_sorted), .true.,&
                "Sorted(sp): all duplicates")
            if(allocated(error)) return
#if STDLIB_HASHMAPS
            output = unique(A, .false.)
            expected_unsorted = [4.0_sp]
            call check(error, all_close(output, expected_unsorted), .true.,&
                "Unsorted(sp): all duplicates")
            if(allocated(error)) return
#endif

            A = [1.0_sp, 2.0_sp, 3.0_sp, 4.0_sp]
            output = unique(A, .true.)
            expected_sorted = [1.0_sp, 2.0_sp, 3.0_sp, 4.0_sp]
            call check(error, all_close(output, expected_sorted), .true.,&
                "Sorted(sp): already unique")
            if(allocated(error)) return
#if STDLIB_HASHMAPS
            output = unique(A, .false.)
            expected_unsorted = [1.0_sp, 2.0_sp, 3.0_sp, 4.0_sp]
            call check(error, all_close(output, expected_unsorted), .true.,&
                "Unsorted(sp): already unique")
            if(allocated(error)) return
#endif

            A = [5.0_sp, 2.0_sp, 3.0_sp, 5.0_sp, 2.0_sp, 1.0_sp]
            output = unique(A, .true.)
            expected_sorted = [1.0_sp, 2.0_sp, 3.0_sp, 5.0_sp]
            call check(error, all_close(output, expected_sorted), .true.,&
                "Sorted(sp): mixed duplicates")
            if(allocated(error)) return
#if STDLIB_HASHMAPS
            output = unique(A, .false.)
            expected_unsorted = [5.0_sp, 2.0_sp, 3.0_sp, 1.0_sp]
            call check(error, all_close(output, expected_unsorted), .true.,&
                "Unsorted(sp): mixed duplicates")
            if(allocated(error)) return
#endif

            A = [-2.0_sp, 5.0_sp, -2.0_sp, 1.0_sp, 5.0_sp]
            output = unique(A, .true.)
            expected_sorted = [-2.0_sp, 1.0_sp, 5.0_sp]
            call check(error, all_close(output, expected_sorted), .true.,&
                "Sorted(sp): negatives and duplicates")
            if(allocated(error)) return
#if STDLIB_HASHMAPS
            output = unique(A, .false.)
            expected_unsorted = [-2.0_sp, 5.0_sp, 1.0_sp]
            call check(error, all_close(output, expected_unsorted), .true.,&
                "Unsorted(sp): negatives and duplicates")
            if(allocated(error)) return
#endif
            ! Tolerance tests
            A = [3.00_sp, 1.05_sp, 1.00_sp, 2.05_sp, 2.00_sp]
            output = unique(A, .true., 0.1_sp)
            expected_sorted = [1.00_sp, 2.00_sp, 3.00_sp]
            call check(error, all_close(output, expected_sorted), .true., &
                "Sorted(sp): basic tolerance case")
            if (allocated(error)) return

            A = [1.00_sp, 2.00_sp, 1.00_sp, 1.0001_sp, 2.00_sp]
            output = unique(A, .true., 0.0_sp)
            expected_sorted = [1.00_sp, 1.0001_sp, 2.00_sp]
            call check(error, all_close(output, expected_sorted), .true., &
                "Sorted(sp): zero tolerance")
            if (allocated(error)) return

            A = [1.18_sp, 1.09_sp, 1.00_sp]
            output = unique(A, .true., 0.1_sp)
            expected_sorted = [1.0_sp, 1.18_sp]
            call check(error, all_close(output, expected_sorted), .true., &
                "Sorted(sp): representative tolerance")
            if (allocated(error)) return

            A = [5.0_sp, 1.0_sp, 3.0_sp]
            output = unique(A, .true., 10.0_sp)
            expected_sorted = [1.0_sp]
            call check(error, all_close(output, expected_sorted), .true., &
                "Sorted(sp): large tolerance")
            if (allocated(error)) return

            A = [1.0_sp, 1.05_sp, 1.005_sp]
            output = unique(A, .true., 0.01_sp)
            expected_sorted = [1.0_sp, 1.05_sp]
            call check(error, all_close(output, expected_sorted), .true., &
                "Sorted(sp): small tolerance")
            if (allocated(error)) return
        end block
        block
            real(dp), allocatable :: A(:)
            real(dp), allocatable :: output(:)
            real(dp), allocatable :: expected_sorted(:)
            real(dp), allocatable :: expected_unsorted(:)
            integer(int8), allocatable :: bytes(:)

            ! Initialize matrix.
            allocate(A(0))
            output = unique(A, .true.)
            allocate(expected_sorted(0))
            call check(error, all_close(output, expected_sorted), .true.,&
                "Sorted(dp): empty array")
            if(allocated(error)) return
#if STDLIB_HASHMAPS
            output = unique(A, .false.)
            allocate(expected_unsorted(0))
            call check(error, all_close(output, expected_unsorted), .true.,&
                "Unsorted(dp): empty array")
            if(allocated(error)) return
#endif

            A = [1.0_dp]
            output = unique(A, .true.)
            expected_sorted = [1.0_dp]
            call check(error, all_close(output, expected_sorted), .true.,&
                "Sorted(dp): single element")
            if(allocated(error)) return
#if STDLIB_HASHMAPS
            output = unique(A, .false.)
            expected_unsorted = [1.0_dp]
            call check(error, all_close(output, expected_unsorted), .true.,&
                "Unsorted(dp): single element")
            if(allocated(error)) return
#endif

            A = [4.0_dp, 4.0_dp, 4.0_dp, 4.0_dp]
            output = unique(A, .true.)
            expected_sorted = [4.0_dp]
            call check(error, all_close(output, expected_sorted), .true.,&
                "Sorted(dp): all duplicates")
            if(allocated(error)) return
#if STDLIB_HASHMAPS
            output = unique(A, .false.)
            expected_unsorted = [4.0_dp]
            call check(error, all_close(output, expected_unsorted), .true.,&
                "Unsorted(dp): all duplicates")
            if(allocated(error)) return
#endif

            A = [1.0_dp, 2.0_dp, 3.0_dp, 4.0_dp]
            output = unique(A, .true.)
            expected_sorted = [1.0_dp, 2.0_dp, 3.0_dp, 4.0_dp]
            call check(error, all_close(output, expected_sorted), .true.,&
                "Sorted(dp): already unique")
            if(allocated(error)) return
#if STDLIB_HASHMAPS
            output = unique(A, .false.)
            expected_unsorted = [1.0_dp, 2.0_dp, 3.0_dp, 4.0_dp]
            call check(error, all_close(output, expected_unsorted), .true.,&
                "Unsorted(dp): already unique")
            if(allocated(error)) return
#endif

            A = [5.0_dp, 2.0_dp, 3.0_dp, 5.0_dp, 2.0_dp, 1.0_dp]
            output = unique(A, .true.)
            expected_sorted = [1.0_dp, 2.0_dp, 3.0_dp, 5.0_dp]
            call check(error, all_close(output, expected_sorted), .true.,&
                "Sorted(dp): mixed duplicates")
            if(allocated(error)) return
#if STDLIB_HASHMAPS
            output = unique(A, .false.)
            expected_unsorted = [5.0_dp, 2.0_dp, 3.0_dp, 1.0_dp]
            call check(error, all_close(output, expected_unsorted), .true.,&
                "Unsorted(dp): mixed duplicates")
            if(allocated(error)) return
#endif

            A = [-2.0_dp, 5.0_dp, -2.0_dp, 1.0_dp, 5.0_dp]
            output = unique(A, .true.)
            expected_sorted = [-2.0_dp, 1.0_dp, 5.0_dp]
            call check(error, all_close(output, expected_sorted), .true.,&
                "Sorted(dp): negatives and duplicates")
            if(allocated(error)) return
#if STDLIB_HASHMAPS
            output = unique(A, .false.)
            expected_unsorted = [-2.0_dp, 5.0_dp, 1.0_dp]
            call check(error, all_close(output, expected_unsorted), .true.,&
                "Unsorted(dp): negatives and duplicates")
            if(allocated(error)) return
#endif
            ! Tolerance tests
            A = [3.00_dp, 1.05_dp, 1.00_dp, 2.05_dp, 2.00_dp]
            output = unique(A, .true., 0.1_dp)
            expected_sorted = [1.00_dp, 2.00_dp, 3.00_dp]
            call check(error, all_close(output, expected_sorted), .true., &
                "Sorted(dp): basic tolerance case")
            if (allocated(error)) return

            A = [1.00_dp, 2.00_dp, 1.00_dp, 1.0001_dp, 2.00_dp]
            output = unique(A, .true., 0.0_dp)
            expected_sorted = [1.00_dp, 1.0001_dp, 2.00_dp]
            call check(error, all_close(output, expected_sorted), .true., &
                "Sorted(dp): zero tolerance")
            if (allocated(error)) return

            A = [1.18_dp, 1.09_dp, 1.00_dp]
            output = unique(A, .true., 0.1_dp)
            expected_sorted = [1.0_dp, 1.18_dp]
            call check(error, all_close(output, expected_sorted), .true., &
                "Sorted(dp): representative tolerance")
            if (allocated(error)) return

            A = [5.0_dp, 1.0_dp, 3.0_dp]
            output = unique(A, .true., 10.0_dp)
            expected_sorted = [1.0_dp]
            call check(error, all_close(output, expected_sorted), .true., &
                "Sorted(dp): large tolerance")
            if (allocated(error)) return

            A = [1.0_dp, 1.05_dp, 1.005_dp]
            output = unique(A, .true., 0.01_dp)
            expected_sorted = [1.0_dp, 1.05_dp]
            call check(error, all_close(output, expected_sorted), .true., &
                "Sorted(dp): small tolerance")
            if (allocated(error)) return
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