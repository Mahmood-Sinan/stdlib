
submodule(stdlib_sorting) stdlib_sorting_unique
    use stdlib_hashmaps, only: chaining_hashmap_type
    use stdlib_hashmap_wrappers, only: key_type
    use stdlib_constants
    implicit none

contains
    module function int8_unique(A, sorted_output) result(output)
        integer(int8), intent(in) :: A(:)
        logical, intent(in) :: sorted_output
        integer(int8), allocatable :: output(:)

        integer(int8), allocatable:: temp(:)

        if(size(A) == 0) then
            allocate(output(0))
            return
        end if
        allocate(temp, source=A)
        if(sorted_output) then
            output = sort_unique(temp)
        else
            output = stable_unique(temp)
        end if
        deallocate(temp)
    end function
    module function int16_unique(A, sorted_output) result(output)
        integer(int16), intent(in) :: A(:)
        logical, intent(in) :: sorted_output
        integer(int16), allocatable :: output(:)

        integer(int16), allocatable:: temp(:)

        if(size(A) == 0) then
            allocate(output(0))
            return
        end if
        allocate(temp, source=A)
        if(sorted_output) then
            output = sort_unique(temp)
        else
            output = stable_unique(temp)
        end if
        deallocate(temp)
    end function
    module function int32_unique(A, sorted_output) result(output)
        integer(int32), intent(in) :: A(:)
        logical, intent(in) :: sorted_output
        integer(int32), allocatable :: output(:)

        integer(int32), allocatable:: temp(:)

        if(size(A) == 0) then
            allocate(output(0))
            return
        end if
        allocate(temp, source=A)
        if(sorted_output) then
            output = sort_unique(temp)
        else
            output = stable_unique(temp)
        end if
        deallocate(temp)
    end function
    module function int64_unique(A, sorted_output) result(output)
        integer(int64), intent(in) :: A(:)
        logical, intent(in) :: sorted_output
        integer(int64), allocatable :: output(:)

        integer(int64), allocatable:: temp(:)

        if(size(A) == 0) then
            allocate(output(0))
            return
        end if
        allocate(temp, source=A)
        if(sorted_output) then
            output = sort_unique(temp)
        else
            output = stable_unique(temp)
        end if
        deallocate(temp)
    end function
    module function sp_unique(A, sorted_output) result(output)
        real(sp), intent(in) :: A(:)
        logical, intent(in) :: sorted_output
        real(sp), allocatable :: output(:)

        real(sp), allocatable:: temp(:)

        if(size(A) == 0) then
            allocate(output(0))
            return
        end if
        allocate(temp, source=A)
        if(sorted_output) then
            output = sort_unique(temp)
        else
            output = stable_unique(temp)
        end if
        deallocate(temp)
    end function
    module function dp_unique(A, sorted_output) result(output)
        real(dp), intent(in) :: A(:)
        logical, intent(in) :: sorted_output
        real(dp), allocatable :: output(:)

        real(dp), allocatable:: temp(:)

        if(size(A) == 0) then
            allocate(output(0))
            return
        end if
        allocate(temp, source=A)
        if(sorted_output) then
            output = sort_unique(temp)
        else
            output = stable_unique(temp)
        end if
        deallocate(temp)
    end function

    module function int8_sort_unique(temp) result(output)
        integer(int8), intent(inout) :: temp(:)
        integer(int8), allocatable :: output(:)

        logical, allocatable :: mask(:)
        integer :: i

        allocate(mask(size(temp)))
        mask(1) = .true.
        call sort(temp)
        do i = 2, size(temp)
            mask(i) = temp(i) /= temp(i-1)
        end do
        output = pack(temp, mask)
        deallocate(mask)
    end function
    module function int16_sort_unique(temp) result(output)
        integer(int16), intent(inout) :: temp(:)
        integer(int16), allocatable :: output(:)

        logical, allocatable :: mask(:)
        integer :: i

        allocate(mask(size(temp)))
        mask(1) = .true.
        call sort(temp)
        do i = 2, size(temp)
            mask(i) = temp(i) /= temp(i-1)
        end do
        output = pack(temp, mask)
        deallocate(mask)
    end function
    module function int32_sort_unique(temp) result(output)
        integer(int32), intent(inout) :: temp(:)
        integer(int32), allocatable :: output(:)

        logical, allocatable :: mask(:)
        integer :: i

        allocate(mask(size(temp)))
        mask(1) = .true.
        call sort(temp)
        do i = 2, size(temp)
            mask(i) = temp(i) /= temp(i-1)
        end do
        output = pack(temp, mask)
        deallocate(mask)
    end function
    module function int64_sort_unique(temp) result(output)
        integer(int64), intent(inout) :: temp(:)
        integer(int64), allocatable :: output(:)

        logical, allocatable :: mask(:)
        integer :: i

        allocate(mask(size(temp)))
        mask(1) = .true.
        call sort(temp)
        do i = 2, size(temp)
            mask(i) = temp(i) /= temp(i-1)
        end do
        output = pack(temp, mask)
        deallocate(mask)
    end function
    module function sp_sort_unique(temp) result(output)
        real(sp), intent(inout) :: temp(:)
        real(sp), allocatable :: output(:)

        logical, allocatable :: mask(:)
        integer :: i

        allocate(mask(size(temp)))
        mask(1) = .true.
        call sort(temp)
        do i = 2, size(temp)
            mask(i) = temp(i) /= temp(i-1)
        end do
        output = pack(temp, mask)
        deallocate(mask)
    end function
    module function dp_sort_unique(temp) result(output)
        real(dp), intent(inout) :: temp(:)
        real(dp), allocatable :: output(:)

        logical, allocatable :: mask(:)
        integer :: i

        allocate(mask(size(temp)))
        mask(1) = .true.
        call sort(temp)
        do i = 2, size(temp)
            mask(i) = temp(i) /= temp(i-1)
        end do
        output = pack(temp, mask)
        deallocate(mask)
    end function

    module function int8_stable_unique(temp) result(output)
        integer(int8), intent(in) :: temp(:)
        integer(int8), allocatable :: output(:)

        type(chaining_hashmap_type) :: map
        logical, allocatable :: mask(:)
        logical :: present
        integer :: i
        integer(int8) :: key(storage_size(temp(1))/8)

        call map%init()
        allocate(mask(size(temp)))
        do i = 1, size(temp)
            key = 0
            key = transfer(temp(i), key)
            call map%key_test(key, present)
            if (.not. present) then
                call map%map_entry(key)
                mask(i) = .true.
            else
                mask(i) = .false.
            end if
        end do
        output = pack(temp, mask)
        deallocate(mask)
    end function
    module function int16_stable_unique(temp) result(output)
        integer(int16), intent(in) :: temp(:)
        integer(int16), allocatable :: output(:)

        type(chaining_hashmap_type) :: map
        logical, allocatable :: mask(:)
        logical :: present
        integer :: i
        integer(int8) :: key(storage_size(temp(1))/8)

        call map%init()
        allocate(mask(size(temp)))
        do i = 1, size(temp)
            key = 0
            key = transfer(temp(i), key)
            call map%key_test(key, present)
            if (.not. present) then
                call map%map_entry(key)
                mask(i) = .true.
            else
                mask(i) = .false.
            end if
        end do
        output = pack(temp, mask)
        deallocate(mask)
    end function
    module function int32_stable_unique(temp) result(output)
        integer(int32), intent(in) :: temp(:)
        integer(int32), allocatable :: output(:)

        type(chaining_hashmap_type) :: map
        logical, allocatable :: mask(:)
        logical :: present
        integer :: i
        integer(int8) :: key(storage_size(temp(1))/8)

        call map%init()
        allocate(mask(size(temp)))
        do i = 1, size(temp)
            key = 0
            key = transfer(temp(i), key)
            call map%key_test(key, present)
            if (.not. present) then
                call map%map_entry(key)
                mask(i) = .true.
            else
                mask(i) = .false.
            end if
        end do
        output = pack(temp, mask)
        deallocate(mask)
    end function
    module function int64_stable_unique(temp) result(output)
        integer(int64), intent(in) :: temp(:)
        integer(int64), allocatable :: output(:)

        type(chaining_hashmap_type) :: map
        logical, allocatable :: mask(:)
        logical :: present
        integer :: i
        integer(int8) :: key(storage_size(temp(1))/8)

        call map%init()
        allocate(mask(size(temp)))
        do i = 1, size(temp)
            key = 0
            key = transfer(temp(i), key)
            call map%key_test(key, present)
            if (.not. present) then
                call map%map_entry(key)
                mask(i) = .true.
            else
                mask(i) = .false.
            end if
        end do
        output = pack(temp, mask)
        deallocate(mask)
    end function
    module function sp_stable_unique(temp) result(output)
        real(sp), intent(in) :: temp(:)
        real(sp), allocatable :: output(:)

        type(chaining_hashmap_type) :: map
        logical, allocatable :: mask(:)
        logical :: present
        integer :: i
        integer(int8) :: key(storage_size(temp(1))/8)

        call map%init()
        allocate(mask(size(temp)))
        do i = 1, size(temp)
            key = 0
            key = transfer(temp(i), key)
            call map%key_test(key, present)
            if (.not. present) then
                call map%map_entry(key)
                mask(i) = .true.
            else
                mask(i) = .false.
            end if
        end do
        output = pack(temp, mask)
        deallocate(mask)
    end function
    module function dp_stable_unique(temp) result(output)
        real(dp), intent(in) :: temp(:)
        real(dp), allocatable :: output(:)

        type(chaining_hashmap_type) :: map
        logical, allocatable :: mask(:)
        logical :: present
        integer :: i
        integer(int8) :: key(storage_size(temp(1))/8)

        call map%init()
        allocate(mask(size(temp)))
        do i = 1, size(temp)
            key = 0
            key = transfer(temp(i), key)
            call map%key_test(key, present)
            if (.not. present) then
                call map%map_entry(key)
                mask(i) = .true.
            else
                mask(i) = .false.
            end if
        end do
        output = pack(temp, mask)
        deallocate(mask)
    end function

end submodule stdlib_sorting_unique