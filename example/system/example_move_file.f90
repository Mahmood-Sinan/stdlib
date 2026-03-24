program example_move_file
    use stdlib_system, only: move_file
    use stdlib_error,  only: state_type
    implicit none

    type(state_type) :: err
    integer :: ftype

    character(len=*), parameter :: src  = "example.txt"
    character(len=*), parameter :: dest = "d1/d2/d3/example.txt"

    ! Step 3: Move the file
    call move_file(src, dest, err)

    if (err%error()) then
        print *, "Error moving file:", err%print()
        stop
    else
        print *, "File moved successfully!"
    end if

end program example_move_file
