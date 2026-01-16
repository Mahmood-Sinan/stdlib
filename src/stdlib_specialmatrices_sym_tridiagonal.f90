submodule (stdlib_specialmatrices) sym_tridiagonal_matrices
    use stdlib_linalg_lapack, only: lagtm

    character(len=*), parameter :: this = "symmetric tridiagonal matrices"
    contains

    !--------------------------------
    !-----                      -----
    !-----     CONSTRUCTORS     -----
    !-----                      -----
    !--------------------------------

    pure module function initialize_sym_tridiagonal_pure_sp(du, dv) result(A)
        !! Construct a `symmetric tridiagonal` matrix from the rank-1 arrays `du` and `dv`.
        real(sp), intent(in) :: du(:), dv(:)
        !! symmetric tridiagonal matrix elements
        type(sym_tridiagonal_sp_type) :: A
        !! Corresponding symmetric tridiagonal matrix.

        ! Internal variables.
        integer(ilp):: n
        type(linalg_state_type) :: err0

        ! Sanity check
        n = size(dv, kind = ilp)
        if (n <= 0) then
            err0 = linalg_state_type(this, LINALG_VALUE_ERROR, "Matrix size needs to be positive, n = ", n, ".")
            call linalg_error_handling(err0)
        endif
        if(size(du, kind=ilp) /= n - 1) then
            err0 = linalg_state_type(this, LINALG_VALUE_ERROR, "Vector du does not have the correct length.")
            call linalg_error_handling(err0)
        endif
        
        ! Description of the matrix.
        A%n = n
        ! Matrix elements.
        A%du = du ; A%dv = dv
    end function

    pure module function initialize_constant_sym_tridiagonal_pure_sp(du, dv, n) result(A)
        !! Construct a `symmetric tridiagonal` matrix with constant elements.
        real(sp), intent(in) :: du, dv
        !! symmetric tridiagonal matrix elements.
        integer(ilp), intent(in) :: n
        !! Matrix dimension.
        type(sym_tridiagonal_sp_type) :: A
        !! Corresponding symmetric tridiagonal matrix.

        ! Internal variables.
        type(linalg_state_type) :: err0

        ! Description of the matrix.
        A%n = n
        if (n <= 0) then
            err0 = linalg_state_type(this, LINALG_VALUE_ERROR, "Matrix size needs to be positive, n = ", n, ".")
            call linalg_error_handling(err0)
        endif
        ! Matrix elements.
        allocate( A%dv(n), source= dv )
        allocate( A%du(n-1), source = du )
    end function

    module function initialize_sym_tridiagonal_impure_sp(du, dv, err) result(A)
        !! Construct a `symmetric tridiagonal` matrix from the rank-1 arrays `du` and `dv`.
        real(sp), intent(in) :: du(:), dv(:)
        !! symmetric tridiagonal matrix elements
        type(linalg_state_type), intent(out) :: err
        !! Error handling.
        type(sym_tridiagonal_sp_type) :: A
        !! Corresponding symmetric tridiagonal matrix.

        ! Internal variables.
        integer(ilp):: n
        type(linalg_state_type) :: err0

        ! Sanity check
        n = size(dv, kind = ilp)
        if (n <= 0) then
            err0 = linalg_state_type(this, LINALG_VALUE_ERROR, "Matrix size needs to be positive, n = ",n, ".")
            call linalg_error_handling(err0, err)
        endif
        if(size(du, kind=ilp) /= n - 1) then
            err0 = linalg_state_type(this, LINALG_VALUE_ERROR, "Vector du does not have the correct length.")
            call linalg_error_handling(err0, err)
        endif
        
        if(err0%ok()) then
            ! Description of the matrix.
            A%n = n
            ! Matrix elements.
            A%du = du ; A%dv = dv
        endif
    end function

    module function initialize_constant_sym_tridiagonal_impure_sp(du, dv, n, err) result(A)
        !! Construct a `symmetric tridiagonal` matrix with constant elements.
        real(sp), intent(in) :: du, dv
        !! symmetric tridiagonal matrix elements.
        integer(ilp), intent(in) :: n
        !! Matrix dimension.
        type(linalg_state_type), intent(out) :: err
        !! Error handling.
        type(sym_tridiagonal_sp_type) :: A
        !! Corresponding symmetric tridiagonal matrix.

        ! Internal variables.
        type(linalg_state_type) :: err0

        if (n <= 0) then
            err0 = linalg_state_type(this, LINALG_VALUE_ERROR, "Matrix size needs to be positive, n = ",n, ".")
            call linalg_error_handling(err0, err)
        endif

        if(err0%ok()) then
            ! Description of the matrix.
            A%n = n
            ! Matrix elements.
            allocate( A%dv(n), source= dv )
            allocate( A%du(n-1), source = du )
        endif
        end function
    pure module function initialize_sym_tridiagonal_pure_dp(du, dv) result(A)
        !! Construct a `symmetric tridiagonal` matrix from the rank-1 arrays `du` and `dv`.
        real(dp), intent(in) :: du(:), dv(:)
        !! symmetric tridiagonal matrix elements
        type(sym_tridiagonal_dp_type) :: A
        !! Corresponding symmetric tridiagonal matrix.

        ! Internal variables.
        integer(ilp):: n
        type(linalg_state_type) :: err0

        ! Sanity check
        n = size(dv, kind = ilp)
        if (n <= 0) then
            err0 = linalg_state_type(this, LINALG_VALUE_ERROR, "Matrix size needs to be positive, n = ", n, ".")
            call linalg_error_handling(err0)
        endif
        if(size(du, kind=ilp) /= n - 1) then
            err0 = linalg_state_type(this, LINALG_VALUE_ERROR, "Vector du does not have the correct length.")
            call linalg_error_handling(err0)
        endif
        
        ! Description of the matrix.
        A%n = n
        ! Matrix elements.
        A%du = du ; A%dv = dv
    end function

    pure module function initialize_constant_sym_tridiagonal_pure_dp(du, dv, n) result(A)
        !! Construct a `symmetric tridiagonal` matrix with constant elements.
        real(dp), intent(in) :: du, dv
        !! symmetric tridiagonal matrix elements.
        integer(ilp), intent(in) :: n
        !! Matrix dimension.
        type(sym_tridiagonal_dp_type) :: A
        !! Corresponding symmetric tridiagonal matrix.

        ! Internal variables.
        type(linalg_state_type) :: err0

        ! Description of the matrix.
        A%n = n
        if (n <= 0) then
            err0 = linalg_state_type(this, LINALG_VALUE_ERROR, "Matrix size needs to be positive, n = ", n, ".")
            call linalg_error_handling(err0)
        endif
        ! Matrix elements.
        allocate( A%dv(n), source= dv )
        allocate( A%du(n-1), source = du )
    end function

    module function initialize_sym_tridiagonal_impure_dp(du, dv, err) result(A)
        !! Construct a `symmetric tridiagonal` matrix from the rank-1 arrays `du` and `dv`.
        real(dp), intent(in) :: du(:), dv(:)
        !! symmetric tridiagonal matrix elements
        type(linalg_state_type), intent(out) :: err
        !! Error handling.
        type(sym_tridiagonal_dp_type) :: A
        !! Corresponding symmetric tridiagonal matrix.

        ! Internal variables.
        integer(ilp):: n
        type(linalg_state_type) :: err0

        ! Sanity check
        n = size(dv, kind = ilp)
        if (n <= 0) then
            err0 = linalg_state_type(this, LINALG_VALUE_ERROR, "Matrix size needs to be positive, n = ",n, ".")
            call linalg_error_handling(err0, err)
        endif
        if(size(du, kind=ilp) /= n - 1) then
            err0 = linalg_state_type(this, LINALG_VALUE_ERROR, "Vector du does not have the correct length.")
            call linalg_error_handling(err0, err)
        endif
        
        if(err0%ok()) then
            ! Description of the matrix.
            A%n = n
            ! Matrix elements.
            A%du = du ; A%dv = dv
        endif
    end function

    module function initialize_constant_sym_tridiagonal_impure_dp(du, dv, n, err) result(A)
        !! Construct a `symmetric tridiagonal` matrix with constant elements.
        real(dp), intent(in) :: du, dv
        !! symmetric tridiagonal matrix elements.
        integer(ilp), intent(in) :: n
        !! Matrix dimension.
        type(linalg_state_type), intent(out) :: err
        !! Error handling.
        type(sym_tridiagonal_dp_type) :: A
        !! Corresponding symmetric tridiagonal matrix.

        ! Internal variables.
        type(linalg_state_type) :: err0

        if (n <= 0) then
            err0 = linalg_state_type(this, LINALG_VALUE_ERROR, "Matrix size needs to be positive, n = ",n, ".")
            call linalg_error_handling(err0, err)
        endif

        if(err0%ok()) then
            ! Description of the matrix.
            A%n = n
            ! Matrix elements.
            allocate( A%dv(n), source= dv )
            allocate( A%du(n-1), source = du )
        endif
        end function
    pure module function initialize_sym_tridiagonal_pure_csp(du, dv) result(A)
        !! Construct a `symmetric tridiagonal` matrix from the rank-1 arrays `du` and `dv`.
        complex(sp), intent(in) :: du(:), dv(:)
        !! symmetric tridiagonal matrix elements
        type(sym_tridiagonal_csp_type) :: A
        !! Corresponding symmetric tridiagonal matrix.

        ! Internal variables.
        integer(ilp):: n
        type(linalg_state_type) :: err0

        ! Sanity check
        n = size(dv, kind = ilp)
        if (n <= 0) then
            err0 = linalg_state_type(this, LINALG_VALUE_ERROR, "Matrix size needs to be positive, n = ", n, ".")
            call linalg_error_handling(err0)
        endif
        if(size(du, kind=ilp) /= n - 1) then
            err0 = linalg_state_type(this, LINALG_VALUE_ERROR, "Vector du does not have the correct length.")
            call linalg_error_handling(err0)
        endif
        
        ! Description of the matrix.
        A%n = n
        ! Matrix elements.
        A%du = du ; A%dv = dv
    end function

    pure module function initialize_constant_sym_tridiagonal_pure_csp(du, dv, n) result(A)
        !! Construct a `symmetric tridiagonal` matrix with constant elements.
        complex(sp), intent(in) :: du, dv
        !! symmetric tridiagonal matrix elements.
        integer(ilp), intent(in) :: n
        !! Matrix dimension.
        type(sym_tridiagonal_csp_type) :: A
        !! Corresponding symmetric tridiagonal matrix.

        ! Internal variables.
        type(linalg_state_type) :: err0

        ! Description of the matrix.
        A%n = n
        if (n <= 0) then
            err0 = linalg_state_type(this, LINALG_VALUE_ERROR, "Matrix size needs to be positive, n = ", n, ".")
            call linalg_error_handling(err0)
        endif
        ! Matrix elements.
        allocate( A%dv(n), source= dv )
        allocate( A%du(n-1), source = du )
    end function

    module function initialize_sym_tridiagonal_impure_csp(du, dv, err) result(A)
        !! Construct a `symmetric tridiagonal` matrix from the rank-1 arrays `du` and `dv`.
        complex(sp), intent(in) :: du(:), dv(:)
        !! symmetric tridiagonal matrix elements
        type(linalg_state_type), intent(out) :: err
        !! Error handling.
        type(sym_tridiagonal_csp_type) :: A
        !! Corresponding symmetric tridiagonal matrix.

        ! Internal variables.
        integer(ilp):: n
        type(linalg_state_type) :: err0

        ! Sanity check
        n = size(dv, kind = ilp)
        if (n <= 0) then
            err0 = linalg_state_type(this, LINALG_VALUE_ERROR, "Matrix size needs to be positive, n = ",n, ".")
            call linalg_error_handling(err0, err)
        endif
        if(size(du, kind=ilp) /= n - 1) then
            err0 = linalg_state_type(this, LINALG_VALUE_ERROR, "Vector du does not have the correct length.")
            call linalg_error_handling(err0, err)
        endif
        
        if(err0%ok()) then
            ! Description of the matrix.
            A%n = n
            ! Matrix elements.
            A%du = du ; A%dv = dv
        endif
    end function

    module function initialize_constant_sym_tridiagonal_impure_csp(du, dv, n, err) result(A)
        !! Construct a `symmetric tridiagonal` matrix with constant elements.
        complex(sp), intent(in) :: du, dv
        !! symmetric tridiagonal matrix elements.
        integer(ilp), intent(in) :: n
        !! Matrix dimension.
        type(linalg_state_type), intent(out) :: err
        !! Error handling.
        type(sym_tridiagonal_csp_type) :: A
        !! Corresponding symmetric tridiagonal matrix.

        ! Internal variables.
        type(linalg_state_type) :: err0

        if (n <= 0) then
            err0 = linalg_state_type(this, LINALG_VALUE_ERROR, "Matrix size needs to be positive, n = ",n, ".")
            call linalg_error_handling(err0, err)
        endif

        if(err0%ok()) then
            ! Description of the matrix.
            A%n = n
            ! Matrix elements.
            allocate( A%dv(n), source= dv )
            allocate( A%du(n-1), source = du )
        endif
        end function
    pure module function initialize_sym_tridiagonal_pure_cdp(du, dv) result(A)
        !! Construct a `symmetric tridiagonal` matrix from the rank-1 arrays `du` and `dv`.
        complex(dp), intent(in) :: du(:), dv(:)
        !! symmetric tridiagonal matrix elements
        type(sym_tridiagonal_cdp_type) :: A
        !! Corresponding symmetric tridiagonal matrix.

        ! Internal variables.
        integer(ilp):: n
        type(linalg_state_type) :: err0

        ! Sanity check
        n = size(dv, kind = ilp)
        if (n <= 0) then
            err0 = linalg_state_type(this, LINALG_VALUE_ERROR, "Matrix size needs to be positive, n = ", n, ".")
            call linalg_error_handling(err0)
        endif
        if(size(du, kind=ilp) /= n - 1) then
            err0 = linalg_state_type(this, LINALG_VALUE_ERROR, "Vector du does not have the correct length.")
            call linalg_error_handling(err0)
        endif
        
        ! Description of the matrix.
        A%n = n
        ! Matrix elements.
        A%du = du ; A%dv = dv
    end function

    pure module function initialize_constant_sym_tridiagonal_pure_cdp(du, dv, n) result(A)
        !! Construct a `symmetric tridiagonal` matrix with constant elements.
        complex(dp), intent(in) :: du, dv
        !! symmetric tridiagonal matrix elements.
        integer(ilp), intent(in) :: n
        !! Matrix dimension.
        type(sym_tridiagonal_cdp_type) :: A
        !! Corresponding symmetric tridiagonal matrix.

        ! Internal variables.
        type(linalg_state_type) :: err0

        ! Description of the matrix.
        A%n = n
        if (n <= 0) then
            err0 = linalg_state_type(this, LINALG_VALUE_ERROR, "Matrix size needs to be positive, n = ", n, ".")
            call linalg_error_handling(err0)
        endif
        ! Matrix elements.
        allocate( A%dv(n), source= dv )
        allocate( A%du(n-1), source = du )
    end function

    module function initialize_sym_tridiagonal_impure_cdp(du, dv, err) result(A)
        !! Construct a `symmetric tridiagonal` matrix from the rank-1 arrays `du` and `dv`.
        complex(dp), intent(in) :: du(:), dv(:)
        !! symmetric tridiagonal matrix elements
        type(linalg_state_type), intent(out) :: err
        !! Error handling.
        type(sym_tridiagonal_cdp_type) :: A
        !! Corresponding symmetric tridiagonal matrix.

        ! Internal variables.
        integer(ilp):: n
        type(linalg_state_type) :: err0

        ! Sanity check
        n = size(dv, kind = ilp)
        if (n <= 0) then
            err0 = linalg_state_type(this, LINALG_VALUE_ERROR, "Matrix size needs to be positive, n = ",n, ".")
            call linalg_error_handling(err0, err)
        endif
        if(size(du, kind=ilp) /= n - 1) then
            err0 = linalg_state_type(this, LINALG_VALUE_ERROR, "Vector du does not have the correct length.")
            call linalg_error_handling(err0, err)
        endif
        
        if(err0%ok()) then
            ! Description of the matrix.
            A%n = n
            ! Matrix elements.
            A%du = du ; A%dv = dv
        endif
    end function

    module function initialize_constant_sym_tridiagonal_impure_cdp(du, dv, n, err) result(A)
        !! Construct a `symmetric tridiagonal` matrix with constant elements.
        complex(dp), intent(in) :: du, dv
        !! symmetric tridiagonal matrix elements.
        integer(ilp), intent(in) :: n
        !! Matrix dimension.
        type(linalg_state_type), intent(out) :: err
        !! Error handling.
        type(sym_tridiagonal_cdp_type) :: A
        !! Corresponding symmetric tridiagonal matrix.

        ! Internal variables.
        type(linalg_state_type) :: err0

        if (n <= 0) then
            err0 = linalg_state_type(this, LINALG_VALUE_ERROR, "Matrix size needs to be positive, n = ",n, ".")
            call linalg_error_handling(err0, err)
        endif

        if(err0%ok()) then
            ! Description of the matrix.
            A%n = n
            ! Matrix elements.
            allocate( A%dv(n), source= dv )
            allocate( A%du(n-1), source = du )
        endif
        end function
    
    !-----------------------------------------
    !-----                               -----
    !-----     MATRIX-VECTOR PRODUCT     -----
    !-----                               -----
    !-----------------------------------------

    !! spmv_sym_tridiag
    module subroutine spmv_sym_tridiag_1d_sp(A, x, y, alpha, beta, op)
        type(sym_tridiagonal_sp_type), intent(in) :: A
        real(sp), intent(in), contiguous, target :: x(:)
        real(sp), intent(inout), contiguous, target :: y(:)
        real(sp), intent(in), optional :: alpha
        real(sp), intent(in), optional :: beta
        character(1), intent(in), optional :: op

        ! Internal variables.
        real(sp) :: alpha_, beta_
        integer(ilp) :: n, nrhs, ldx, ldy
        character(1) :: op_
        logical :: is_alpha_special, is_beta_special

        real(sp), pointer :: xmat(:, :), ymat(:, :)

        ! Deal with optional arguments.
        alpha_ = 1.0_sp ; if (present(alpha)) alpha_ = alpha
        beta_  = 0.0_sp ; if (present(beta))  beta_  = beta
        op_    = "N"        ; if (present(op))    op_    = op
        is_alpha_special = (alpha_ ==  1.0_sp  .or. alpha_ ==  0.0_sp  .or. alpha_ == -1.0_sp)
        is_beta_special  = (beta_  ==  1.0_sp  .or. beta_  ==  0.0_sp  .or. beta_  == -1.0_sp)

        ! Prepare Lapack arguments.
        n = A%n ; ldx = n ; ldy = n ;
        nrhs =  1 

        ! Pointer trick.
        xmat(1:n, 1:nrhs) => x ; ymat(1:n, 1:nrhs) => y
        if(is_alpha_special .and. is_beta_special) then
            call lagtm(op_, n, nrhs, alpha_, A%du, A%dv, A%du, xmat, ldx, beta_, ymat, ldy)
        else
            call glagtm(op_, n, nrhs, alpha_, A%du, A%dv, A%du, xmat, ldx, beta_, ymat, ldy)
        end if
    end subroutine
    module subroutine spmv_sym_tridiag_2d_sp(A, x, y, alpha, beta, op)
        type(sym_tridiagonal_sp_type), intent(in) :: A
        real(sp), intent(in), contiguous, target :: x(:,:)
        real(sp), intent(inout), contiguous, target :: y(:,:)
        real(sp), intent(in), optional :: alpha
        real(sp), intent(in), optional :: beta
        character(1), intent(in), optional :: op

        ! Internal variables.
        real(sp) :: alpha_, beta_
        integer(ilp) :: n, nrhs, ldx, ldy
        character(1) :: op_
        logical :: is_alpha_special, is_beta_special


        ! Deal with optional arguments.
        alpha_ = 1.0_sp ; if (present(alpha)) alpha_ = alpha
        beta_  = 0.0_sp ; if (present(beta))  beta_  = beta
        op_    = "N"        ; if (present(op))    op_    = op
        is_alpha_special = (alpha_ ==  1.0_sp  .or. alpha_ ==  0.0_sp  .or. alpha_ == -1.0_sp)
        is_beta_special  = (beta_  ==  1.0_sp  .or. beta_  ==  0.0_sp  .or. beta_  == -1.0_sp)

        ! Prepare Lapack arguments.
        n = A%n ; ldx = n ; ldy = n ;
        nrhs =  size(x, dim=2, kind=ilp) 

        if(is_alpha_special .and. is_beta_special) then
            call lagtm(op_, n, nrhs, alpha_, A%du, A%dv, A%du, x, ldx, beta_, y, ldy)
        else
            call glagtm(op_, n, nrhs, alpha_, A%du, A%dv, A%du, x, ldx, beta_, y, ldy)
        end if
    end subroutine
    module subroutine spmv_sym_tridiag_1d_dp(A, x, y, alpha, beta, op)
        type(sym_tridiagonal_dp_type), intent(in) :: A
        real(dp), intent(in), contiguous, target :: x(:)
        real(dp), intent(inout), contiguous, target :: y(:)
        real(dp), intent(in), optional :: alpha
        real(dp), intent(in), optional :: beta
        character(1), intent(in), optional :: op

        ! Internal variables.
        real(dp) :: alpha_, beta_
        integer(ilp) :: n, nrhs, ldx, ldy
        character(1) :: op_
        logical :: is_alpha_special, is_beta_special

        real(dp), pointer :: xmat(:, :), ymat(:, :)

        ! Deal with optional arguments.
        alpha_ = 1.0_dp ; if (present(alpha)) alpha_ = alpha
        beta_  = 0.0_dp ; if (present(beta))  beta_  = beta
        op_    = "N"        ; if (present(op))    op_    = op
        is_alpha_special = (alpha_ ==  1.0_dp  .or. alpha_ ==  0.0_dp  .or. alpha_ == -1.0_dp)
        is_beta_special  = (beta_  ==  1.0_dp  .or. beta_  ==  0.0_dp  .or. beta_  == -1.0_dp)

        ! Prepare Lapack arguments.
        n = A%n ; ldx = n ; ldy = n ;
        nrhs =  1 

        ! Pointer trick.
        xmat(1:n, 1:nrhs) => x ; ymat(1:n, 1:nrhs) => y
        if(is_alpha_special .and. is_beta_special) then
            call lagtm(op_, n, nrhs, alpha_, A%du, A%dv, A%du, xmat, ldx, beta_, ymat, ldy)
        else
            call glagtm(op_, n, nrhs, alpha_, A%du, A%dv, A%du, xmat, ldx, beta_, ymat, ldy)
        end if
    end subroutine
    module subroutine spmv_sym_tridiag_2d_dp(A, x, y, alpha, beta, op)
        type(sym_tridiagonal_dp_type), intent(in) :: A
        real(dp), intent(in), contiguous, target :: x(:,:)
        real(dp), intent(inout), contiguous, target :: y(:,:)
        real(dp), intent(in), optional :: alpha
        real(dp), intent(in), optional :: beta
        character(1), intent(in), optional :: op

        ! Internal variables.
        real(dp) :: alpha_, beta_
        integer(ilp) :: n, nrhs, ldx, ldy
        character(1) :: op_
        logical :: is_alpha_special, is_beta_special


        ! Deal with optional arguments.
        alpha_ = 1.0_dp ; if (present(alpha)) alpha_ = alpha
        beta_  = 0.0_dp ; if (present(beta))  beta_  = beta
        op_    = "N"        ; if (present(op))    op_    = op
        is_alpha_special = (alpha_ ==  1.0_dp  .or. alpha_ ==  0.0_dp  .or. alpha_ == -1.0_dp)
        is_beta_special  = (beta_  ==  1.0_dp  .or. beta_  ==  0.0_dp  .or. beta_  == -1.0_dp)

        ! Prepare Lapack arguments.
        n = A%n ; ldx = n ; ldy = n ;
        nrhs =  size(x, dim=2, kind=ilp) 

        if(is_alpha_special .and. is_beta_special) then
            call lagtm(op_, n, nrhs, alpha_, A%du, A%dv, A%du, x, ldx, beta_, y, ldy)
        else
            call glagtm(op_, n, nrhs, alpha_, A%du, A%dv, A%du, x, ldx, beta_, y, ldy)
        end if
    end subroutine
    module subroutine spmv_sym_tridiag_1d_csp(A, x, y, alpha, beta, op)
        type(sym_tridiagonal_csp_type), intent(in) :: A
        complex(sp), intent(in), contiguous, target :: x(:)
        complex(sp), intent(inout), contiguous, target :: y(:)
        complex(sp), intent(in), optional :: alpha
        complex(sp), intent(in), optional :: beta
        character(1), intent(in), optional :: op

        ! Internal variables.
        complex(sp) :: alpha_, beta_
        integer(ilp) :: n, nrhs, ldx, ldy
        character(1) :: op_

        complex(sp), pointer :: xmat(:, :), ymat(:, :)

        ! Deal with optional arguments.
        alpha_ = 1.0_sp ; if (present(alpha)) alpha_ = alpha
        beta_  = 0.0_sp ; if (present(beta))  beta_  = beta
        op_    = "N"        ; if (present(op))    op_    = op

        ! Prepare Lapack arguments.
        n = A%n ; ldx = n ; ldy = n ;
        nrhs =  1 

        ! Pointer trick.
        xmat(1:n, 1:nrhs) => x ; ymat(1:n, 1:nrhs) => y
        call glagtm(op_, n, nrhs, alpha_, A%du, A%dv, A%du, xmat, ldx, beta_, ymat, ldy)
    end subroutine
    module subroutine spmv_sym_tridiag_2d_csp(A, x, y, alpha, beta, op)
        type(sym_tridiagonal_csp_type), intent(in) :: A
        complex(sp), intent(in), contiguous, target :: x(:,:)
        complex(sp), intent(inout), contiguous, target :: y(:,:)
        complex(sp), intent(in), optional :: alpha
        complex(sp), intent(in), optional :: beta
        character(1), intent(in), optional :: op

        ! Internal variables.
        complex(sp) :: alpha_, beta_
        integer(ilp) :: n, nrhs, ldx, ldy
        character(1) :: op_


        ! Deal with optional arguments.
        alpha_ = 1.0_sp ; if (present(alpha)) alpha_ = alpha
        beta_  = 0.0_sp ; if (present(beta))  beta_  = beta
        op_    = "N"        ; if (present(op))    op_    = op

        ! Prepare Lapack arguments.
        n = A%n ; ldx = n ; ldy = n ;
        nrhs =  size(x, dim=2, kind=ilp) 

        call glagtm(op_, n, nrhs, alpha_, A%du, A%dv, A%du, x, ldx, beta_, y, ldy)
    end subroutine
    module subroutine spmv_sym_tridiag_1d_cdp(A, x, y, alpha, beta, op)
        type(sym_tridiagonal_cdp_type), intent(in) :: A
        complex(dp), intent(in), contiguous, target :: x(:)
        complex(dp), intent(inout), contiguous, target :: y(:)
        complex(dp), intent(in), optional :: alpha
        complex(dp), intent(in), optional :: beta
        character(1), intent(in), optional :: op

        ! Internal variables.
        complex(dp) :: alpha_, beta_
        integer(ilp) :: n, nrhs, ldx, ldy
        character(1) :: op_

        complex(dp), pointer :: xmat(:, :), ymat(:, :)

        ! Deal with optional arguments.
        alpha_ = 1.0_dp ; if (present(alpha)) alpha_ = alpha
        beta_  = 0.0_dp ; if (present(beta))  beta_  = beta
        op_    = "N"        ; if (present(op))    op_    = op

        ! Prepare Lapack arguments.
        n = A%n ; ldx = n ; ldy = n ;
        nrhs =  1 

        ! Pointer trick.
        xmat(1:n, 1:nrhs) => x ; ymat(1:n, 1:nrhs) => y
        call glagtm(op_, n, nrhs, alpha_, A%du, A%dv, A%du, xmat, ldx, beta_, ymat, ldy)
    end subroutine
    module subroutine spmv_sym_tridiag_2d_cdp(A, x, y, alpha, beta, op)
        type(sym_tridiagonal_cdp_type), intent(in) :: A
        complex(dp), intent(in), contiguous, target :: x(:,:)
        complex(dp), intent(inout), contiguous, target :: y(:,:)
        complex(dp), intent(in), optional :: alpha
        complex(dp), intent(in), optional :: beta
        character(1), intent(in), optional :: op

        ! Internal variables.
        complex(dp) :: alpha_, beta_
        integer(ilp) :: n, nrhs, ldx, ldy
        character(1) :: op_


        ! Deal with optional arguments.
        alpha_ = 1.0_dp ; if (present(alpha)) alpha_ = alpha
        beta_  = 0.0_dp ; if (present(beta))  beta_  = beta
        op_    = "N"        ; if (present(op))    op_    = op

        ! Prepare Lapack arguments.
        n = A%n ; ldx = n ; ldy = n ;
        nrhs =  size(x, dim=2, kind=ilp) 

        call glagtm(op_, n, nrhs, alpha_, A%du, A%dv, A%du, x, ldx, beta_, y, ldy)
    end subroutine

    !-------------------------------------
    !-----                           -----
    !-----     UTILITY FUNCTIONS     -----
    !-----                           -----
    !-------------------------------------
    pure module function sym_tridiagonal_to_dense_sp(A) result(B)
        !! Convert a `symmetric tridiagonal` matrix to its dense representation.
        type(sym_tridiagonal_sp_type), intent(in) :: A
        !! Input matrix.
        real(sp), allocatable :: B(:, :)
        !! Corresponding dense matrix.

        ! Internal variables.
        integer(ilp) :: i

        associate (n => A%n)
            allocate(B(n,n), source=zero_sp)
            B(1,1) = A%dv(1) ; B(1,2) = A%du(1)
            do concurrent (i = 2: n - 1)
                B(i, i - 1) = A%du(i - 1)
                B(i, i) = A%dv(i)
                B(i, i + 1) = A%du(i)
            enddo
            B(n , n -1) = A%du(n - 1); B(n, n) = A%dv(n)
        end associate
    end function
    pure module function sym_tridiagonal_to_dense_dp(A) result(B)
        !! Convert a `symmetric tridiagonal` matrix to its dense representation.
        type(sym_tridiagonal_dp_type), intent(in) :: A
        !! Input matrix.
        real(dp), allocatable :: B(:, :)
        !! Corresponding dense matrix.

        ! Internal variables.
        integer(ilp) :: i

        associate (n => A%n)
            allocate(B(n,n), source=zero_dp)
            B(1,1) = A%dv(1) ; B(1,2) = A%du(1)
            do concurrent (i = 2: n - 1)
                B(i, i - 1) = A%du(i - 1)
                B(i, i) = A%dv(i)
                B(i, i + 1) = A%du(i)
            enddo
            B(n , n -1) = A%du(n - 1); B(n, n) = A%dv(n)
        end associate
    end function
    pure module function sym_tridiagonal_to_dense_csp(A) result(B)
        !! Convert a `symmetric tridiagonal` matrix to its dense representation.
        type(sym_tridiagonal_csp_type), intent(in) :: A
        !! Input matrix.
        complex(sp), allocatable :: B(:, :)
        !! Corresponding dense matrix.

        ! Internal variables.
        integer(ilp) :: i

        associate (n => A%n)
            allocate(B(n,n), source=zero_csp)
            B(1,1) = A%dv(1) ; B(1,2) = A%du(1)
            do concurrent (i = 2: n - 1)
                B(i, i - 1) = A%du(i - 1)
                B(i, i) = A%dv(i)
                B(i, i + 1) = A%du(i)
            enddo
            B(n , n -1) = A%du(n - 1); B(n, n) = A%dv(n)
        end associate
    end function
    pure module function sym_tridiagonal_to_dense_cdp(A) result(B)
        !! Convert a `symmetric tridiagonal` matrix to its dense representation.
        type(sym_tridiagonal_cdp_type), intent(in) :: A
        !! Input matrix.
        complex(dp), allocatable :: B(:, :)
        !! Corresponding dense matrix.

        ! Internal variables.
        integer(ilp) :: i

        associate (n => A%n)
            allocate(B(n,n), source=zero_cdp)
            B(1,1) = A%dv(1) ; B(1,2) = A%du(1)
            do concurrent (i = 2: n - 1)
                B(i, i - 1) = A%du(i - 1)
                B(i, i) = A%dv(i)
                B(i, i + 1) = A%du(i)
            enddo
            B(n , n -1) = A%du(n - 1); B(n, n) = A%dv(n)
        end associate
    end function

    pure module function transpose_sym_tridiagonal_sp(A) result(B)
        type(sym_tridiagonal_sp_type), intent(in) :: A
        !! Input matrix.
        type(sym_tridiagonal_sp_type) :: B
        B = sym_tridiagonal(A%du, A%dv)
    end function
    pure module function transpose_sym_tridiagonal_dp(A) result(B)
        type(sym_tridiagonal_dp_type), intent(in) :: A
        !! Input matrix.
        type(sym_tridiagonal_dp_type) :: B
        B = sym_tridiagonal(A%du, A%dv)
    end function
    pure module function transpose_sym_tridiagonal_csp(A) result(B)
        type(sym_tridiagonal_csp_type), intent(in) :: A
        !! Input matrix.
        type(sym_tridiagonal_csp_type) :: B
        B = sym_tridiagonal(A%du, A%dv)
    end function
    pure module function transpose_sym_tridiagonal_cdp(A) result(B)
        type(sym_tridiagonal_cdp_type), intent(in) :: A
        !! Input matrix.
        type(sym_tridiagonal_cdp_type) :: B
        B = sym_tridiagonal(A%du, A%dv)
    end function

    pure module function hermitian_sym_tridiagonal_sp(A) result(B)
        type(sym_tridiagonal_sp_type), intent(in) :: A
        !! Input matrix.
        type(sym_tridiagonal_sp_type) :: B
        B = sym_tridiagonal(A%du, A%dv)
    end function
    pure module function hermitian_sym_tridiagonal_dp(A) result(B)
        type(sym_tridiagonal_dp_type), intent(in) :: A
        !! Input matrix.
        type(sym_tridiagonal_dp_type) :: B
        B = sym_tridiagonal(A%du, A%dv)
    end function
    pure module function hermitian_sym_tridiagonal_csp(A) result(B)
        type(sym_tridiagonal_csp_type), intent(in) :: A
        !! Input matrix.
        type(sym_tridiagonal_csp_type) :: B
        B = sym_tridiagonal(conjg(A%du), conjg(A%dv))
    end function
    pure module function hermitian_sym_tridiagonal_cdp(A) result(B)
        type(sym_tridiagonal_cdp_type), intent(in) :: A
        !! Input matrix.
        type(sym_tridiagonal_cdp_type) :: B
        B = sym_tridiagonal(conjg(A%du), conjg(A%dv))
    end function

    pure module function scalar_multiplication_sym_tridiagonal_sp(alpha, A) result(B)
        real(sp), intent(in) :: alpha
        type(sym_tridiagonal_sp_type), intent(in) :: A
        type(sym_tridiagonal_sp_type) :: B
        B = sym_tridiagonal(A%du, A%dv)
        B%du = alpha*B%du; B%dv = alpha*B%dv;
    end function

    pure module function scalar_multiplication_bis_sym_tridiagonal_sp(A, alpha) result(B)
        type(sym_tridiagonal_sp_type), intent(in) :: A
        real(sp), intent(in) :: alpha
        type(sym_tridiagonal_sp_type) :: B
        B = sym_tridiagonal(A%du, A%dv)
        B%du = alpha*B%du; B%dv = alpha*B%dv;
    end function
    pure module function scalar_multiplication_sym_tridiagonal_dp(alpha, A) result(B)
        real(dp), intent(in) :: alpha
        type(sym_tridiagonal_dp_type), intent(in) :: A
        type(sym_tridiagonal_dp_type) :: B
        B = sym_tridiagonal(A%du, A%dv)
        B%du = alpha*B%du; B%dv = alpha*B%dv;
    end function

    pure module function scalar_multiplication_bis_sym_tridiagonal_dp(A, alpha) result(B)
        type(sym_tridiagonal_dp_type), intent(in) :: A
        real(dp), intent(in) :: alpha
        type(sym_tridiagonal_dp_type) :: B
        B = sym_tridiagonal(A%du, A%dv)
        B%du = alpha*B%du; B%dv = alpha*B%dv;
    end function
    pure module function scalar_multiplication_sym_tridiagonal_csp(alpha, A) result(B)
        complex(sp), intent(in) :: alpha
        type(sym_tridiagonal_csp_type), intent(in) :: A
        type(sym_tridiagonal_csp_type) :: B
        B = sym_tridiagonal(A%du, A%dv)
        B%du = alpha*B%du; B%dv = alpha*B%dv;
    end function

    pure module function scalar_multiplication_bis_sym_tridiagonal_csp(A, alpha) result(B)
        type(sym_tridiagonal_csp_type), intent(in) :: A
        complex(sp), intent(in) :: alpha
        type(sym_tridiagonal_csp_type) :: B
        B = sym_tridiagonal(A%du, A%dv)
        B%du = alpha*B%du; B%dv = alpha*B%dv;
    end function
    pure module function scalar_multiplication_sym_tridiagonal_cdp(alpha, A) result(B)
        complex(dp), intent(in) :: alpha
        type(sym_tridiagonal_cdp_type), intent(in) :: A
        type(sym_tridiagonal_cdp_type) :: B
        B = sym_tridiagonal(A%du, A%dv)
        B%du = alpha*B%du; B%dv = alpha*B%dv;
    end function

    pure module function scalar_multiplication_bis_sym_tridiagonal_cdp(A, alpha) result(B)
        type(sym_tridiagonal_cdp_type), intent(in) :: A
        complex(dp), intent(in) :: alpha
        type(sym_tridiagonal_cdp_type) :: B
        B = sym_tridiagonal(A%du, A%dv)
        B%du = alpha*B%du; B%dv = alpha*B%dv;
    end function

    pure module function matrix_add_sym_tridiagonal_sp(A, B) result(C)
        type(sym_tridiagonal_sp_type), intent(in) :: A
        type(sym_tridiagonal_sp_type), intent(in) :: B
        type(sym_tridiagonal_sp_type) :: C
        C = sym_tridiagonal(A%du, A%dv)
        C%du = C%du + B%du; C%dv = C%dv + B%dv
    end function

    pure module function matrix_sub_sym_tridiagonal_sp(A, B) result(C)
        type(sym_tridiagonal_sp_type), intent(in) :: A
        type(sym_tridiagonal_sp_type), intent(in) :: B
        type(sym_tridiagonal_sp_type) :: C
        C = sym_tridiagonal(A%du, A%dv)
        C%du = C%du - B%du; C%dv = C%dv - B%dv
    end function
    pure module function matrix_add_sym_tridiagonal_dp(A, B) result(C)
        type(sym_tridiagonal_dp_type), intent(in) :: A
        type(sym_tridiagonal_dp_type), intent(in) :: B
        type(sym_tridiagonal_dp_type) :: C
        C = sym_tridiagonal(A%du, A%dv)
        C%du = C%du + B%du; C%dv = C%dv + B%dv
    end function

    pure module function matrix_sub_sym_tridiagonal_dp(A, B) result(C)
        type(sym_tridiagonal_dp_type), intent(in) :: A
        type(sym_tridiagonal_dp_type), intent(in) :: B
        type(sym_tridiagonal_dp_type) :: C
        C = sym_tridiagonal(A%du, A%dv)
        C%du = C%du - B%du; C%dv = C%dv - B%dv
    end function
    pure module function matrix_add_sym_tridiagonal_csp(A, B) result(C)
        type(sym_tridiagonal_csp_type), intent(in) :: A
        type(sym_tridiagonal_csp_type), intent(in) :: B
        type(sym_tridiagonal_csp_type) :: C
        C = sym_tridiagonal(A%du, A%dv)
        C%du = C%du + B%du; C%dv = C%dv + B%dv
    end function

    pure module function matrix_sub_sym_tridiagonal_csp(A, B) result(C)
        type(sym_tridiagonal_csp_type), intent(in) :: A
        type(sym_tridiagonal_csp_type), intent(in) :: B
        type(sym_tridiagonal_csp_type) :: C
        C = sym_tridiagonal(A%du, A%dv)
        C%du = C%du - B%du; C%dv = C%dv - B%dv
    end function
    pure module function matrix_add_sym_tridiagonal_cdp(A, B) result(C)
        type(sym_tridiagonal_cdp_type), intent(in) :: A
        type(sym_tridiagonal_cdp_type), intent(in) :: B
        type(sym_tridiagonal_cdp_type) :: C
        C = sym_tridiagonal(A%du, A%dv)
        C%du = C%du + B%du; C%dv = C%dv + B%dv
    end function

    pure module function matrix_sub_sym_tridiagonal_cdp(A, B) result(C)
        type(sym_tridiagonal_cdp_type), intent(in) :: A
        type(sym_tridiagonal_cdp_type), intent(in) :: B
        type(sym_tridiagonal_cdp_type) :: C
        C = sym_tridiagonal(A%du, A%dv)
        C%du = C%du - B%du; C%dv = C%dv - B%dv
    end function

end submodule