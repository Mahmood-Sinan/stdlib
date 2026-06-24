submodule(stdlib_spatial) stdlib_spatial_kabsch_umeyama
    use stdlib_linalg, only: svd, det
    use stdlib_intrinsics, only: stdlib_sum_kahan, stdlib_dot_product_kahan, kahan_kernel
    use stdlib_error, only: error_stop
    use stdlib_optval, only: optval
    use stdlib_linalg_lapack, only: gemm, gemv

contains
    module subroutine kabsch_umeyama_sp(P, Q, R, t, c, rmsd, W, scale)
        real(sp), intent(in) :: P(:, :)
        !! Target point set (d × N)
        real(sp), intent(in) :: Q(:, :)
        !! Reference point set (d × N)
        real(sp), intent(out) :: R(:, :)
        !! Optimal rotation matrix (d × d)
        real(sp), intent(out) :: t(:)
        !! Translation vector (d)
        real(sp), intent(out) :: c
        !! Scale factor
        real(sp), intent(out) :: rmsd
        !! Root-mean-square deviation
        real(sp), intent(in), optional :: W(:)
        !! Optional weights
        logical, intent(in), optional :: scale
        !! Enable scaling

        ! Internal variables.
        integer :: i, j, point, d, N
        real(sp), allocatable :: covariance(:,:), U(:,:), Vt(:,:), vec(:), tmp_N(:), tmp_d(:), c_P(:), c_Q(:)
        real(sp) :: sum_w, variance_p
        real(sp), allocatable :: S(:)
        real(sp) :: temp
        logical :: scale_
        logical :: reflect_
        real(sp) :: rmsd_err

        scale_ = optval(scale, .true.)
        ! Dimension checks
        d = size(P,dim=1)
        N = size(P,dim=2)
        if(any(shape(P)/=shape(Q)) .or. any(shape(R)/=[d,d]) .or. size(t)/=d) then
            call error_stop("array sizes do not match")
        end if
        if (present(W)) then
            if (size(W) /= N) then
                call error_stop("array sizes do not match")
            end if
        end if

        if(present(W)) then
            sum_w = stdlib_sum_kahan(W)
        else
            sum_w = real(N, kind = sp)
        end if
        !> leave opportunity for future discussion on how to add spmd reduction needed here to reduce sum_w before division

        sum_w = one_sp / sum_w
        if(sum_w<zero_sp) call error_stop("Invalid weights: sum of weights must be positive")

        allocate(c_P(d), c_Q(d), tmp_N(N), source=zero_sp)

        ! Compute centroids of P and Q
        if(present(W)) then
            do i = 1, d
                tmp_N(:) = W(:) * P(i,:)
                c_P(i) = stdlib_sum_kahan(tmp_N)
                tmp_N(:) = W(:) * Q(i,:)
                c_Q(i) = stdlib_sum_kahan(tmp_N)
            end do
        else
            c_P = stdlib_sum_kahan(P, dim=2)
            c_Q = stdlib_sum_kahan(Q, dim=2)
        end if
        c_P = c_P * sum_w
        c_Q = c_Q * sum_w

        ! Compute covariance matrix H = (P - c_P) * (Q - c_Q)^T and variance of P
        allocate(covariance(d,d), source=zero_sp)
        allocate(tmp_d(d), source=zero_sp)
        variance_p = zero_sp

        if (present(W)) then
            do point = 1, N
                tmp_d = P(:, point) - c_P(:)
                tmp_N(point) = stdlib_dot_product_kahan(tmp_d, tmp_d)
            end do
            tmp_N(:) = W(:) * tmp_N(:)
            variance_p = stdlib_sum_kahan(tmp_N)
            do j = 1, d
                do i = 1, d
                    tmp_N(:) = W(:) * (P(i,:) - c_P(i)) * (Q(j,:) - c_Q(j))
                    covariance(i,j) = stdlib_sum_kahan(tmp_N)
                end do
            end do
        else
            ! Calculate variance by the formula (1/n)*sigma(P - c_P)^2
            do point = 1, N
                tmp_d = P(:, point) - c_P(:)
                tmp_N(point) = stdlib_dot_product_kahan(tmp_d, tmp_d)
            end do
            variance_p = stdlib_sum_kahan(tmp_N)
            do j = 1, d
                do i = 1, d
                    covariance(i,j) = stdlib_dot_product_kahan((P(i,:) - c_P(i)),(Q(j,:) - c_Q(j)))
                end do
            end do
        end if

        covariance = covariance * sum_w
        variance_p = variance_p * sum_w

        allocate(U(d,d), source=zero_sp)
        allocate(Vt(d,d), source=zero_sp)
        allocate(S(d), source=zero_sp)

        ! SVD of covariance matrix H -> H = U * S * Vt
        call svd(covariance, S, U, Vt)

        ! Check for reflections in case of real entries.
        reflect_ = det(matmul(U,Vt)) < zero_sp
        if(reflect_) Vt(d,:) = -Vt(d,:)

        ! Optimal rotation matrix.
        call gemm(transa='N', transb='N', m=d,n=d,k=d, alpha=one_sp, a=U,lda=d, b=Vt, ldb=d, beta=zero_sp, c=R, ldc=d)

        ! Scaling factor
        c = one_sp
        if(scale_) then
            if(reflect_) then
                c = sum(S(1:d-1)) - S(d)
            else
                c = sum(S(1:d))
            end if
            c = variance_p / c
        end if

        ! Translation vector t = c_P - c*R*c_Q
        t = c_P
        call gemv(trans='N', m=d, n=d, alpha=-c, A=R, lda=d, x=c_Q, incx=1, beta=one_sp, y=t, incy=1)

        ! Compute RMSD
        allocate(vec(d), source=zero_sp)
        rmsd = zero_sp
        rmsd_err = zero_sp
        do point = 1, N
            ! Calculate the k^th difference vector by the formula vec_k = c*R*Q_k + t - P_k
            vec = t
            call gemv(trans='N', m=d, n=d, alpha=c, A=R, lda=d, x=Q(:, point), incx=1, beta=one_sp, y=vec, incy=1)
            vec = vec - P(:,point)
            temp = stdlib_dot_product_kahan(vec,vec)
            call kahan_kernel(temp, rmsd, rmsd_err)
        end do
        rmsd = sqrt(rmsd * sum_w)
    end subroutine
    module subroutine kabsch_umeyama_dp(P, Q, R, t, c, rmsd, W, scale)
        real(dp), intent(in) :: P(:, :)
        !! Target point set (d × N)
        real(dp), intent(in) :: Q(:, :)
        !! Reference point set (d × N)
        real(dp), intent(out) :: R(:, :)
        !! Optimal rotation matrix (d × d)
        real(dp), intent(out) :: t(:)
        !! Translation vector (d)
        real(dp), intent(out) :: c
        !! Scale factor
        real(dp), intent(out) :: rmsd
        !! Root-mean-square deviation
        real(dp), intent(in), optional :: W(:)
        !! Optional weights
        logical, intent(in), optional :: scale
        !! Enable scaling

        ! Internal variables.
        integer :: i, j, point, d, N
        real(dp), allocatable :: covariance(:,:), U(:,:), Vt(:,:), vec(:), tmp_N(:), tmp_d(:), c_P(:), c_Q(:)
        real(dp) :: sum_w, variance_p
        real(dp), allocatable :: S(:)
        real(dp) :: temp
        logical :: scale_
        logical :: reflect_
        real(dp) :: rmsd_err

        scale_ = optval(scale, .true.)
        ! Dimension checks
        d = size(P,dim=1)
        N = size(P,dim=2)
        if(any(shape(P)/=shape(Q)) .or. any(shape(R)/=[d,d]) .or. size(t)/=d) then
            call error_stop("array sizes do not match")
        end if
        if (present(W)) then
            if (size(W) /= N) then
                call error_stop("array sizes do not match")
            end if
        end if

        if(present(W)) then
            sum_w = stdlib_sum_kahan(W)
        else
            sum_w = real(N, kind = dp)
        end if
        !> leave opportunity for future discussion on how to add spmd reduction needed here to reduce sum_w before division

        sum_w = one_dp / sum_w
        if(sum_w<zero_dp) call error_stop("Invalid weights: sum of weights must be positive")

        allocate(c_P(d), c_Q(d), tmp_N(N), source=zero_dp)

        ! Compute centroids of P and Q
        if(present(W)) then
            do i = 1, d
                tmp_N(:) = W(:) * P(i,:)
                c_P(i) = stdlib_sum_kahan(tmp_N)
                tmp_N(:) = W(:) * Q(i,:)
                c_Q(i) = stdlib_sum_kahan(tmp_N)
            end do
        else
            c_P = stdlib_sum_kahan(P, dim=2)
            c_Q = stdlib_sum_kahan(Q, dim=2)
        end if
        c_P = c_P * sum_w
        c_Q = c_Q * sum_w

        ! Compute covariance matrix H = (P - c_P) * (Q - c_Q)^T and variance of P
        allocate(covariance(d,d), source=zero_dp)
        allocate(tmp_d(d), source=zero_dp)
        variance_p = zero_dp

        if (present(W)) then
            do point = 1, N
                tmp_d = P(:, point) - c_P(:)
                tmp_N(point) = stdlib_dot_product_kahan(tmp_d, tmp_d)
            end do
            tmp_N(:) = W(:) * tmp_N(:)
            variance_p = stdlib_sum_kahan(tmp_N)
            do j = 1, d
                do i = 1, d
                    tmp_N(:) = W(:) * (P(i,:) - c_P(i)) * (Q(j,:) - c_Q(j))
                    covariance(i,j) = stdlib_sum_kahan(tmp_N)
                end do
            end do
        else
            ! Calculate variance by the formula (1/n)*sigma(P - c_P)^2
            do point = 1, N
                tmp_d = P(:, point) - c_P(:)
                tmp_N(point) = stdlib_dot_product_kahan(tmp_d, tmp_d)
            end do
            variance_p = stdlib_sum_kahan(tmp_N)
            do j = 1, d
                do i = 1, d
                    covariance(i,j) = stdlib_dot_product_kahan((P(i,:) - c_P(i)),(Q(j,:) - c_Q(j)))
                end do
            end do
        end if

        covariance = covariance * sum_w
        variance_p = variance_p * sum_w

        allocate(U(d,d), source=zero_dp)
        allocate(Vt(d,d), source=zero_dp)
        allocate(S(d), source=zero_dp)

        ! SVD of covariance matrix H -> H = U * S * Vt
        call svd(covariance, S, U, Vt)

        ! Check for reflections in case of real entries.
        reflect_ = det(matmul(U,Vt)) < zero_dp
        if(reflect_) Vt(d,:) = -Vt(d,:)

        ! Optimal rotation matrix.
        call gemm(transa='N', transb='N', m=d,n=d,k=d, alpha=one_dp, a=U,lda=d, b=Vt, ldb=d, beta=zero_dp, c=R, ldc=d)

        ! Scaling factor
        c = one_dp
        if(scale_) then
            if(reflect_) then
                c = sum(S(1:d-1)) - S(d)
            else
                c = sum(S(1:d))
            end if
            c = variance_p / c
        end if

        ! Translation vector t = c_P - c*R*c_Q
        t = c_P
        call gemv(trans='N', m=d, n=d, alpha=-c, A=R, lda=d, x=c_Q, incx=1, beta=one_dp, y=t, incy=1)

        ! Compute RMSD
        allocate(vec(d), source=zero_dp)
        rmsd = zero_dp
        rmsd_err = zero_dp
        do point = 1, N
            ! Calculate the k^th difference vector by the formula vec_k = c*R*Q_k + t - P_k
            vec = t
            call gemv(trans='N', m=d, n=d, alpha=c, A=R, lda=d, x=Q(:, point), incx=1, beta=one_dp, y=vec, incy=1)
            vec = vec - P(:,point)
            temp = stdlib_dot_product_kahan(vec,vec)
            call kahan_kernel(temp, rmsd, rmsd_err)
        end do
        rmsd = sqrt(rmsd * sum_w)
    end subroutine
    module subroutine kabsch_umeyama_csp(P, Q, R, t, c, rmsd, W, scale)
        complex(sp), intent(in) :: P(:, :)
        !! Target point set (d × N)
        complex(sp), intent(in) :: Q(:, :)
        !! Reference point set (d × N)
        complex(sp), intent(out) :: R(:, :)
        !! Optimal rotation matrix (d × d)
        complex(sp), intent(out) :: t(:)
        !! Translation vector (d)
        complex(sp), intent(out) :: c
        !! Scale factor
        real(sp), intent(out) :: rmsd
        !! Root-mean-square deviation
        real(sp), intent(in), optional :: W(:)
        !! Optional weights
        logical, intent(in), optional :: scale
        !! Enable scaling

        ! Internal variables.
        integer :: i, j, point, d, N
        complex(sp), allocatable :: covariance(:,:), U(:,:), Vt(:,:), vec(:), tmp_N(:), tmp_d(:), c_P(:), c_Q(:)
        real(sp) :: sum_w, variance_p
        real(sp), allocatable :: S(:)
        complex(sp) :: temp
        real(sp) :: rtemp
        logical :: scale_
        real(sp) :: rmsd_err

        scale_ = optval(scale, .true.)
        ! Dimension checks
        d = size(P,dim=1)
        N = size(P,dim=2)
        if(any(shape(P)/=shape(Q)) .or. any(shape(R)/=[d,d]) .or. size(t)/=d) then
            call error_stop("array sizes do not match")
        end if
        if (present(W)) then
            if (size(W) /= N) then
                call error_stop("array sizes do not match")
            end if
        end if

        if(present(W)) then
            sum_w = stdlib_sum_kahan(W)
        else
            sum_w = real(N, kind = sp)
        end if
        !> leave opportunity for future discussion on how to add spmd reduction needed here to reduce sum_w before division

        sum_w = one_sp / sum_w
        if(sum_w<zero_sp) call error_stop("Invalid weights: sum of weights must be positive")

        allocate(c_P(d), c_Q(d), tmp_N(N), source=zero_csp)

        ! Compute centroids of P and Q
        if(present(W)) then
            do i = 1, d
                tmp_N(:) = W(:) * P(i,:)
                c_P(i) = stdlib_sum_kahan(tmp_N)
                tmp_N(:) = W(:) * Q(i,:)
                c_Q(i) = stdlib_sum_kahan(tmp_N)
            end do
        else
            c_P = stdlib_sum_kahan(P, dim=2)
            c_Q = stdlib_sum_kahan(Q, dim=2)
        end if
        c_P = c_P * sum_w
        c_Q = c_Q * sum_w

        ! Compute covariance matrix H = (P - c_P) * (Q - c_Q)^T and variance of P
        allocate(covariance(d,d), source=zero_csp)
        allocate(tmp_d(d), source=zero_csp)
        variance_p = zero_sp

        if (present(W)) then
            do point = 1, N
                tmp_d = P(:, point) - c_P(:)
                tmp_N(point) = stdlib_dot_product_kahan(tmp_d, tmp_d)
            end do
            tmp_N(:) = W(:) * tmp_N(:)
            variance_p = stdlib_sum_kahan(tmp_N)
            do j = 1, d
                do i = 1, d
                    tmp_N(:) = W(:) * (P(i,:) - c_P(i)) * conjg(Q(j,:) - c_Q(j))
                    covariance(i,j) = stdlib_sum_kahan(tmp_N)
                end do
            end do
        else
            ! Calculate variance by the formula (1/n)*sigma(P - c_P)^2
            do point = 1, N
                tmp_d = P(:, point) - c_P(:)
                tmp_N(point) = stdlib_dot_product_kahan(tmp_d, tmp_d)
            end do
            variance_p = stdlib_sum_kahan(tmp_N)
            do j = 1, d
                do i = 1, d
                    covariance(i,j) = stdlib_dot_product_kahan((Q(j,:) - c_Q(j)), (P(i,:) - c_P(i)))
                end do
            end do
        end if

        covariance = covariance * sum_w
        variance_p = variance_p * sum_w

        allocate(U(d,d), source=zero_csp)
        allocate(Vt(d,d), source=zero_csp)
        allocate(S(d), source=zero_sp)

        ! SVD of covariance matrix H -> H = U * S * Vt
        call svd(covariance, S, U, Vt)

        ! Check for reflections in case of real entries.

        ! Optimal rotation matrix.
        call gemm(transa='N', transb='N', m=d,n=d,k=d, alpha=one_csp, a=U,lda=d, b=Vt, ldb=d, beta=zero_csp, c=R, ldc=d)

        ! Scaling factor
        c = one_csp
        if(scale_) then
            c = sum(S(1:d))
            c = variance_p / c
        end if

        ! Translation vector t = c_P - c*R*c_Q
        t = c_P
        call gemv(trans='N', m=d, n=d, alpha=-c, A=R, lda=d, x=c_Q, incx=1, beta=one_csp, y=t, incy=1)

        ! Compute RMSD
        allocate(vec(d), source=zero_csp)
        rmsd = zero_sp
        rmsd_err = zero_sp
        do point = 1, N
            ! Calculate the k^th difference vector by the formula vec_k = c*R*Q_k + t - P_k
            vec = t
            call gemv(trans='N', m=d, n=d, alpha=c, A=R, lda=d, x=Q(:, point), incx=1, beta=one_csp, y=vec, incy=1)
            vec = vec - P(:,point)
            temp = stdlib_dot_product_kahan(vec,vec)
            rtemp = real(temp, kind=sp)
            call kahan_kernel(rtemp, rmsd, rmsd_err)
        end do
        rmsd = sqrt(rmsd * sum_w)
    end subroutine
    module subroutine kabsch_umeyama_cdp(P, Q, R, t, c, rmsd, W, scale)
        complex(dp), intent(in) :: P(:, :)
        !! Target point set (d × N)
        complex(dp), intent(in) :: Q(:, :)
        !! Reference point set (d × N)
        complex(dp), intent(out) :: R(:, :)
        !! Optimal rotation matrix (d × d)
        complex(dp), intent(out) :: t(:)
        !! Translation vector (d)
        complex(dp), intent(out) :: c
        !! Scale factor
        real(dp), intent(out) :: rmsd
        !! Root-mean-square deviation
        real(dp), intent(in), optional :: W(:)
        !! Optional weights
        logical, intent(in), optional :: scale
        !! Enable scaling

        ! Internal variables.
        integer :: i, j, point, d, N
        complex(dp), allocatable :: covariance(:,:), U(:,:), Vt(:,:), vec(:), tmp_N(:), tmp_d(:), c_P(:), c_Q(:)
        real(dp) :: sum_w, variance_p
        real(dp), allocatable :: S(:)
        complex(dp) :: temp
        real(dp) :: rtemp
        logical :: scale_
        real(dp) :: rmsd_err

        scale_ = optval(scale, .true.)
        ! Dimension checks
        d = size(P,dim=1)
        N = size(P,dim=2)
        if(any(shape(P)/=shape(Q)) .or. any(shape(R)/=[d,d]) .or. size(t)/=d) then
            call error_stop("array sizes do not match")
        end if
        if (present(W)) then
            if (size(W) /= N) then
                call error_stop("array sizes do not match")
            end if
        end if

        if(present(W)) then
            sum_w = stdlib_sum_kahan(W)
        else
            sum_w = real(N, kind = dp)
        end if
        !> leave opportunity for future discussion on how to add spmd reduction needed here to reduce sum_w before division

        sum_w = one_dp / sum_w
        if(sum_w<zero_dp) call error_stop("Invalid weights: sum of weights must be positive")

        allocate(c_P(d), c_Q(d), tmp_N(N), source=zero_cdp)

        ! Compute centroids of P and Q
        if(present(W)) then
            do i = 1, d
                tmp_N(:) = W(:) * P(i,:)
                c_P(i) = stdlib_sum_kahan(tmp_N)
                tmp_N(:) = W(:) * Q(i,:)
                c_Q(i) = stdlib_sum_kahan(tmp_N)
            end do
        else
            c_P = stdlib_sum_kahan(P, dim=2)
            c_Q = stdlib_sum_kahan(Q, dim=2)
        end if
        c_P = c_P * sum_w
        c_Q = c_Q * sum_w

        ! Compute covariance matrix H = (P - c_P) * (Q - c_Q)^T and variance of P
        allocate(covariance(d,d), source=zero_cdp)
        allocate(tmp_d(d), source=zero_cdp)
        variance_p = zero_dp

        if (present(W)) then
            do point = 1, N
                tmp_d = P(:, point) - c_P(:)
                tmp_N(point) = stdlib_dot_product_kahan(tmp_d, tmp_d)
            end do
            tmp_N(:) = W(:) * tmp_N(:)
            variance_p = stdlib_sum_kahan(tmp_N)
            do j = 1, d
                do i = 1, d
                    tmp_N(:) = W(:) * (P(i,:) - c_P(i)) * conjg(Q(j,:) - c_Q(j))
                    covariance(i,j) = stdlib_sum_kahan(tmp_N)
                end do
            end do
        else
            ! Calculate variance by the formula (1/n)*sigma(P - c_P)^2
            do point = 1, N
                tmp_d = P(:, point) - c_P(:)
                tmp_N(point) = stdlib_dot_product_kahan(tmp_d, tmp_d)
            end do
            variance_p = stdlib_sum_kahan(tmp_N)
            do j = 1, d
                do i = 1, d
                    covariance(i,j) = stdlib_dot_product_kahan((Q(j,:) - c_Q(j)), (P(i,:) - c_P(i)))
                end do
            end do
        end if

        covariance = covariance * sum_w
        variance_p = variance_p * sum_w

        allocate(U(d,d), source=zero_cdp)
        allocate(Vt(d,d), source=zero_cdp)
        allocate(S(d), source=zero_dp)

        ! SVD of covariance matrix H -> H = U * S * Vt
        call svd(covariance, S, U, Vt)

        ! Check for reflections in case of real entries.

        ! Optimal rotation matrix.
        call gemm(transa='N', transb='N', m=d,n=d,k=d, alpha=one_cdp, a=U,lda=d, b=Vt, ldb=d, beta=zero_cdp, c=R, ldc=d)

        ! Scaling factor
        c = one_cdp
        if(scale_) then
            c = sum(S(1:d))
            c = variance_p / c
        end if

        ! Translation vector t = c_P - c*R*c_Q
        t = c_P
        call gemv(trans='N', m=d, n=d, alpha=-c, A=R, lda=d, x=c_Q, incx=1, beta=one_cdp, y=t, incy=1)

        ! Compute RMSD
        allocate(vec(d), source=zero_cdp)
        rmsd = zero_dp
        rmsd_err = zero_dp
        do point = 1, N
            ! Calculate the k^th difference vector by the formula vec_k = c*R*Q_k + t - P_k
            vec = t
            call gemv(trans='N', m=d, n=d, alpha=c, A=R, lda=d, x=Q(:, point), incx=1, beta=one_cdp, y=vec, incy=1)
            vec = vec - P(:,point)
            temp = stdlib_dot_product_kahan(vec,vec)
            rtemp = real(temp, kind=dp)
            call kahan_kernel(rtemp, rmsd, rmsd_err)
        end do
        rmsd = sqrt(rmsd * sum_w)
    end subroutine
end submodule stdlib_spatial_kabsch_umeyama