submodule(stdlib_spatial) stdlib_spatial_kabsch
    use stdlib_linalg, only: svd, det, eye
    use stdlib_intrinsics, only: stdlib_sum_kahan, stdlib_dot_product_kahan, stdlib_matmul, stdlib_sum, kahan_kernel

contains
    module subroutine kabsch_sp(P, Q, R, t, c, rmsd, W, scale)
        !> Reference point set (d × N)
        real(sp), intent(in) :: P(:, :)
        !> Target point set (d × N)
        real(sp), intent(in) :: Q(:, :)
        !> Optimal rotation matrix (d × d)
        real(sp), intent(out) :: R(:,:)
        !> Translation vector (d)
        real(sp), intent(out) :: t(:)
        !> Scale factor
        real(sp), intent(out) :: c
        !> Root-mean-square deviation
        real(sp), intent(out) :: rmsd
        !> Optional weights
        real(sp), intent(in), optional :: W(:)
        !> Enable scaling
        logical, intent(in), optional :: scale

        ! Internal variables.
        integer(ilp) :: i, j, point, d, N
        real(sp), allocatable :: covariance(:,:), U(:,:), Vt(:,:), B(:,:), vec(:), tmp_N(:), tmp_d(:), c_P(:), c_Q(:)
        real(sp) ::  vp, vq
        real(sp) :: sum_w, variance_p
        real(sp), allocatable :: S(:)
        logical :: scale_
        real(sp) :: rmsd_err


        ! Dimension checks
        if(size(P,dim=1)/=size(Q,dim=1) .or. size(P,dim=1)/=size(R,dim=1) .or. size(P,dim=1)/=size(R,dim=2) &
                    .or. size(P,dim=1)/=size(t)) then
            call error_stop("array sizes do not match")
        end if
        if(size(P,dim=2)/=size(Q,dim=2)) then
            call error_stop("array sizes do not match")
        end if
        if (present(W)) then
            if (size(W) /= size(P,dim=2)) then
                call error_stop("array sizes do not match")
            end if
        end if
        d = size(P,dim=1)
        N = size(P,dim=2)
        scale_ = .true.
        if(present(scale)) scale_ = scale

        sum_w = one_sp / N
        if(present(W)) sum_w = one_sp / stdlib_sum_kahan(W)

        allocate(c_P(d), source=zero_sp)
        allocate(c_Q(d), source=zero_sp)

        ! Compute centroids of P and Q
        if(present(W)) then
            do i = 1, d
                c_P(i) = stdlib_dot_product_kahan(w,P(i, :))
                c_Q(i) = stdlib_dot_product_kahan(w,Q(i, :))
            end do
        else
            do i = 1, d
                c_P(i) = stdlib_sum_kahan(P(i, :))
                c_Q(i) = stdlib_sum_kahan(Q(i, :))
            end do
        end if
        c_P = c_P * sum_w
        c_Q = c_Q * sum_w

        ! Compute covariance matrix H = (P - c_P) * (Q - c_Q)^T and variance of P
        allocate(covariance(d,d), source=zero_sp)
        allocate(tmp_N(N), source=zero_sp)
        allocate(tmp_d(d), source=zero_sp)
        variance_p = zero_sp

        if (present(W)) then
            do point = 1, N
                tmp_d = P(:, point) - c_P(:)
                tmp_N(point) = stdlib_dot_product_kahan(tmp_d, tmp_d)
            end do
            variance_p = stdlib_dot_product_kahan(w, tmp_N)
            do j = 1, d
                do i = 1, d
                    tmp_N(:) = (P(i,:) - c_P(i)) * (Q(j,:) - c_Q(j))
                    covariance(i,j) = stdlib_dot_product_kahan(w, tmp_N)
                end do
            end do
        else
            do point = 1, N
                tmp_d = P(:, point) - c_P(:)
                tmp_N(point) = stdlib_dot_product_kahan(tmp_d, tmp_d)
            end do
            variance_p = stdlib_sum_kahan(tmp_N)
            do j = 1, d
                do i = 1, d
                    tmp_N(:) = (P(i,:) - c_P(i)) * (Q(j,:) - c_Q(j))
                    covariance(i,j) = stdlib_sum_kahan(tmp_N)
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

        ! Optimal rotation matrix.
        do i = 1,d
            do j = 1,d
                R(i,j) = stdlib_dot_product_kahan(U(i,:), Vt(:, j))
            end do
        end do

        ! Scaling factor
        c = variance_p / (sum(S(1:d)))
        if (.not. scale_) c = one_sp

        ! Translation vector t = c_P - c*R*c_Q
        do i = 1, d
            t(i) = c_P(i) - c * stdlib_dot_product_kahan(R(i,1:d), c_Q(1:d))
        end do

        ! Compute RMSD
        allocate(vec(d), source=zero_sp)
        rmsd = zero_sp
        rmsd_err = zero_sp
        do point = 1, N
            ! Calculate the k^th difference vector by the formula vec_k = c*R*Q_k + t - P_k
            do i = 1, d
                vec(i) = c * stdlib_dot_product_kahan(R(i,1:d), Q(1:d,point))
            end do
            vec(1:d) = vec(1:d) + t(1:d) - P(1:d,point)
            call kahan_kernel(real(stdlib_dot_product_kahan(vec,vec), kind=sp), rmsd, rmsd_err)
        end do
        rmsd = sqrt(rmsd * sum_w)
    end subroutine
    module subroutine kabsch_dp(P, Q, R, t, c, rmsd, W, scale)
        !> Reference point set (d × N)
        real(dp), intent(in) :: P(:, :)
        !> Target point set (d × N)
        real(dp), intent(in) :: Q(:, :)
        !> Optimal rotation matrix (d × d)
        real(dp), intent(out) :: R(:,:)
        !> Translation vector (d)
        real(dp), intent(out) :: t(:)
        !> Scale factor
        real(dp), intent(out) :: c
        !> Root-mean-square deviation
        real(dp), intent(out) :: rmsd
        !> Optional weights
        real(dp), intent(in), optional :: W(:)
        !> Enable scaling
        logical, intent(in), optional :: scale

        ! Internal variables.
        integer(ilp) :: i, j, point, d, N
        real(dp), allocatable :: covariance(:,:), U(:,:), Vt(:,:), B(:,:), vec(:), tmp_N(:), tmp_d(:), c_P(:), c_Q(:)
        real(dp) ::  vp, vq
        real(dp) :: sum_w, variance_p
        real(dp), allocatable :: S(:)
        logical :: scale_
        real(dp) :: rmsd_err


        ! Dimension checks
        if(size(P,dim=1)/=size(Q,dim=1) .or. size(P,dim=1)/=size(R,dim=1) .or. size(P,dim=1)/=size(R,dim=2) &
                    .or. size(P,dim=1)/=size(t)) then
            call error_stop("array sizes do not match")
        end if
        if(size(P,dim=2)/=size(Q,dim=2)) then
            call error_stop("array sizes do not match")
        end if
        if (present(W)) then
            if (size(W) /= size(P,dim=2)) then
                call error_stop("array sizes do not match")
            end if
        end if
        d = size(P,dim=1)
        N = size(P,dim=2)
        scale_ = .true.
        if(present(scale)) scale_ = scale

        sum_w = one_dp / N
        if(present(W)) sum_w = one_dp / stdlib_sum_kahan(W)

        allocate(c_P(d), source=zero_dp)
        allocate(c_Q(d), source=zero_dp)

        ! Compute centroids of P and Q
        if(present(W)) then
            do i = 1, d
                c_P(i) = stdlib_dot_product_kahan(w,P(i, :))
                c_Q(i) = stdlib_dot_product_kahan(w,Q(i, :))
            end do
        else
            do i = 1, d
                c_P(i) = stdlib_sum_kahan(P(i, :))
                c_Q(i) = stdlib_sum_kahan(Q(i, :))
            end do
        end if
        c_P = c_P * sum_w
        c_Q = c_Q * sum_w

        ! Compute covariance matrix H = (P - c_P) * (Q - c_Q)^T and variance of P
        allocate(covariance(d,d), source=zero_dp)
        allocate(tmp_N(N), source=zero_dp)
        allocate(tmp_d(d), source=zero_dp)
        variance_p = zero_dp

        if (present(W)) then
            do point = 1, N
                tmp_d = P(:, point) - c_P(:)
                tmp_N(point) = stdlib_dot_product_kahan(tmp_d, tmp_d)
            end do
            variance_p = stdlib_dot_product_kahan(w, tmp_N)
            do j = 1, d
                do i = 1, d
                    tmp_N(:) = (P(i,:) - c_P(i)) * (Q(j,:) - c_Q(j))
                    covariance(i,j) = stdlib_dot_product_kahan(w, tmp_N)
                end do
            end do
        else
            do point = 1, N
                tmp_d = P(:, point) - c_P(:)
                tmp_N(point) = stdlib_dot_product_kahan(tmp_d, tmp_d)
            end do
            variance_p = stdlib_sum_kahan(tmp_N)
            do j = 1, d
                do i = 1, d
                    tmp_N(:) = (P(i,:) - c_P(i)) * (Q(j,:) - c_Q(j))
                    covariance(i,j) = stdlib_sum_kahan(tmp_N)
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

        ! Optimal rotation matrix.
        do i = 1,d
            do j = 1,d
                R(i,j) = stdlib_dot_product_kahan(U(i,:), Vt(:, j))
            end do
        end do

        ! Scaling factor
        c = variance_p / (sum(S(1:d)))
        if (.not. scale_) c = one_dp

        ! Translation vector t = c_P - c*R*c_Q
        do i = 1, d
            t(i) = c_P(i) - c * stdlib_dot_product_kahan(R(i,1:d), c_Q(1:d))
        end do

        ! Compute RMSD
        allocate(vec(d), source=zero_dp)
        rmsd = zero_dp
        rmsd_err = zero_dp
        do point = 1, N
            ! Calculate the k^th difference vector by the formula vec_k = c*R*Q_k + t - P_k
            do i = 1, d
                vec(i) = c * stdlib_dot_product_kahan(R(i,1:d), Q(1:d,point))
            end do
            vec(1:d) = vec(1:d) + t(1:d) - P(1:d,point)
            call kahan_kernel(real(stdlib_dot_product_kahan(vec,vec), kind=dp), rmsd, rmsd_err)
        end do
        rmsd = sqrt(rmsd * sum_w)
    end subroutine
    module subroutine kabsch_csp(P, Q, R, t, c, rmsd, W, scale)
        !> Reference point set (d × N)
        complex(sp), intent(in) :: P(:, :)
        !> Target point set (d × N)
        complex(sp), intent(in) :: Q(:, :)
        !> Optimal rotation matrix (d × d)
        complex(sp), intent(out) :: R(:,:)
        !> Translation vector (d)
        complex(sp), intent(out) :: t(:)
        !> Scale factor
        complex(sp), intent(out) :: c
        !> Root-mean-square deviation
        real(sp), intent(out) :: rmsd
        !> Optional weights
        complex(sp), intent(in), optional :: W(:)
        !> Enable scaling
        logical, intent(in), optional :: scale

        ! Internal variables.
        integer(ilp) :: i, j, point, d, N
        complex(sp), allocatable :: covariance(:,:), U(:,:), Vt(:,:), B(:,:), vec(:), tmp_N(:), tmp_d(:), c_P(:), c_Q(:)
        complex(sp) ::  vp, vq
        real(sp) :: sum_w, variance_p
        real(sp), allocatable :: S(:)
        logical :: scale_
        real(sp) :: rmsd_err


        ! Dimension checks
        if(size(P,dim=1)/=size(Q,dim=1) .or. size(P,dim=1)/=size(R,dim=1) .or. size(P,dim=1)/=size(R,dim=2) &
                    .or. size(P,dim=1)/=size(t)) then
            call error_stop("array sizes do not match")
        end if
        if(size(P,dim=2)/=size(Q,dim=2)) then
            call error_stop("array sizes do not match")
        end if
        if (present(W)) then
            if (size(W) /= size(P,dim=2)) then
                call error_stop("array sizes do not match")
            end if
        end if
        d = size(P,dim=1)
        N = size(P,dim=2)
        scale_ = .true.
        if(present(scale)) scale_ = scale

        sum_w = one_csp / N
        if(present(W)) sum_w = one_csp / stdlib_sum_kahan(W)

        allocate(c_P(d), source=zero_csp)
        allocate(c_Q(d), source=zero_csp)

        ! Compute centroids of P and Q
        if(present(W)) then
            do i = 1, d
                c_P(i) = stdlib_dot_product_kahan(w,P(i, :))
                c_Q(i) = stdlib_dot_product_kahan(w,Q(i, :))
            end do
        else
            do i = 1, d
                c_P(i) = stdlib_sum_kahan(P(i, :))
                c_Q(i) = stdlib_sum_kahan(Q(i, :))
            end do
        end if
        c_P = c_P * sum_w
        c_Q = c_Q * sum_w

        ! Compute covariance matrix H = (P - c_P) * (Q - c_Q)^T and variance of P
        allocate(covariance(d,d), source=zero_csp)
        allocate(tmp_N(N), source=zero_csp)
        allocate(tmp_d(d), source=zero_csp)
        variance_p = zero_csp

        if (present(W)) then
            do point = 1, N
                tmp_d = P(:, point) - c_P(:)
                tmp_N(point) = stdlib_dot_product_kahan(tmp_d, tmp_d)
            end do
            variance_p = stdlib_dot_product_kahan(w, tmp_N)
            do j = 1, d
                do i = 1, d
                    tmp_N(:) = (P(i,:) - c_P(i)) * conjg(Q(j,:) - c_Q(j))
                    covariance(i,j) = stdlib_dot_product_kahan(w, tmp_N)
                end do
            end do
        else
            do point = 1, N
                tmp_d = P(:, point) - c_P(:)
                tmp_N(point) = stdlib_dot_product_kahan(tmp_d, tmp_d)
            end do
            variance_p = stdlib_sum_kahan(tmp_N)
            do j = 1, d
                do i = 1, d
                    tmp_N(:) = (P(i,:) - c_P(i)) * conjg(Q(j,:) - c_Q(j))
                    covariance(i,j) = stdlib_sum_kahan(tmp_N)
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

        ! Optimal rotation matrix.
        do i = 1,d
            do j = 1,d
                R(i,j) = stdlib_dot_product_kahan(conjg(U(i,:)), Vt(:, j))
            end do
        end do

        ! Scaling factor
        c = variance_p / (sum(S(1:d)))
        if (.not. scale_) c = one_csp

        ! Translation vector t = c_P - c*R*c_Q
        do i = 1, d
            t(i) = c_P(i) - c * stdlib_dot_product_kahan(conjg(R(i,1:d)), c_Q(1:d))
        end do

        ! Compute RMSD
        allocate(vec(d), source=zero_csp)
        rmsd = zero_csp
        rmsd_err = zero_csp
        do point = 1, N
            ! Calculate the k^th difference vector by the formula vec_k = c*R*Q_k + t - P_k
            do i = 1, d
                vec(i) = c * stdlib_dot_product_kahan(conjg(R(i,1:d)), Q(1:d,point))
            end do
            vec(1:d) = vec(1:d) + t(1:d) - P(1:d,point)
            call kahan_kernel(real(stdlib_dot_product_kahan(vec,vec), kind=sp), rmsd, rmsd_err)
        end do
        rmsd = sqrt(rmsd * sum_w)
    end subroutine
    module subroutine kabsch_cdp(P, Q, R, t, c, rmsd, W, scale)
        !> Reference point set (d × N)
        complex(dp), intent(in) :: P(:, :)
        !> Target point set (d × N)
        complex(dp), intent(in) :: Q(:, :)
        !> Optimal rotation matrix (d × d)
        complex(dp), intent(out) :: R(:,:)
        !> Translation vector (d)
        complex(dp), intent(out) :: t(:)
        !> Scale factor
        complex(dp), intent(out) :: c
        !> Root-mean-square deviation
        real(dp), intent(out) :: rmsd
        !> Optional weights
        complex(dp), intent(in), optional :: W(:)
        !> Enable scaling
        logical, intent(in), optional :: scale

        ! Internal variables.
        integer(ilp) :: i, j, point, d, N
        complex(dp), allocatable :: covariance(:,:), U(:,:), Vt(:,:), B(:,:), vec(:), tmp_N(:), tmp_d(:), c_P(:), c_Q(:)
        complex(dp) ::  vp, vq
        real(dp) :: sum_w, variance_p
        real(dp), allocatable :: S(:)
        logical :: scale_
        real(dp) :: rmsd_err


        ! Dimension checks
        if(size(P,dim=1)/=size(Q,dim=1) .or. size(P,dim=1)/=size(R,dim=1) .or. size(P,dim=1)/=size(R,dim=2) &
                    .or. size(P,dim=1)/=size(t)) then
            call error_stop("array sizes do not match")
        end if
        if(size(P,dim=2)/=size(Q,dim=2)) then
            call error_stop("array sizes do not match")
        end if
        if (present(W)) then
            if (size(W) /= size(P,dim=2)) then
                call error_stop("array sizes do not match")
            end if
        end if
        d = size(P,dim=1)
        N = size(P,dim=2)
        scale_ = .true.
        if(present(scale)) scale_ = scale

        sum_w = one_cdp / N
        if(present(W)) sum_w = one_cdp / stdlib_sum_kahan(W)

        allocate(c_P(d), source=zero_cdp)
        allocate(c_Q(d), source=zero_cdp)

        ! Compute centroids of P and Q
        if(present(W)) then
            do i = 1, d
                c_P(i) = stdlib_dot_product_kahan(w,P(i, :))
                c_Q(i) = stdlib_dot_product_kahan(w,Q(i, :))
            end do
        else
            do i = 1, d
                c_P(i) = stdlib_sum_kahan(P(i, :))
                c_Q(i) = stdlib_sum_kahan(Q(i, :))
            end do
        end if
        c_P = c_P * sum_w
        c_Q = c_Q * sum_w

        ! Compute covariance matrix H = (P - c_P) * (Q - c_Q)^T and variance of P
        allocate(covariance(d,d), source=zero_cdp)
        allocate(tmp_N(N), source=zero_cdp)
        allocate(tmp_d(d), source=zero_cdp)
        variance_p = zero_cdp

        if (present(W)) then
            do point = 1, N
                tmp_d = P(:, point) - c_P(:)
                tmp_N(point) = stdlib_dot_product_kahan(tmp_d, tmp_d)
            end do
            variance_p = stdlib_dot_product_kahan(w, tmp_N)
            do j = 1, d
                do i = 1, d
                    tmp_N(:) = (P(i,:) - c_P(i)) * conjg(Q(j,:) - c_Q(j))
                    covariance(i,j) = stdlib_dot_product_kahan(w, tmp_N)
                end do
            end do
        else
            do point = 1, N
                tmp_d = P(:, point) - c_P(:)
                tmp_N(point) = stdlib_dot_product_kahan(tmp_d, tmp_d)
            end do
            variance_p = stdlib_sum_kahan(tmp_N)
            do j = 1, d
                do i = 1, d
                    tmp_N(:) = (P(i,:) - c_P(i)) * conjg(Q(j,:) - c_Q(j))
                    covariance(i,j) = stdlib_sum_kahan(tmp_N)
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

        ! Optimal rotation matrix.
        do i = 1,d
            do j = 1,d
                R(i,j) = stdlib_dot_product_kahan(conjg(U(i,:)), Vt(:, j))
            end do
        end do

        ! Scaling factor
        c = variance_p / (sum(S(1:d)))
        if (.not. scale_) c = one_cdp

        ! Translation vector t = c_P - c*R*c_Q
        do i = 1, d
            t(i) = c_P(i) - c * stdlib_dot_product_kahan(conjg(R(i,1:d)), c_Q(1:d))
        end do

        ! Compute RMSD
        allocate(vec(d), source=zero_cdp)
        rmsd = zero_cdp
        rmsd_err = zero_cdp
        do point = 1, N
            ! Calculate the k^th difference vector by the formula vec_k = c*R*Q_k + t - P_k
            do i = 1, d
                vec(i) = c * stdlib_dot_product_kahan(conjg(R(i,1:d)), Q(1:d,point))
            end do
            vec(1:d) = vec(1:d) + t(1:d) - P(1:d,point)
            call kahan_kernel(real(stdlib_dot_product_kahan(vec,vec), kind=dp), rmsd, rmsd_err)
        end do
        rmsd = sqrt(rmsd * sum_w)
    end subroutine
end submodule stdlib_spatial_kabsch