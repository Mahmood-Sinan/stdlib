submodule(stdlib_spatial) stdlib_spatial_kabsch

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
        real(sp), allocatable :: c_P(:), c_Q(:), covariance(:,:), U(:,:), Vt(:,:), B(:,:), vec(:)
        real(sp) ::  vp, vq
        real(sp) :: sum_w, variance_p
        real(sp), allocatable :: S(:)
        logical :: scale_

        ! Dimension checks
        if(size(P,dim=1)/=size(Q,dim=1) .or. size(P,dim=1)/=size(R,dim=1) .or. size(P,dim=1)/=size(R,dim=2) &
                    .or. size(P,dim=1)/=size(t)) then
            error stop "array sizes do not match"
        end if
        if(size(P,dim=2)/=size(Q,dim=2)) then
            error stop "array sizes do not match"
        end if
        if (present(W)) then
            if (size(W) /= size(P,dim=2)) then
                error stop "array sizes do not match"
            end if
        end if
        d = size(P,dim=1)
        N = size(P,dim=2)
        scale_ = .true.
        if(present(scale)) scale_ = scale

        sum_w = one_sp / N
        if(present(W)) sum_w = one_sp / stdlib_sum(W)

        allocate(c_P(d), source=zero_sp)
        allocate(c_Q(d), source=zero_sp)

        ! Compute centroids of P and Q
        if(present(W)) then
            do point = 1, N
                c_P(1:d) = c_P(1:d) + P(1:d,point)*W(point)
                c_Q(1:d) = c_Q(1:d) + Q(1:d,point)*W(point)
            end do
        else
            do point = 1, N
                c_P(1:d) = c_P(1:d) + P(1:d,point)
                c_Q(1:d) = c_Q(1:d) + Q(1:d,point)
            end do
        end if
        c_P = c_P * sum_w
        c_Q = c_Q * sum_w

        ! Compute covariance matrix H = (P - c_P) * (Q - c_Q)^T and variance of P
        allocate(covariance(d,d), source=zero_sp)
        variance_p = zero_sp

        if (present(W)) then
            do point = 1, N
                do i = 1, d
                    vp = P(i,point) - c_P(i)
                    variance_p = variance_p + vp*vp * W(point)
                end do
                do j = 1, d
                    vq = Q(j,point) - c_Q(j)
                    do i = 1, d
                        vp = P(i,point) - c_P(i)
                        covariance(i,j) = covariance(i,j) + vp*vq * W(point)
                    end do
                end do
            end do
        else
            do point = 1, N
                do i = 1, d
                    vp = P(i,point) - c_P(i)
                    variance_p = variance_p + vp*vp
                end do
                do j = 1, d
                    vq = Q(j,point) - c_Q(j)
                    do i = 1, d
                        vp = P(i,point) - c_P(i)
                        covariance(i,j) = covariance(i,j) + vp*vq
                    end do
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

        allocate(B(d,d), source=zero_sp)

        ! Validate right-handed coordinate system
        R = matmul(U, Vt)
        B = eye(d,d)
        B(d,d) = sign(one_sp, det(R))

        ! Optimal rotation matrix.
        R = matmul(U,matmul(B,Vt))

        ! Scaling factor
        c = variance_p / (sum(S(1:d-1))+B(d,d)*S(d))
        if (.not. scale_) c = one_sp

        ! Translation vector
        t = c_P - c*matmul(R , c_Q)

        ! Compute RMSD
        allocate(vec(d), source=zero_sp)
        rmsd = zero_sp
        do i = 1, N
            vec(1:d) = c * matmul(R, Q(1:d,i))
            vec(1:d) = vec(1:d) + t(1:d)
            vec(1:d) = vec(1:d) - P(1:d,i)
            rmsd = rmsd + dot_product(vec, vec)
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
        real(dp), allocatable :: c_P(:), c_Q(:), covariance(:,:), U(:,:), Vt(:,:), B(:,:), vec(:)
        real(dp) ::  vp, vq
        real(dp) :: sum_w, variance_p
        real(dp), allocatable :: S(:)
        logical :: scale_

        ! Dimension checks
        if(size(P,dim=1)/=size(Q,dim=1) .or. size(P,dim=1)/=size(R,dim=1) .or. size(P,dim=1)/=size(R,dim=2) &
                    .or. size(P,dim=1)/=size(t)) then
            error stop "array sizes do not match"
        end if
        if(size(P,dim=2)/=size(Q,dim=2)) then
            error stop "array sizes do not match"
        end if
        if (present(W)) then
            if (size(W) /= size(P,dim=2)) then
                error stop "array sizes do not match"
            end if
        end if
        d = size(P,dim=1)
        N = size(P,dim=2)
        scale_ = .true.
        if(present(scale)) scale_ = scale

        sum_w = one_dp / N
        if(present(W)) sum_w = one_dp / stdlib_sum(W)

        allocate(c_P(d), source=zero_dp)
        allocate(c_Q(d), source=zero_dp)

        ! Compute centroids of P and Q
        if(present(W)) then
            do point = 1, N
                c_P(1:d) = c_P(1:d) + P(1:d,point)*W(point)
                c_Q(1:d) = c_Q(1:d) + Q(1:d,point)*W(point)
            end do
        else
            do point = 1, N
                c_P(1:d) = c_P(1:d) + P(1:d,point)
                c_Q(1:d) = c_Q(1:d) + Q(1:d,point)
            end do
        end if
        c_P = c_P * sum_w
        c_Q = c_Q * sum_w

        ! Compute covariance matrix H = (P - c_P) * (Q - c_Q)^T and variance of P
        allocate(covariance(d,d), source=zero_dp)
        variance_p = zero_dp

        if (present(W)) then
            do point = 1, N
                do i = 1, d
                    vp = P(i,point) - c_P(i)
                    variance_p = variance_p + vp*vp * W(point)
                end do
                do j = 1, d
                    vq = Q(j,point) - c_Q(j)
                    do i = 1, d
                        vp = P(i,point) - c_P(i)
                        covariance(i,j) = covariance(i,j) + vp*vq * W(point)
                    end do
                end do
            end do
        else
            do point = 1, N
                do i = 1, d
                    vp = P(i,point) - c_P(i)
                    variance_p = variance_p + vp*vp
                end do
                do j = 1, d
                    vq = Q(j,point) - c_Q(j)
                    do i = 1, d
                        vp = P(i,point) - c_P(i)
                        covariance(i,j) = covariance(i,j) + vp*vq
                    end do
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

        allocate(B(d,d), source=zero_dp)

        ! Validate right-handed coordinate system
        R = matmul(U, Vt)
        B = eye(d,d)
        B(d,d) = sign(one_dp, det(R))

        ! Optimal rotation matrix.
        R = matmul(U,matmul(B,Vt))

        ! Scaling factor
        c = variance_p / (sum(S(1:d-1))+B(d,d)*S(d))
        if (.not. scale_) c = one_dp

        ! Translation vector
        t = c_P - c*matmul(R , c_Q)

        ! Compute RMSD
        allocate(vec(d), source=zero_dp)
        rmsd = zero_dp
        do i = 1, N
            vec(1:d) = c * matmul(R, Q(1:d,i))
            vec(1:d) = vec(1:d) + t(1:d)
            vec(1:d) = vec(1:d) - P(1:d,i)
            rmsd = rmsd + dot_product(vec, vec)
        end do
        rmsd = sqrt(rmsd * sum_w)
    end subroutine
end submodule stdlib_spatial_kabsch
