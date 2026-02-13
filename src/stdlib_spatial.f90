module stdlib_spatial
    use stdlib_linalg_constants
    use stdlib_constants
    use stdlib_error, only: error_stop
    implicit none
    private
    public :: kabsch_umeyama

    interface kabsch_umeyama
        !-----------------------------------------------------------------------
        !> Compute the optimal similarity transform (Kabsch–Umeyama):
        !>
        !>     P ≈ c * R * Q + t
        !>
        !> where:
        !>   - R is an orthogonal rotation matrix
        !>   - c is an optional scale factor
        !>   - t is a translation vector
        !>
        !> The transformation minimizes the RMSD between corresponding columns
        !> of P and Q, optionally using weights.
        !-----------------------------------------------------------------------
        module subroutine kabsch_umeyama_sp(P, Q, R, t, c, rmsd, W, scale)
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
        end subroutine
        module subroutine kabsch_umeyama_dp(P, Q, R, t, c, rmsd, W, scale)
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
        end subroutine
        module subroutine kabsch_umeyama_csp(P, Q, R, t, c, rmsd, W, scale)
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
        end subroutine
        module subroutine kabsch_umeyama_cdp(P, Q, R, t, c, rmsd, W, scale)
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
        end subroutine
    end interface
end module stdlib_spatial