module stdlib_spatial
    use stdlib_kinds, only: sp, dp, xdp, qp
    use stdlib_constants
    implicit none
    private
    public :: kabsch_umeyama

    interface kabsch_umeyama
        !! ([Specifications](../page/specs/stdlib_spatial.html#kabsch_umeyama))
        !! This interface computes the optimal similarity transform (Kabsch–Umeyama):
        !! \[
        !!      P \approx c \, R \, Q + t
        !! \]
        !! The transformation minimizes the RMSD between corresponding columns
        !! of P and Q, optionally using weights and with optional scaling.
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
        end subroutine
    end interface
end module stdlib_spatial