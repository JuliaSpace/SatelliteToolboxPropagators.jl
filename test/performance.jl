## Description #############################################################################
#
# Performance tests using AllocCheck.jl.
#
############################################################################################

@testset "Allocation Check" begin
    _D = ForwardDiff.Dual{ForwardDiff.Tag{Nothing, Float64}, Float64, 6}

    # The kernels and the Jacobians must not allocate on any Julia version. The fitting
    # functions are checked against a maximum number of allocation sites, which are all in
    # the verbose printing path. Julia 1.12 compiles that path with many more sites, so the
    # limits below only hold on the previous releases.
    check_fit_allocs = VERSION < v"1.12"
    check_fit_allocs ||
        @warn "The fitting allocation sites are only checked on Julia 1.10 and 1.11."

    # == Two-Body ======================================================================

    @testset "twobody_init!" begin
        @test length(
            check_allocs(
                (tbd, orb₀) -> twobody_init!(tbd, orb₀),
                (
                    TwoBodyPropagator{Float64, Float64},
                    KeplerianElements{TrueAnomaly, Float64, Float64},
                ),
            ),
        ) == 0
    end

    @testset "twobody!" begin
        @test length(
            check_allocs(
                (tbd, t) -> twobody!(tbd, t), (TwoBodyPropagator{Float64, Float64}, Float64)
            ),
        ) == 0
    end

    @testset "_mean_elements_jacobian (FiniteDiffJacobian)" begin
        @test length(
            check_allocs(
                (tbd, Δt, x₁, y₁) -> begin
                    SatelliteToolboxPropagators._mean_elements_jacobian(
                        FiniteDiffJacobian(), tbd, Δt, x₁, y₁
                    )
                end,
                (
                    TwoBodyPropagator{Float64, Float64},
                    Float64,
                    SVector{6, Float64},
                    SVector{6, Float64},
                ),
            ),
        ) == 0
    end

    @testset "_mean_elements_jacobian (ForwardDiffJacobian)" begin
        @test length(
            check_allocs(
                (tbd, tbd_ad, Δt, x₁, y₁) -> begin
                    SatelliteToolboxPropagators._mean_elements_jacobian(
                        ForwardDiffJacobian(), tbd, Δt, x₁, y₁; pd_ad = tbd_ad
                    )
                end,
                (
                    TwoBodyPropagator{Float64, Float64},
                    TwoBodyPropagator{Float64, _D},
                    Float64,
                    SVector{6, Float64},
                    SVector{6, Float64},
                ),
            ),
        ) == 0
    end

    check_fit_allocs && @testset "fit_twobody_mean_elements! (FiniteDiffJacobian)" begin
        @test length(
            check_allocs(
                (tbd, vjd, vr_i, vv_i) -> begin
                    fit_twobody_mean_elements!(
                        tbd,
                        vjd,
                        vr_i,
                        vv_i;
                        jacobian_method = FiniteDiffJacobian(),
                        verbose = false,
                    )
                end,
                (
                    TwoBodyPropagator{Float64, Float64},
                    Vector{Float64},
                    Vector{SVector{3, Float64}},
                    Vector{SVector{3, Float64}},
                ),
            ),
        ) <= 14
    end

    check_fit_allocs && @testset "fit_twobody_mean_elements! (ForwardDiffJacobian)" begin
        @test length(
            check_allocs(
                (tbd, vjd, vr_i, vv_i) -> begin
                    fit_twobody_mean_elements!(
                        tbd,
                        vjd,
                        vr_i,
                        vv_i;
                        jacobian_method = ForwardDiffJacobian(),
                        verbose = false,
                    )
                end,
                (
                    TwoBodyPropagator{Float64, Float64},
                    Vector{Float64},
                    Vector{SVector{3, Float64}},
                    Vector{SVector{3, Float64}},
                ),
            ),
        ) <= 14
    end

    # == J2 ============================================================================

    @testset "j2_init!" begin
        @test length(
            check_allocs(
                (j2d, orb₀) -> j2_init!(j2d, orb₀),
                (
                    J2Propagator{Float64, Float64},
                    KeplerianElements{TrueAnomaly, Float64, Float64},
                ),
            ),
        ) == 0
    end

    @testset "j2!" begin
        @test length(
            check_allocs((j2d, t) -> j2!(j2d, t), (J2Propagator{Float64, Float64}, Float64))
        ) == 0
    end

    @testset "_mean_elements_jacobian (FiniteDiffJacobian)" begin
        @test length(
            check_allocs(
                (j2d, Δt, x₁, y₁) -> begin
                    SatelliteToolboxPropagators._mean_elements_jacobian(
                        FiniteDiffJacobian(), j2d, Δt, x₁, y₁
                    )
                end,
                (
                    J2Propagator{Float64, Float64},
                    Float64,
                    SVector{6, Float64},
                    SVector{6, Float64},
                ),
            ),
        ) == 0
    end

    @testset "_mean_elements_jacobian (ForwardDiffJacobian)" begin
        @test length(
            check_allocs(
                (j2d, j2d_ad, Δt, x₁, y₁) -> begin
                    SatelliteToolboxPropagators._mean_elements_jacobian(
                        ForwardDiffJacobian(), j2d, Δt, x₁, y₁; pd_ad = j2d_ad
                    )
                end,
                (
                    J2Propagator{Float64, Float64},
                    J2Propagator{Float64, _D},
                    Float64,
                    SVector{6, Float64},
                    SVector{6, Float64},
                ),
            ),
        ) == 0
    end

    check_fit_allocs && @testset "fit_j2_mean_elements! (FiniteDiffJacobian)" begin
        @test length(
            check_allocs(
                (j2d, vjd, vr_i, vv_i) -> begin
                    fit_j2_mean_elements!(
                        j2d,
                        vjd,
                        vr_i,
                        vv_i;
                        jacobian_method = FiniteDiffJacobian(),
                        verbose = false,
                    )
                end,
                (
                    J2Propagator{Float64, Float64},
                    Vector{Float64},
                    Vector{SVector{3, Float64}},
                    Vector{SVector{3, Float64}},
                ),
            ),
        ) <= 14
    end

    check_fit_allocs && @testset "fit_j2_mean_elements! (ForwardDiffJacobian)" begin
        @test length(
            check_allocs(
                (j2d, vjd, vr_i, vv_i) -> begin
                    fit_j2_mean_elements!(
                        j2d,
                        vjd,
                        vr_i,
                        vv_i;
                        jacobian_method = ForwardDiffJacobian(),
                        verbose = false,
                    )
                end,
                (
                    J2Propagator{Float64, Float64},
                    Vector{Float64},
                    Vector{SVector{3, Float64}},
                    Vector{SVector{3, Float64}},
                ),
            ),
        ) <= 14
    end

    # == J2 Osculating =================================================================

    @testset "j2osc_init!" begin
        @test length(
            check_allocs(
                (j2oscd, orb₀) -> j2osc_init!(j2oscd, orb₀),
                (
                    J2OsculatingPropagator{Float64, Float64},
                    KeplerianElements{TrueAnomaly, Float64, Float64},
                ),
            ),
        ) == 0
    end

    @testset "j2osc!" begin
        @test length(
            check_allocs(
                (j2oscd, t) -> j2osc!(j2oscd, t),
                (J2OsculatingPropagator{Float64, Float64}, Float64),
            ),
        ) == 0
    end

    @testset "_mean_elements_jacobian (FiniteDiffJacobian)" begin
        @test length(
            check_allocs(
                (j2oscd, Δt, x₁, y₁) -> begin
                    SatelliteToolboxPropagators._mean_elements_jacobian(
                        FiniteDiffJacobian(), j2oscd, Δt, x₁, y₁
                    )
                end,
                (
                    J2OsculatingPropagator{Float64, Float64},
                    Float64,
                    SVector{6, Float64},
                    SVector{6, Float64},
                ),
            ),
        ) == 0
    end

    @testset "_mean_elements_jacobian (ForwardDiffJacobian)" begin
        @test length(
            check_allocs(
                (j2oscd, j2oscd_ad, Δt, x₁, y₁) -> begin
                    SatelliteToolboxPropagators._mean_elements_jacobian(
                        ForwardDiffJacobian(), j2oscd, Δt, x₁, y₁; pd_ad = j2oscd_ad
                    )
                end,
                (
                    J2OsculatingPropagator{Float64, Float64},
                    J2OsculatingPropagator{Float64, _D},
                    Float64,
                    SVector{6, Float64},
                    SVector{6, Float64},
                ),
            ),
        ) == 0
    end

    check_fit_allocs && @testset "fit_j2osc_mean_elements! (FiniteDiffJacobian)" begin
        @test length(
            check_allocs(
                (j2oscd, vjd, vr_i, vv_i) -> begin
                    fit_j2osc_mean_elements!(
                        j2oscd,
                        vjd,
                        vr_i,
                        vv_i;
                        jacobian_method = FiniteDiffJacobian(),
                        verbose = false,
                    )
                end,
                (
                    J2OsculatingPropagator{Float64, Float64},
                    Vector{Float64},
                    Vector{SVector{3, Float64}},
                    Vector{SVector{3, Float64}},
                ),
            ),
        ) <= 15
    end

    check_fit_allocs && @testset "fit_j2osc_mean_elements! (ForwardDiffJacobian)" begin
        @test length(
            check_allocs(
                (j2oscd, vjd, vr_i, vv_i) -> begin
                    fit_j2osc_mean_elements!(
                        j2oscd,
                        vjd,
                        vr_i,
                        vv_i;
                        jacobian_method = ForwardDiffJacobian(),
                        verbose = false,
                    )
                end,
                (
                    J2OsculatingPropagator{Float64, Float64},
                    Vector{Float64},
                    Vector{SVector{3, Float64}},
                    Vector{SVector{3, Float64}},
                ),
            ),
        ) <= 15
    end

    # == J4 ============================================================================

    @testset "j4_init!" begin
        @test length(
            check_allocs(
                (j4d, orb₀) -> j4_init!(j4d, orb₀),
                (
                    J4Propagator{Float64, Float64},
                    KeplerianElements{TrueAnomaly, Float64, Float64},
                ),
            ),
        ) == 0
    end

    @testset "j4!" begin
        @test length(
            check_allocs((j4d, t) -> j4!(j4d, t), (J4Propagator{Float64, Float64}, Float64))
        ) == 0
    end

    @testset "_mean_elements_jacobian (FiniteDiffJacobian)" begin
        @test length(
            check_allocs(
                (j4d, Δt, x₁, y₁) -> begin
                    SatelliteToolboxPropagators._mean_elements_jacobian(
                        FiniteDiffJacobian(), j4d, Δt, x₁, y₁
                    )
                end,
                (
                    J4Propagator{Float64, Float64},
                    Float64,
                    SVector{6, Float64},
                    SVector{6, Float64},
                ),
            ),
        ) == 0
    end

    @testset "_mean_elements_jacobian (ForwardDiffJacobian)" begin
        @test length(
            check_allocs(
                (j4d, j4d_ad, Δt, x₁, y₁) -> begin
                    SatelliteToolboxPropagators._mean_elements_jacobian(
                        ForwardDiffJacobian(), j4d, Δt, x₁, y₁; pd_ad = j4d_ad
                    )
                end,
                (
                    J4Propagator{Float64, Float64},
                    J4Propagator{Float64, _D},
                    Float64,
                    SVector{6, Float64},
                    SVector{6, Float64},
                ),
            ),
        ) == 0
    end

    check_fit_allocs && @testset "fit_j4_mean_elements! (FiniteDiffJacobian)" begin
        @test length(
            check_allocs(
                (j4d, vjd, vr_i, vv_i) -> begin
                    fit_j4_mean_elements!(
                        j4d,
                        vjd,
                        vr_i,
                        vv_i;
                        jacobian_method = FiniteDiffJacobian(),
                        verbose = false,
                    )
                end,
                (
                    J4Propagator{Float64, Float64},
                    Vector{Float64},
                    Vector{SVector{3, Float64}},
                    Vector{SVector{3, Float64}},
                ),
            ),
        ) <= 14
    end

    check_fit_allocs && @testset "fit_j4_mean_elements! (ForwardDiffJacobian)" begin
        @test length(
            check_allocs(
                (j4d, vjd, vr_i, vv_i) -> begin
                    fit_j4_mean_elements!(
                        j4d,
                        vjd,
                        vr_i,
                        vv_i;
                        jacobian_method = ForwardDiffJacobian(),
                        verbose = false,
                    )
                end,
                (
                    J4Propagator{Float64, Float64},
                    Vector{Float64},
                    Vector{SVector{3, Float64}},
                    Vector{SVector{3, Float64}},
                ),
            ),
        ) <= 14
    end

    # == J4 Osculating =================================================================

    @testset "j4osc_init!" begin
        @test length(
            check_allocs(
                (j4oscd, orb₀) -> j4osc_init!(j4oscd, orb₀),
                (
                    J4OsculatingPropagator{Float64, Float64},
                    KeplerianElements{TrueAnomaly, Float64, Float64},
                ),
            ),
        ) == 0
    end

    @testset "j4osc!" begin
        @test length(
            check_allocs(
                (j4oscd, t) -> j4osc!(j4oscd, t),
                (J4OsculatingPropagator{Float64, Float64}, Float64),
            ),
        ) == 0
    end

    @testset "_mean_elements_jacobian (FiniteDiffJacobian)" begin
        @test length(
            check_allocs(
                (j4oscd, Δt, x₁, y₁) -> begin
                    SatelliteToolboxPropagators._mean_elements_jacobian(
                        FiniteDiffJacobian(), j4oscd, Δt, x₁, y₁
                    )
                end,
                (
                    J4OsculatingPropagator{Float64, Float64},
                    Float64,
                    SVector{6, Float64},
                    SVector{6, Float64},
                ),
            ),
        ) == 0
    end

    @testset "_mean_elements_jacobian (ForwardDiffJacobian)" begin
        @test length(
            check_allocs(
                (j4oscd, j4oscd_ad, Δt, x₁, y₁) -> begin
                    SatelliteToolboxPropagators._mean_elements_jacobian(
                        ForwardDiffJacobian(), j4oscd, Δt, x₁, y₁; pd_ad = j4oscd_ad
                    )
                end,
                (
                    J4OsculatingPropagator{Float64, Float64},
                    J4OsculatingPropagator{Float64, _D},
                    Float64,
                    SVector{6, Float64},
                    SVector{6, Float64},
                ),
            ),
        ) == 0
    end

    check_fit_allocs && @testset "fit_j4osc_mean_elements! (FiniteDiffJacobian)" begin
        @test length(
            check_allocs(
                (j4oscd, vjd, vr_i, vv_i) -> begin
                    fit_j4osc_mean_elements!(
                        j4oscd,
                        vjd,
                        vr_i,
                        vv_i;
                        jacobian_method = FiniteDiffJacobian(),
                        verbose = false,
                    )
                end,
                (
                    J4OsculatingPropagator{Float64, Float64},
                    Vector{Float64},
                    Vector{SVector{3, Float64}},
                    Vector{SVector{3, Float64}},
                ),
            ),
        ) <= 15
    end

    check_fit_allocs && @testset "fit_j4osc_mean_elements! (ForwardDiffJacobian)" begin
        @test length(
            check_allocs(
                (j4oscd, vjd, vr_i, vv_i) -> begin
                    fit_j4osc_mean_elements!(
                        j4oscd,
                        vjd,
                        vr_i,
                        vv_i;
                        jacobian_method = ForwardDiffJacobian(),
                        verbose = false,
                    )
                end,
                (
                    J4OsculatingPropagator{Float64, Float64},
                    Vector{Float64},
                    Vector{SVector{3, Float64}},
                    Vector{SVector{3, Float64}},
                ),
            ),
        ) <= 15
    end
end
