## Description #############################################################################
#
# Tests related to performance and memory allocations.
#
############################################################################################

@testset "Aqua.jl" begin
    Aqua.test_all(
        SatelliteToolboxPropagators;
        ambiguities  = (recursive = false),
        deps_compat  = (check_extras = false),
    )
end

if VERSION >= v"1.12"
    @warn "JET.jl test skipped on Julia 1.12+ due to MethodTableView incompatibility"
else
    @testset "JET Testing" begin
        rep = JET.test_package(SatelliteToolboxPropagators; toplevel_logger=nothing, target_modules=(@__MODULE__,))
    end
end

if VERSION >= v"1.12"
    @warn "Allocation Check skipped on Julia 1.12+ as it is falsely flagging rem2pi internals"
else
    @testset "Allocation Check" begin

        _D = ForwardDiff.Dual{ForwardDiff.Tag{Nothing, Float64}, Float64, 6}

        # == Two-Body ======================================================================

        @testset "twobody_init!" begin
            @test length(
                check_allocs(
                    (tbd, orb₀) -> twobody_init!(tbd, orb₀),
                    (TwoBodyPropagator{Float64, Float64}, KeplerianElements{Float64, Float64})
                )
            ) == 0
        end

        @testset "twobody!" begin
            @test length(
                check_allocs(
                    (tbd, t) -> twobody!(tbd, t),
                    (TwoBodyPropagator{Float64, Float64}, Float64)
                )
            ) == 0
        end

        # == J2 ============================================================================

        @testset "j2_init!" begin
            @test length(
                check_allocs(
                    (j2d, orb₀) -> j2_init!(j2d, orb₀),
                    (J2Propagator{Float64, Float64}, KeplerianElements{Float64, Float64})
                )
            ) == 0
        end

        @testset "j2!" begin
            @test length(
                check_allocs(
                    (j2d, t) -> j2!(j2d, t),
                    (J2Propagator{Float64, Float64}, Float64)
                )
            ) == 0
        end

        @testset "_j2_jacobian (FiniteDiffJacobian)" begin
            @test length(
                check_allocs(
                    (j2d, Δt, x₁, y₁) -> begin
                        SatelliteToolboxPropagators._j2_jacobian(
                            FiniteDiffJacobian(), j2d, Δt, x₁, y₁
                        )
                    end,
                    (J2Propagator{Float64, Float64}, Float64, SVector{6, Float64}, SVector{6, Float64})
                )
            ) == 0
        end

        @testset "_j2_jacobian (ForwardDiffJacobian)" begin
            @test length(
                check_allocs(
                    (j2d, j2d_ad, Δt, x₁, y₁) -> begin
                        SatelliteToolboxPropagators._j2_jacobian(
                            ForwardDiffJacobian(), j2d, Δt, x₁, y₁;
                            j2d_ad = j2d_ad
                        )
                    end,
                    (J2Propagator{Float64, Float64}, J2Propagator{Float64, _D}, Float64, SVector{6, Float64}, SVector{6, Float64})
                )
            ) == 0
        end

        @testset "fit_j2_mean_elements! (FiniteDiffJacobian)" begin
            @test length(
                check_allocs(
                    (j2d, vjd, vr_i, vv_i) -> begin
                        fit_j2_mean_elements!(
                            j2d, vjd, vr_i, vv_i;
                            jacobian_method = FiniteDiffJacobian(),
                            verbose = false,
                        )
                    end,
                    (
                        J2Propagator{Float64, Float64},
                        Vector{Float64},
                        Vector{SVector{3, Float64}},
                        Vector{SVector{3, Float64}},
                    )
                )
            ) <= 13
        end

        @testset "fit_j2_mean_elements! (ForwardDiffJacobian)" begin
            @test length(
                check_allocs(
                    (j2d, vjd, vr_i, vv_i) -> begin
                        fit_j2_mean_elements!(
                            j2d, vjd, vr_i, vv_i;
                            jacobian_method = ForwardDiffJacobian(),
                            verbose = false,
                        )
                    end,
                    (
                        J2Propagator{Float64, Float64},
                        Vector{Float64},
                        Vector{SVector{3, Float64}},
                        Vector{SVector{3, Float64}},
                    )
                )
            ) <= 14
        end

        # == J2 Osculating =================================================================

        @testset "j2osc_init!" begin
            @test length(
                check_allocs(
                    (j2oscd, orb₀) -> j2osc_init!(j2oscd, orb₀),
                    (J2OsculatingPropagator{Float64, Float64}, KeplerianElements{Float64, Float64})
                )
            ) == 0
        end

        @testset "j2osc!" begin
            @test length(
                check_allocs(
                    (j2oscd, t) -> j2osc!(j2oscd, t),
                    (J2OsculatingPropagator{Float64, Float64}, Float64)
                )
            ) == 0
        end

        @testset "_j2osc_jacobian (FiniteDiffJacobian)" begin
            @test length(
                check_allocs(
                    (j2oscd, Δt, x₁, y₁) -> begin
                        SatelliteToolboxPropagators._j2osc_jacobian(
                            FiniteDiffJacobian(), j2oscd, Δt, x₁, y₁
                        )
                    end,
                    (J2OsculatingPropagator{Float64, Float64}, Float64, SVector{6, Float64}, SVector{6, Float64})
                )
            ) == 0
        end

        @testset "_j2osc_jacobian (ForwardDiffJacobian)" begin
            @test length(
                check_allocs(
                    (j2oscd, j2oscd_ad, Δt, x₁, y₁) -> begin
                        SatelliteToolboxPropagators._j2osc_jacobian(
                            ForwardDiffJacobian(), j2oscd, Δt, x₁, y₁;
                            j2oscd_ad = j2oscd_ad
                        )
                    end,
                    (J2OsculatingPropagator{Float64, Float64}, J2OsculatingPropagator{Float64, _D}, Float64, SVector{6, Float64}, SVector{6, Float64})
                )
            ) == 0
        end

        @testset "fit_j2osc_mean_elements! (FiniteDiffJacobian)" begin
            @test length(
                check_allocs(
                    (j2oscd, vjd, vr_i, vv_i) -> begin
                        fit_j2osc_mean_elements!(
                            j2oscd, vjd, vr_i, vv_i;
                            jacobian_method = FiniteDiffJacobian(),
                            verbose = false,
                        )
                    end,
                    (
                        J2OsculatingPropagator{Float64, Float64},
                        Vector{Float64},
                        Vector{SVector{3, Float64}},
                        Vector{SVector{3, Float64}},
                    )
                )
            ) <= 13
        end

        @testset "fit_j2osc_mean_elements! (ForwardDiffJacobian)" begin
            @test length(
                check_allocs(
                    (j2oscd, vjd, vr_i, vv_i) -> begin
                        fit_j2osc_mean_elements!(
                            j2oscd, vjd, vr_i, vv_i;
                            jacobian_method = ForwardDiffJacobian(),
                            verbose = false,
                        )
                    end,
                    (
                        J2OsculatingPropagator{Float64, Float64},
                        Vector{Float64},
                        Vector{SVector{3, Float64}},
                        Vector{SVector{3, Float64}},
                    )
                )
            ) <= 15
        end

        # == J4 ============================================================================

        @testset "j4_init!" begin
            @test length(
                check_allocs(
                    (j4d, orb₀) -> j4_init!(j4d, orb₀),
                    (J4Propagator{Float64, Float64}, KeplerianElements{Float64, Float64})
                )
            ) == 0
        end

        @testset "j4!" begin
            @test length(
                check_allocs(
                    (j4d, t) -> j4!(j4d, t),
                    (J4Propagator{Float64, Float64}, Float64)
                )
            ) == 0
        end

        @testset "_j4_jacobian (FiniteDiffJacobian)" begin
            @test length(
                check_allocs(
                    (j4d, Δt, x₁, y₁) -> begin
                        SatelliteToolboxPropagators._j4_jacobian(
                            FiniteDiffJacobian(), j4d, Δt, x₁, y₁
                        )
                    end,
                    (J4Propagator{Float64, Float64}, Float64, SVector{6, Float64}, SVector{6, Float64})
                )
            ) == 0
        end

        @testset "_j4_jacobian (ForwardDiffJacobian)" begin
            @test length(
                check_allocs(
                    (j4d, j4d_ad, Δt, x₁, y₁) -> begin
                        SatelliteToolboxPropagators._j4_jacobian(
                            ForwardDiffJacobian(), j4d, Δt, x₁, y₁;
                            j4d_ad = j4d_ad
                        )
                    end,
                    (J4Propagator{Float64, Float64}, J4Propagator{Float64, _D}, Float64, SVector{6, Float64}, SVector{6, Float64})
                )
            ) == 0
        end

        @testset "fit_j4_mean_elements! (FiniteDiffJacobian)" begin
            @test length(
                check_allocs(
                    (j4d, vjd, vr_i, vv_i) -> begin
                        fit_j4_mean_elements!(
                            j4d, vjd, vr_i, vv_i;
                            jacobian_method = FiniteDiffJacobian(),
                            verbose = false,
                        )
                    end,
                    (
                        J4Propagator{Float64, Float64},
                        Vector{Float64},
                        Vector{SVector{3, Float64}},
                        Vector{SVector{3, Float64}},
                    )
                )
            ) <= 13
        end

        @testset "fit_j4_mean_elements! (ForwardDiffJacobian)" begin
            @test length(
                check_allocs(
                    (j4d, vjd, vr_i, vv_i) -> begin
                        fit_j4_mean_elements!(
                            j4d, vjd, vr_i, vv_i;
                            jacobian_method = ForwardDiffJacobian(),
                            verbose = false,
                        )
                    end,
                    (
                        J4Propagator{Float64, Float64},
                        Vector{Float64},
                        Vector{SVector{3, Float64}},
                        Vector{SVector{3, Float64}},
                    )
                )
            ) <= 14
        end

        # == J4 Osculating =================================================================

        @testset "j4osc_init!" begin
            @test length(
                check_allocs(
                    (j4oscd, orb₀) -> j4osc_init!(j4oscd, orb₀),
                    (J4OsculatingPropagator{Float64, Float64}, KeplerianElements{Float64, Float64})
                )
            ) == 0
        end

        @testset "j4osc!" begin
            @test length(
                check_allocs(
                    (j4oscd, t) -> j4osc!(j4oscd, t),
                    (J4OsculatingPropagator{Float64, Float64}, Float64)
                )
            ) == 0
        end

        @testset "_j4osc_jacobian (FiniteDiffJacobian)" begin
            @test length(
                check_allocs(
                    (j4oscd, Δt, x₁, y₁) -> begin
                        SatelliteToolboxPropagators._j4osc_jacobian(
                            FiniteDiffJacobian(), j4oscd, Δt, x₁, y₁
                        )
                    end,
                    (J4OsculatingPropagator{Float64, Float64}, Float64, SVector{6, Float64}, SVector{6, Float64})
                )
            ) == 0
        end

        @testset "_j4osc_jacobian (ForwardDiffJacobian)" begin
            @test length(
                check_allocs(
                    (j4oscd, j4oscd_ad, Δt, x₁, y₁) -> begin
                        SatelliteToolboxPropagators._j4osc_jacobian(
                            ForwardDiffJacobian(), j4oscd, Δt, x₁, y₁;
                            j4oscd_ad = j4oscd_ad
                        )
                    end,
                    (J4OsculatingPropagator{Float64, Float64}, J4OsculatingPropagator{Float64, _D}, Float64, SVector{6, Float64}, SVector{6, Float64})
                )
            ) == 0
        end

        @testset "fit_j4osc_mean_elements! (FiniteDiffJacobian)" begin
            @test length(
                check_allocs(
                    (j4oscd, vjd, vr_i, vv_i) -> begin
                        fit_j4osc_mean_elements!(
                            j4oscd, vjd, vr_i, vv_i;
                            jacobian_method = FiniteDiffJacobian(),
                            verbose = false,
                        )
                    end,
                    (
                        J4OsculatingPropagator{Float64, Float64},
                        Vector{Float64},
                        Vector{SVector{3, Float64}},
                        Vector{SVector{3, Float64}},
                    )
                )
            ) <= 13
        end

        @testset "fit_j4osc_mean_elements! (ForwardDiffJacobian)" begin
            @test length(
                check_allocs(
                    (j4oscd, vjd, vr_i, vv_i) -> begin
                        fit_j4osc_mean_elements!(
                            j4oscd, vjd, vr_i, vv_i;
                            jacobian_method = ForwardDiffJacobian(),
                            verbose = false,
                        )
                    end,
                    (
                        J4OsculatingPropagator{Float64, Float64},
                        Vector{Float64},
                        Vector{SVector{3, Float64}},
                        Vector{SVector{3, Float64}},
                    )
                )
            ) <= 15
        end

    end
end
