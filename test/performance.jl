## Description #############################################################################
#
# Tests related to memory allocations.
#
############################################################################################

@testset "Allocation Check" begin
    _D = ForwardDiff.Dual{ForwardDiff.Tag{Nothing, Float64}, Float64, 6}

    # The kernels and the Jacobians must not allocate on any Julia version.

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
end

@testset "Fitting Allocations" begin
    # The fitting allocates only a fixed set of workspaces. Hence, the number of allocations
    # must not depend on the number of measurements or iterations. We check this at runtime
    # because the static analysis of AllocCheck also follows the code paths used only when
    # `verbose = true`, whose allocations inside Base I/O internals vary with the platform,
    # the loaded packages, and the AllocCheck version.

    orb_input = KeplerianElements(
        DateTime("2023-01-01") |> datetime2julian,
        7130.982e3,
        0.001111,
        98.405 |> deg2rad,
        90 |> deg2rad,
        200 |> deg2rad,
        45 |> deg2rad,
    )

    function fit_allocations(
        propagator::Val, jacobian_method, num_measurements, num_iterations
    )
        orbp = Propagators.init(propagator, orb_input)
        vt   = range(0, 6000; length = num_measurements)
        ret  = Propagators.propagate!.(orbp, vt)
        vr_i = first.(ret)
        vv_i = last.(ret)
        vjd  = collect(orb_input.t .+ vt ./ 86400)

        # Setting `atol` and `rtol` to 0 disables the convergence check, so the algorithm
        # performs exactly `num_iterations` iterations.
        fit!() = Propagators.fit_mean_elements!(
            orbp,
            vjd,
            vr_i,
            vv_i;
            atol            = 0,
            rtol            = 0,
            jacobian_method = jacobian_method,
            max_iterations  = num_iterations,
            verbose         = false,
        )

        # Compile the function before measuring the allocations.
        fit!()

        return @allocations fit!()
    end

    for propagator in (Val(:J2), Val(:J2osc), Val(:J4), Val(:J4osc), Val(:TwoBody)),
        jacobian_method in (FiniteDiffJacobian(), ForwardDiffJacobian())

        num_allocations = fit_allocations(propagator, jacobian_method, 50, 2)

        @test fit_allocations(propagator, jacobian_method, 200, 2) == num_allocations
        @test fit_allocations(propagator, jacobian_method, 50, 6) == num_allocations
    end
end
