## Description #############################################################################
#
#   Tests related to the two-body orbit propagator.
#
## References ##############################################################################
#
# [1] Vallado, D. A (2013). Fundamentals of Astrodynamics and Applications. Microcosm Press,
#     Hawthorn, CA, USA.
#
############################################################################################

############################################################################################
#                                       Test Results                                       #
############################################################################################
#
# == Scenario 01 ===========================================================================
#
# Example 2-4. Solving Kepler's Problem [1, p. 94-95].
#
#   Initial position:
#
#       r0_ijk = + 1131.340 i - 2282.343 j + 6672.423 k [km]
#       v0_ijk = - 5.64305  i + 4.30333  j + 2.42879  k [km/s]
#
#    After 40 min., the solution of Kepler's problem leads to the following
#    position:
#
#       rf_ijk = - 4219.7527 i + 4363.0292 j - 3958.7666 k [km]
#       vf_ijk = - 3.689866  i - 1.916735  j - 6.112511  k [km/s]
#
############################################################################################

@testset "Two-Body Orbit Propagator" verbose = true begin
    jd₀ = date_to_jd(1986, 6, 19, 18, 35, 0)
    jd₁ = date_to_jd(1986, 6, 19, 19, 15, 0)

    # == Constructor =======================================================================

    @testset "Constructor" begin
        orb = KeplerianElements{MeanAnomaly}(0.0, 8000.0e3, 0.0, 0.0, 0.0, 0.0, 0.0)
        tbd = TwoBodyPropagator{Float64, Float64}(orb, orb, 0, 0, 0)

        # Test some random fields.
        @test tbd.Δt == 0
        @test tbd.orb₀ == orb
        @test tbd.orbk == orb
        @test tbd.μ == 0
        @test tbd.Δt == 0
        @test tbd.n₀ == 0

        orb = KeplerianElements(0.0, 8_000_000, 0, 0, 0, 0, 0)
        orbp = Propagators.init(Val(:TwoBody), orb)
        @test orbp.tbd.orb₀ isa KeplerianElements{MeanAnomaly, Float64, Float64}
    end

    # == General API Functions =============================================================

    @testset "General API Functions" begin
        orb = KeplerianElements{MeanAnomaly}(0.0, 8000.0e3, 0.0, 0.0, 0.0, 0.0, 0.0)
        orbp = Propagators.init(Val(:TwoBody), orb)
        @test Propagators.name(orbp) == "Two-Body Orbit Propagator"
    end

    # == Float64 ===========================================================================

    @testset "Float64" begin
        T = Float64

        orb = rv_to_kepler(
            [1131340.0, -2282343.0, 6672423.0],
            [-5643.05, 4303.33, 2428.79],
            date_to_jd(1986, 6, 19, 18, 35, 0),
        )

        orbp = Propagators.init(Val(:TwoBody), orb)

        orbk = Propagators.mean_elements(orbp)
        @test orbk isa KeplerianElements{MeanAnomaly, Float64, Float64}
        @test orbk.t ≈ orb.t
        @test orbk.a ≈ orb.a
        @test orbk.e ≈ orb.e
        @test orbk.i ≈ orb.i
        @test orbk.Ω ≈ orb.Ω
        @test orbk.ω ≈ orb.ω
        @test orbk.f ≈ orb.f

        r, v = Propagators.step!(orbp, 40 * 60)

        @test r[1] / 1000 ≈ -4219.7527 atol = 1e-3
        @test r[2] / 1000 ≈ +4363.0292 atol = 1e-3
        @test r[3] / 1000 ≈ -3958.7666 atol = 1e-3
        @test eltype(r) == T

        @test v[1] / 1000 ≈ +3.689866 atol = 1e-6
        @test v[2] / 1000 ≈ -1.916735 atol = 1e-6
        @test v[3] / 1000 ≈ -6.112511 atol = 1e-6
        @test eltype(v) == T

        r, v = Propagators.propagate_to_epoch!(orbp, jd₁)

        @test r[1] / 1000 ≈ -4219.7527 atol = 1e-3
        @test r[2] / 1000 ≈ +4363.0292 atol = 1e-3
        @test r[3] / 1000 ≈ -3958.7666 atol = 1e-3
        @test eltype(r) == T

        @test v[1] / 1000 ≈ +3.689866 atol = 1e-6
        @test v[2] / 1000 ≈ -1.916735 atol = 1e-6
        @test v[3] / 1000 ≈ -6.112511 atol = 1e-6
        @test eltype(v) == T

        # Test in-place initialization.
        orbp = OrbitPropagatorTwoBody(TwoBodyPropagator{Float64, T}())
        orbp.tbd.μ = TBC_M0
        Propagators.init!(orbp, orb)

        r, v = Propagators.step!(orbp, 40 * 60)

        @test r[1] / 1000 ≈ -4219.7527 atol = 1e-3
        @test r[2] / 1000 ≈ +4363.0292 atol = 1e-3
        @test r[3] / 1000 ≈ -3958.7666 atol = 1e-3
        @test eltype(r) == T

        @test v[1] / 1000 ≈ +3.689866 atol = 1e-6
        @test v[2] / 1000 ≈ -1.916735 atol = 1e-6
        @test v[3] / 1000 ≈ -6.112511 atol = 1e-6
        @test eltype(v) == T

        # Test simultaneous initialization and propagation.
        r, v, orbp = Propagators.propagate(Val(:TwoBody), 40 * 60, orb)

        orbk = Propagators.mean_elements(orbp)
        @test orbk isa KeplerianElements{MeanAnomaly, Float64, Float64}

        @test r[1] / 1000 ≈ -4219.7527 atol = 1e-3
        @test r[2] / 1000 ≈ +4363.0292 atol = 1e-3
        @test r[3] / 1000 ≈ -3958.7666 atol = 1e-3
        @test eltype(r) == T

        @test v[1] / 1000 ≈ +3.689866 atol = 1e-6
        @test v[2] / 1000 ≈ -1.916735 atol = 1e-6
        @test v[3] / 1000 ≈ -6.112511 atol = 1e-6
        @test eltype(v) == T

        r, v, orbp = Propagators.propagate_to_epoch(Val(:TwoBody), jd₁, orb)

        orbk = Propagators.mean_elements(orbp)
        @test orbk isa KeplerianElements{MeanAnomaly, Float64, Float64}

        @test r[1] / 1000 ≈ -4219.7527 atol = 1e-3
        @test r[2] / 1000 ≈ +4363.0292 atol = 1e-3
        @test r[3] / 1000 ≈ -3958.7666 atol = 1e-3
        @test eltype(r) == T

        @test v[1] / 1000 ≈ +3.689866 atol = 1e-6
        @test v[2] / 1000 ≈ -1.916735 atol = 1e-6
        @test v[3] / 1000 ≈ -6.112511 atol = 1e-6
        @test eltype(v) == T

        r, v, tbd = twobody(40 * 60, orb)

        @test tbd isa TwoBodyPropagator{Float64, Float64}

        @test r[1] / 1000 ≈ -4219.7527 atol = 1e-3
        @test r[2] / 1000 ≈ +4363.0292 atol = 1e-3
        @test r[3] / 1000 ≈ -3958.7666 atol = 1e-3
        @test eltype(r) == T

        @test v[1] / 1000 ≈ +3.689866 atol = 1e-6
        @test v[2] / 1000 ≈ -1.916735 atol = 1e-6
        @test v[3] / 1000 ≈ -6.112511 atol = 1e-6
        @test eltype(v) == T
    end

    # == Float32 ===========================================================================

    @testset "Float32" begin
        T = Float32

        # We saw rounding problems when converting from state vector to Keplerian elements
        # in GitHub actions. The source was the processor that does not have FMA (fused
        # multiply add). Hence, we will obtain the Keplerian elements using `Float64` but
        # propagate the orbit using `Float32`.
        orb = rv_to_kepler(
            [1131340.0, -2282343.0, 6672423.0],
            [-5643.05, 4303.33, 2428.79],
            date_to_jd(1986, 6, 19, 18, 35, 0),
        )

        orbp = Propagators.init(Val(:TwoBody), orb; m0 = TBC_M0_F32)

        orbk = Propagators.mean_elements(orbp)
        @test orbk isa KeplerianElements{MeanAnomaly, Float64, Float32}
        @test orbk.t ≈ orb.t
        @test orbk.a ≈ orb.a
        @test orbk.e ≈ orb.e
        @test orbk.i ≈ orb.i
        @test orbk.Ω ≈ orb.Ω
        @test orbk.ω ≈ orb.ω
        @test orbk.f ≈ orb.f

        r, v = Propagators.step!(orbp, 40 * 60)

        @test r[1] / 1000 ≈ -4219.7527 atol = 5e-1
        @test r[2] / 1000 ≈ +4363.0292 atol = 5e-1
        @test r[3] / 1000 ≈ -3958.7666 atol = 5e-1
        @test eltype(r) == T

        @test v[1] / 1000 ≈ +3.689866 atol = 5e-3
        @test v[2] / 1000 ≈ -1.916735 atol = 5e-3
        @test v[3] / 1000 ≈ -6.112511 atol = 5e-3
        @test eltype(v) == T

        r, v = Propagators.propagate_to_epoch!(orbp, jd₁)

        @test r[1] / 1000 ≈ -4219.7527 atol = 5e-1
        @test r[2] / 1000 ≈ +4363.0292 atol = 5e-1
        @test r[3] / 1000 ≈ -3958.7666 atol = 5e-1
        @test eltype(r) == T

        @test v[1] / 1000 ≈ +3.689866 atol = 5e-3
        @test v[2] / 1000 ≈ -1.916735 atol = 5e-3
        @test v[3] / 1000 ≈ -6.112511 atol = 5e-3
        @test eltype(v) == T

        # Test in-place initialization.
        orbp = OrbitPropagatorTwoBody(TwoBodyPropagator{Float64, T}())
        orbp.tbd.μ = TBC_M0_F32
        Propagators.init!(orbp, orb)

        r, v = Propagators.step!(orbp, 40 * 60)

        @test r[1] / 1000 ≈ -4219.7527 atol = 5e-1
        @test r[2] / 1000 ≈ +4363.0292 atol = 5e-1
        @test r[3] / 1000 ≈ -3958.7666 atol = 5e-1
        @test eltype(r) == T

        @test v[1] / 1000 ≈ +3.689866 atol = 5e-3
        @test v[2] / 1000 ≈ -1.916735 atol = 5e-3
        @test v[3] / 1000 ≈ -6.112511 atol = 5e-3
        @test eltype(v) == T

        # Test simultaneous initialization and propagation.
        r, v, orbp = Propagators.propagate(Val(:TwoBody), 40 * 60, orb; m0 = TBC_M0_F32)

        orbk = Propagators.mean_elements(orbp)
        @test orbk isa KeplerianElements{MeanAnomaly, Float64, Float32}

        @test r[1] / 1000 ≈ -4219.7527 atol = 5e-1
        @test r[2] / 1000 ≈ +4363.0292 atol = 5e-1
        @test r[3] / 1000 ≈ -3958.7666 atol = 5e-1
        @test eltype(r) == T

        @test v[1] / 1000 ≈ +3.689866 atol = 5e-3
        @test v[2] / 1000 ≈ -1.916735 atol = 5e-3
        @test v[3] / 1000 ≈ -6.112511 atol = 5e-3
        @test eltype(v) == T

        r, v, orbp = Propagators.propagate_to_epoch(
            Val(:TwoBody), jd₁, orb; m0 = TBC_M0_F32
        )

        orbk = Propagators.mean_elements(orbp)
        @test orbk isa KeplerianElements{MeanAnomaly, Float64, Float32}

        @test r[1] / 1000 ≈ -4219.7527 atol = 5e-1
        @test r[2] / 1000 ≈ +4363.0292 atol = 5e-1
        @test r[3] / 1000 ≈ -3958.7666 atol = 5e-1
        @test eltype(r) == T

        @test v[1] / 1000 ≈ +3.689866 atol = 5e-3
        @test v[2] / 1000 ≈ -1.916735 atol = 5e-3
        @test v[3] / 1000 ≈ -6.112511 atol = 5e-3
        @test eltype(v) == T

        r, v, tbd = twobody(40 * 60, orb; m0 = TBC_M0_F32)

        @test tbd isa TwoBodyPropagator{Float64, Float32}

        @test r[1] / 1000 ≈ -4219.7527 atol = 5e-1
        @test r[2] / 1000 ≈ +4363.0292 atol = 5e-1
        @test r[3] / 1000 ≈ -3958.7666 atol = 5e-1
        @test eltype(r) == T

        @test v[1] / 1000 ≈ +3.689866 atol = 5e-3
        @test v[2] / 1000 ≈ -1.916735 atol = 5e-3
        @test v[3] / 1000 ≈ -6.112511 atol = 5e-3
        @test eltype(v) == T
    end
end

@testset "Copying Structure" verbose = true begin
    for (T, tbc) in ((Float64, TBC_M0), (Float32, TBC_M0_F32))
        @testset "$T" begin
            jd₀ = date_to_jd(2023, 1, 1, 0, 0, 0)

            orb = KeplerianElements(
                jd₀,
                T(8000e3),
                T(0.015),
                T(28.5) |> deg2rad,
                T(100) |> deg2rad,
                T(400) |> deg2rad,
                T(45) |> deg2rad,
            )

            orbp = Propagators.init(Val(:TwoBody), orb; m0 = tbc)
            new_orbp = copy(orbp)

            for f in fieldnames(typeof(orbp.tbd))
                @test getfield(new_orbp.tbd, f) == getfield(orbp.tbd, f)
            end

            new_orbp.tbd.Δt = 1000
            @test new_orbp.tbd.Δt != orbp.tbd.Δt
        end
    end
end

@testset "Fitting Mean Elements for the Two-Body Propagator" verbose = true begin
    # We just need to run the two-body propagator, obtain the state vectors, convert them to
    # mean Keplerian elements, and compare with the original mean elements.

    orb_input = KeplerianElements(
        DateTime("2023-01-01") |> datetime2julian,
        7130.982e3,
        0.001111,
        98.405 |> deg2rad,
        90 |> deg2rad,
        200 |> deg2rad,
        45 |> deg2rad,
    )

    # Generate the osculating elements.
    orbp = Propagators.init(Val(:TwoBody), orb_input)
    ret  = Propagators.propagate!.(orbp, 0:10:12_000)
    vr_i = first.(ret)
    vv_i = last.(ret)
    vjd  = Propagators.epoch(orbp) .+ (0:10:12_000) ./ 86400

    @testset "Without Initial Guess" begin
        for jacobian_method in (FiniteDiffJacobian(), ForwardDiffJacobian())
            orb, P, stats = redirect_stdout(devnull) do
                Propagators.fit_mean_elements(
                    Val(:TwoBody),
                    vjd,
                    vr_i,
                    vv_i;
                    jacobian_method     = jacobian_method,
                    mean_elements_epoch = vjd[begin],
                )
            end

            # The statistics must describe a converged fit with a negligible residue,
            # since the measurements were generated by the same propagator.
            @test P isa SMatrix{6, 6, Float64}
            @test stats.converged === true
            @test 1 <= stats.iterations <= 50
            @test stats.position_rmse isa Float64
            @test stats.velocity_rmse isa Float64
            @test stats.total_rmse isa Float64
            @test 0 <= stats.position_rmse < 1
            @test 0 <= stats.velocity_rmse < 1e-3
            @test 0 <= stats.total_rmse < 1

            @test orb isa KeplerianElements{MeanAnomaly, Float64, Float64}
            @test P isa SMatrix{6, 6, Float64}
            @test orb.t ≈ orb_input.t
            @test orb.a ≈ orb_input.a
            @test orb.e ≈ orb_input.e
            @test orb.i ≈ orb_input.i
            @test orb.Ω ≈ orb_input.Ω
            @test orb.ω ≈ orb_input.ω
            @test orb.f ≈ orb_input.f atol = 1e-6

            orb, P = redirect_stdout(devnull) do
                Propagators.fit_mean_elements!(
                    orbp,
                    vjd,
                    vr_i,
                    vv_i;
                    jacobian_method     = jacobian_method,
                    mean_elements_epoch = vjd[begin],
                )
            end

            @test orb isa KeplerianElements{MeanAnomaly, Float64, Float64}
            @test orb.t ≈ orb_input.t
            @test orb.a ≈ orb_input.a
            @test orb.e ≈ orb_input.e
            @test orb.i ≈ orb_input.i
            @test orb.Ω ≈ orb_input.Ω
            @test orb.ω ≈ orb_input.ω
            @test orb.f ≈ orb_input.f atol = 1e-6

            # The propagator must be initialized with the fitted elements.
            @test Propagators.mean_elements(orbp).a ≈ orb.a
        end
    end

    @testset "With Initial Guess" begin
        orb, ~ = redirect_stdout(devnull) do
            Propagators.fit_mean_elements(
                Val(:TwoBody),
                vjd,
                vr_i,
                vv_i;
                max_iterations = 3,
                mean_elements_epoch = vjd[begin],
                initial_guess = orb_input,
            )
        end

        @test orb.t ≈ orb_input.t
        @test orb.a ≈ orb_input.a
        @test orb.e ≈ orb_input.e
        @test orb.i ≈ orb_input.i
        @test orb.Ω ≈ orb_input.Ω
        @test orb.ω ≈ orb_input.ω
        @test orb.f ≈ orb_input.f atol = 1e-6
    end

    @testset "Without Initial Guess and Updating the Epoch" begin
        orb, ~ = redirect_stdout(devnull) do
            Propagators.fit_mean_elements(
                Val(:TwoBody), vjd, vr_i, vv_i; mean_elements_epoch = vjd[begin] + 1
            )
        end

        # The two-body propagator only changes the mean anomaly.
        @test orb.t ≈ orb_input.t + 1
        @test orb.a ≈ orb_input.a
        @test orb.e ≈ orb_input.e
        @test orb.i ≈ orb_input.i
        @test orb.Ω ≈ orb_input.Ω
        @test orb.ω ≈ orb_input.ω

        n₀ = √(TBC_M0 / orb_input.a^3)
        M₁ = mod(mean_anomaly(orb_input) + n₀ * 86400, 2π)
        @test orb.anomaly ≈ M₁ atol = 1e-6
    end

    @testset "Constants Keyword" begin
        # The number type of the constants selects the type of the fit.
        orb, P, stats = redirect_stdout(devnull) do
            Propagators.fit_mean_elements(Val(:TwoBody), vjd, vr_i, vv_i; m0 = TBC_M0_F32)
        end

        @test orb isa KeplerianElements{MeanAnomaly, Float64, Float32}
        @test P isa SMatrix{6, 6, Float32}
        @test stats.position_rmse isa Float32
        @test orb.a ≈ orb_input.a rtol = 1e-3

        # The allocating function must use the selected constants, leading to the same
        # result obtained with a propagator initialized with them.
        orb_alt, P_alt = fit_twobody_mean_elements(vjd, vr_i, vv_i; m0 = 3.986e14, verbose = false)
        pd = twobody_init(orb_input; m0 = 3.986e14)
        orb_ref, P_ref = fit_twobody_mean_elements!(pd, vjd, vr_i, vv_i; verbose = false)

        @test orb_alt == orb_ref
        @test P_alt == P_ref
    end

    @testset "Errors" begin
        @test_throws ArgumentError Propagators.fit_mean_elements(
            Val(:TwoBody), vjd[1:(end - 1)], vr_i, vv_i
        )
        @test_throws ArgumentError Propagators.fit_mean_elements(
            Val(:TwoBody), vjd, vr_i[1:(end - 1)], vv_i
        )
        @test_throws ArgumentError Propagators.fit_mean_elements(
            Val(:TwoBody), vjd, vr_i, vv_i[1:(end - 1)]
        )
        @test_throws ArgumentError Propagators.fit_mean_elements(
            Val(:TwoBody), vjd, vr_i, vv_i; weight_vector = [1, 2, 3, 4, 5]
        )

        # == Invalid maximum number of iterations ==========================================

        @test_throws ArgumentError Propagators.fit_mean_elements(
            Val(:TwoBody), vjd, vr_i, vv_i; max_iterations = 0
        )
    end
end

@testset "Update Two-Body Mean Elements Epoch" begin
    orb_input = KeplerianElements(
        DateTime("2023-01-01") |> datetime2julian,
        7130.982e3,
        0.001111,
        98.405 |> deg2rad,
        90 |> deg2rad,
        200 |> deg2rad,
        45 |> deg2rad,
    )

    # Only the mean anomaly changes in the two-body propagator, advancing by the mean
    # motion times the elapsed time.
    n₀ = √(TBC_M0 / orb_input.a^3)
    M₁ = mod(mean_anomaly(orb_input) + n₀ * 86400, 2π)

    for new_epoch in (DateTime("2023-01-02"), orb_input.t + 1)
        orb = update_twobody_mean_elements_epoch(orb_input, new_epoch)

        @test orb isa KeplerianElements{MeanAnomaly, Float64, Float64}
        @test orb.t ≈ orb_input.t + 1
        @test orb.a == orb_input.a
        @test orb.e == orb_input.e
        @test orb.i == orb_input.i
        @test orb.Ω == orb_input.Ω
        @test orb.ω == orb_input.ω
        @test orb.anomaly ≈ M₁ atol = 1e-9

        tbd = twobody_init(orb_input)
        orb = update_twobody_mean_elements_epoch!(tbd, orb_input, new_epoch)

        @test orb.anomaly ≈ M₁ atol = 1e-9
        @test tbd.orb₀ == orb
    end

    # The function must also work when the element type is not `Float64`.
    orb_input_f32 = convert(KeplerianElements{TrueAnomaly, Float64, Float32}, orb_input)
    orb_f32 = update_twobody_mean_elements_epoch(orb_input_f32, DateTime("2023-01-02"))

    @test orb_f32 isa KeplerianElements{MeanAnomaly, Float64, Float32}
    @test orb_f32.anomaly ≈ M₁ atol = 1e-3
    # The keyword selects the constants, leading to the same result obtained with a
    # propagator initialized with them.
    new_epoch = DateTime("2023-01-02")
    orb_alt   = update_twobody_mean_elements_epoch(orb_input, new_epoch; m0 = 3.986e14)
    pd        = twobody_init(orb_input; m0 = 3.986e14)

    @test orb_alt == update_twobody_mean_elements_epoch!(pd, orb_input, new_epoch)
    @test orb_alt != update_twobody_mean_elements_epoch(orb_input, new_epoch)

end
