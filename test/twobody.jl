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
        orb = KeplerianElements(0.0, 8000.0e3, 0.0, 0.0, 0.0, 0.0, 0.0)
        tbd = TwoBodyPropagator{Float64, Float64}(orb, orb, 0, 0, 0, 0)

        # Test some random fields.
        @test tbd.Δt   == 0
        @test tbd.orb₀ == orb
        @test tbd.orbk == orb
        @test tbd.μ    == 0
        @test tbd.Δt   == 0
        @test tbd.M₀   == 0
        @test tbd.n₀   == 0
    end

    # == General API Functions =============================================================

    @testset "General API Functions" begin
        orb = KeplerianElements(0.0, 8000.0e3, 0.0, 0.0, 0.0, 0.0, 0.0)
        orbp = Propagators.init(Val(:TwoBody), orb)
        @test Propagators.name(orbp) == "Two-Body Orbit Propagator"
    end

    # == Float64 ===========================================================================

    @testset "Float64" begin
        T = Float64

        orb = rv_to_kepler(
            [1131340.0, -2282343.0, 6672423.0],
            [-5643.05, 4303.33, 2428.79],
            date_to_jd(1986, 6, 19, 18, 35, 0)
        )

        orbp = Propagators.init(Val(:TwoBody), orb)

        orbk = Propagators.mean_elements(orbp)
        @test orbk isa KeplerianElements{Float64, Float64}
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
        orbp.tbd.μ = tbc_m0
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
        @test orbk isa KeplerianElements{Float64, Float64}

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
        @test orbk isa KeplerianElements{Float64, Float64}

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
            date_to_jd(1986, 6, 19, 18, 35, 0)
        )

        orbp = Propagators.init(Val(:TwoBody), orb; m0 = tbc_m0_f32)

        orbk = Propagators.mean_elements(orbp)
        @test orbk isa KeplerianElements{Float64, Float32}
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
        orbp.tbd.μ = tbc_m0_f32
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
        r, v, orbp = Propagators.propagate(Val(:TwoBody), 40 * 60, orb; m0 = tbc_m0_f32)

        orbk = Propagators.mean_elements(orbp)
        @test orbk isa KeplerianElements{Float64, Float32}

        @test r[1] / 1000 ≈ -4219.7527 atol = 5e-1
        @test r[2] / 1000 ≈ +4363.0292 atol = 5e-1
        @test r[3] / 1000 ≈ -3958.7666 atol = 5e-1
        @test eltype(r) == T

        @test v[1] / 1000 ≈ +3.689866 atol = 5e-3
        @test v[2] / 1000 ≈ -1.916735 atol = 5e-3
        @test v[3] / 1000 ≈ -6.112511 atol = 5e-3
        @test eltype(v) == T

        r, v, orbp = Propagators.propagate_to_epoch(Val(:TwoBody), jd₁, orb; m0 = tbc_m0_f32)

        orbk = Propagators.mean_elements(orbp)
        @test orbk isa KeplerianElements{Float64, Float32}

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
    for (T, tbc) in ((Float64, tbc_m0), (Float32, tbc_m0_f32))
        @testset "$T" begin
            jd₀ = date_to_jd(2023, 1, 1, 0, 0, 0)

            orb = KeplerianElements(
                jd₀,
                T(8000e3),
                T(0.015),
                T(28.5) |> deg2rad,
                T(100)  |> deg2rad,
                T(400)  |> deg2rad,
                T(45)   |> deg2rad
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
