## Description #############################################################################
#
#   Tests related to the J2 orbit propagator.
#
############################################################################################

############################################################################################
#                                       Test Results                                       #
############################################################################################
#
# == Scenario 01 ===========================================================================
#
# Source: https://github.com/JuliaSpace/SatelliteToolbox.jl/issues/91
#
# The STK provided the following result for the J2 orbit propagator.
#
# Initial mean orbital elements:
#
#            Orbit epoch : 2023-01-01T00:00:00.000
#        Semi-major axis : 8000.000 km
#           Eccentricity : 0.015
#            Inclination : 28.5°
#                   RAAN : 100°
#    Argument of perigee : 200°
#   Initial mean anomaly : 45°
#
#
# State vector on 2023-01-05T00:00:00.000
#
#   r = [-6849.654348, -2253.059809, 3574.529667] km
#   v = [    1.656142,    -6.699518,   -1.266334] km/s
#
############################################################################################

@testset "J2 Orbit Propagator" verbose = true begin
    jd₀ = date_to_jd(2023, 1, 1, 0, 0, 0)
    jd₁ = date_to_jd(2023, 1, 5, 0, 0, 0)

    # == Constructor =======================================================================

    @testset "Constructor" begin
        orb = KeplerianElements(0.0, 8000.0e3, 0.0, 0.0, 0.0, 0.0, 0.0)
        j2d = J2Propagator{Float64, Float64}(orb, orb, j2c_egm2008, 0, 0, 0, 0, 0)

        # Test some random fields.
        @test j2d.Δt  == 0
        @test j2d.n̄   == 0
        @test j2d.j2c == j2c_egm2008
    end

    # == General API Functions =============================================================

    @testset "General API Functions" begin
        orb = KeplerianElements(0.0, 8000.0e3, 0.0, 0.0, 0.0, 0.0, 0.0)
        orbp = Propagators.init(Val(:J2), orb)
        @test Propagators.name(orbp) == "J2 Orbit Propagator"
    end

    # == Float64 ===========================================================================

    @testset "Float64" begin
        T = Float64

        orb = KeplerianElements(
            jd₀,
            T(8000e3),
            T(0.015),
            T(28.5) |> deg2rad,
            T(100)  |> deg2rad,
            T(200)  |> deg2rad,
            T(45)   |> deg2rad
        )

        orbp = Propagators.init(Val(:J2), orb; j2c = j2c_egm2008)

        orbk = Propagators.mean_elements(orbp)
        @test orbk isa KeplerianElements{Float64, Float64}
        @test orbk.t ≈ orb.t
        @test orbk.a ≈ orb.a
        @test orbk.e ≈ orb.e
        @test orbk.i ≈ orb.i
        @test orbk.Ω ≈ orb.Ω
        @test orbk.ω ≈ orb.ω
        @test orbk.f ≈ orb.f

        r, v = Propagators.step!(orbp, (jd₁ - jd₀) * 86400)

        @test r[1] / 1000 ≈ -6849.654348 atol = 5e-3
        @test r[2] / 1000 ≈ -2253.059809 atol = 5e-3
        @test r[3] / 1000 ≈ +3574.529667 atol = 5e-3
        @test eltype(r) == T

        @test v[1] / 1000 ≈ +1.656142 atol = 1e-5
        @test v[2] / 1000 ≈ -6.699518 atol = 1e-5
        @test v[3] / 1000 ≈ -1.266334 atol = 1e-5
        @test eltype(v) == T

        r, v = Propagators.propagate_to_epoch!(orbp, jd₁)

        @test r[1] / 1000 ≈ -6849.654348 atol = 5e-3
        @test r[2] / 1000 ≈ -2253.059809 atol = 5e-3
        @test r[3] / 1000 ≈ +3574.529667 atol = 5e-3
        @test eltype(r) == T

        @test v[1] / 1000 ≈ +1.656142 atol = 1e-5
        @test v[2] / 1000 ≈ -6.699518 atol = 1e-5
        @test v[3] / 1000 ≈ -1.266334 atol = 1e-5
        @test eltype(v) == T

        # Test in-place initialization.
        orbp = OrbitPropagatorJ2(J2Propagator{Float64, T}())
        orbp.j2d.j2c = j2c_egm2008
        Propagators.init!(orbp, orb)

        r, v = Propagators.step!(orbp, (jd₁ - jd₀) * 86400)

        @test r[1] / 1000 ≈ -6849.654348 atol = 5e-3
        @test r[2] / 1000 ≈ -2253.059809 atol = 5e-3
        @test r[3] / 1000 ≈ +3574.529667 atol = 5e-3
        @test eltype(r) == T

        @test v[1] / 1000 ≈ +1.656142 atol = 1e-5
        @test v[2] / 1000 ≈ -6.699518 atol = 1e-5
        @test v[3] / 1000 ≈ -1.266334 atol = 1e-5
        @test eltype(v) == T

        # Test simultaneous initialization and propagation.
        r, v, orbp = Propagators.propagate(Val(:J2), (jd₁ - jd₀) * 86400, orb)

        orbk = Propagators.mean_elements(orbp)
        @test orbk isa KeplerianElements{Float64, Float64}

        @test r[1] / 1000 ≈ -6849.654348 atol = 5e-3
        @test r[2] / 1000 ≈ -2253.059809 atol = 5e-3
        @test r[3] / 1000 ≈ +3574.529667 atol = 5e-3
        @test eltype(r) == T

        @test v[1] / 1000 ≈ +1.656142 atol = 1e-5
        @test v[2] / 1000 ≈ -6.699518 atol = 1e-5
        @test v[3] / 1000 ≈ -1.266334 atol = 1e-5
        @test eltype(v) == T

        r, v, orbp = Propagators.propagate_to_epoch(Val(:J2), jd₁, orb)

        orbk = Propagators.mean_elements(orbp)
        @test orbk isa KeplerianElements{Float64, Float64}

        @test r[1] / 1000 ≈ -6849.654348 atol = 5e-3
        @test r[2] / 1000 ≈ -2253.059809 atol = 5e-3
        @test r[3] / 1000 ≈ +3574.529667 atol = 5e-3
        @test eltype(r) == T

        @test v[1] / 1000 ≈ +1.656142 atol = 1e-5
        @test v[2] / 1000 ≈ -6.699518 atol = 1e-5
        @test v[3] / 1000 ≈ -1.266334 atol = 1e-5
        @test eltype(v) == T

        r, v, j2d = j2((jd₁ - jd₀) * 86400, orb)

        @test j2d isa J2Propagator{Float64, Float64}

        @test r[1] / 1000 ≈ -6849.654348 atol = 5e-3
        @test r[2] / 1000 ≈ -2253.059809 atol = 5e-3
        @test r[3] / 1000 ≈ +3574.529667 atol = 5e-3
        @test eltype(r) == T

        @test v[1] / 1000 ≈ +1.656142 atol = 1e-5
        @test v[2] / 1000 ≈ -6.699518 atol = 1e-5
        @test v[3] / 1000 ≈ -1.266334 atol = 1e-5
        @test eltype(v) == T
    end

    # == Float32 ===========================================================================

    @testset "Float32" verbose = true begin
        T = Float32

        orb = KeplerianElements(
            jd₀,
            T(8000e3),
            T(0.015),
            T(28.5) |> deg2rad,
            T(100)  |> deg2rad,
            T(200)  |> deg2rad,
            T(45)   |> deg2rad
        )

        orbp = Propagators.init(Val(:J2), orb; j2c = j2c_egm2008_f32)

        orbk = Propagators.mean_elements(orbp)
        @test orbk isa KeplerianElements{Float64, Float32}
        @test orbk.t ≈ orb.t
        @test orbk.a ≈ orb.a
        @test orbk.e ≈ orb.e
        @test orbk.i ≈ orb.i
        @test orbk.Ω ≈ orb.Ω
        @test orbk.ω ≈ orb.ω
        @test orbk.f ≈ orb.f

        r, v = Propagators.step!(orbp, (jd₁ - jd₀) * 86400)

        @test r[1] / 1000 ≈ -6849.654348 atol = 5e-1
        @test r[2] / 1000 ≈ -2253.059809 atol = 5e-1
        @test r[3] / 1000 ≈ +3574.529667 atol = 5e-1
        @test eltype(r) == T

        @test v[1] / 1000 ≈ +1.656142 atol = 1e-3
        @test v[2] / 1000 ≈ -6.699518 atol = 1e-3
        @test v[3] / 1000 ≈ -1.266334 atol = 1e-3
        @test eltype(v) == T

        r, v = Propagators.propagate_to_epoch!(orbp, jd₁)

        @test r[1] / 1000 ≈ -6849.654348 atol = 5e-1
        @test r[2] / 1000 ≈ -2253.059809 atol = 5e-1
        @test r[3] / 1000 ≈ +3574.529667 atol = 5e-1
        @test eltype(r) == T

        @test v[1] / 1000 ≈ +1.656142 atol = 1e-3
        @test v[2] / 1000 ≈ -6.699518 atol = 1e-3
        @test v[3] / 1000 ≈ -1.266334 atol = 1e-3
        @test eltype(v) == T

        # Test in-place initialization.
        orbp = OrbitPropagatorJ2(J2Propagator{Float64, T}())
        orbp.j2d.j2c = j2c_egm2008_f32
        Propagators.init!(orbp, orb)

        r, v = Propagators.step!(orbp, (jd₁ - jd₀) * 86400)

        @test r[1] / 1000 ≈ -6849.654348 atol = 5e-1
        @test r[2] / 1000 ≈ -2253.059809 atol = 5e-1
        @test r[3] / 1000 ≈ +3574.529667 atol = 5e-1
        @test eltype(r) == T

        @test v[1] / 1000 ≈ +1.656142 atol = 1e-3
        @test v[2] / 1000 ≈ -6.699518 atol = 1e-3
        @test v[3] / 1000 ≈ -1.266334 atol = 1e-3
        @test eltype(v) == T

        # Test simultaneous initialization and propagation.
        r, v, orbp = Propagators.propagate(Val(:J2), (jd₁ - jd₀) * 86400, orb; j2c = j2c_egm2008_f32)

        orbk = Propagators.mean_elements(orbp)
        @test orbk isa KeplerianElements{Float64, Float32}

        @test r[1] / 1000 ≈ -6849.654348 atol = 5e-1
        @test r[2] / 1000 ≈ -2253.059809 atol = 5e-1
        @test r[3] / 1000 ≈ +3574.529667 atol = 5e-1
        @test eltype(r) == T

        @test v[1] / 1000 ≈ +1.656142 atol = 1e-3
        @test v[2] / 1000 ≈ -6.699518 atol = 1e-3
        @test v[3] / 1000 ≈ -1.266334 atol = 1e-3
        @test eltype(v) == T

        r, v, orbp = Propagators.propagate_to_epoch(Val(:J2), jd₁, orb; j2c = j2c_egm2008_f32)

        orbk = Propagators.mean_elements(orbp)
        @test orbk isa KeplerianElements{Float64, Float32}

        @test r[1] / 1000 ≈ -6849.654348 atol = 5e-1
        @test r[2] / 1000 ≈ -2253.059809 atol = 5e-1
        @test r[3] / 1000 ≈ +3574.529667 atol = 5e-1
        @test eltype(r) == T

        @test v[1] / 1000 ≈ +1.656142 atol = 1e-3
        @test v[2] / 1000 ≈ -6.699518 atol = 1e-3
        @test v[3] / 1000 ≈ -1.266334 atol = 1e-3
        @test eltype(v) == T

        r, v, j2d = j2((jd₁ - jd₀) * 86400, orb; j2c = j2c_egm2008_f32)

        @test j2d isa J2Propagator{Float64, Float32}

        @test r[1] / 1000 ≈ -6849.654348 atol = 5e-1
        @test r[2] / 1000 ≈ -2253.059809 atol = 5e-1
        @test r[3] / 1000 ≈ +3574.529667 atol = 5e-1
        @test eltype(r) == T

        @test v[1] / 1000 ≈ +1.656142 atol = 1e-3
        @test v[2] / 1000 ≈ -6.699518 atol = 1e-3
        @test v[3] / 1000 ≈ -1.266334 atol = 1e-3
        @test eltype(v) == T
    end
end

@testset "Fitting Mean Elements for the J2 Propagator" verbose = true begin
    # We just need to run the J2 propagator, obtain the stave vector, convert them to mean
    # Keplerian elements, and compare with the original mean elements.

    orb_input = KeplerianElements(
        DateTime("2023-01-01") |> datetime2julian,
        7130.982e3,
        0.001111,
        98.405 |> deg2rad,
        90     |> deg2rad,
        200    |> deg2rad,
        45     |> deg2rad
    )

    # Generate the osculating elements.
    orbp = Propagators.init(Val(:J2), orb_input)
    ret  = Propagators.propagate!.(orbp, 0:10:12_000)
    vr_i = first.(ret)
    vv_i = last.(ret)
    vjd  = Propagators.epoch(orbp) .+ (0:10:12_000) ./ 86400

    @testset "Without Initial Guess" begin
        # Obtain the mean elements.
        orb, ~ = redirect_stdout(devnull) do
            Propagators.fit_mean_elements(
                Val(:J2),
                vjd,
                vr_i,
                vv_i;
                mean_elements_epoch = vjd[begin],
            )
        end

        @test orb.t ≈ orb_input.t
        @test orb.a ≈ orb_input.a
        @test orb.e ≈ orb_input.e
        @test orb.i ≈ orb_input.i
        @test orb.Ω ≈ orb_input.Ω
        @test orb.ω ≈ orb_input.ω
        @test orb.f ≈ orb_input.f atol = 1e-6

        # Obtain the mean elements.
        orb, ~ = redirect_stdout(devnull) do
            Propagators.fit_mean_elements!(
                orbp,
                vjd,
                vr_i,
                vv_i;
                mean_elements_epoch = vjd[begin],
            )
        end

        @test orb.t ≈ orb_input.t
        @test orb.a ≈ orb_input.a
        @test orb.e ≈ orb_input.e
        @test orb.i ≈ orb_input.i
        @test orb.Ω ≈ orb_input.Ω
        @test orb.ω ≈ orb_input.ω
        @test orb.f ≈ orb_input.f atol = 1e-6

        # Test with very low perturbation in the Jacobian.
        orb, ~ = redirect_stdout(devnull) do
            Propagators.fit_mean_elements(
                Val(:J2),
                vjd,
                vr_i,
                vv_i;
                mean_elements_epoch = vjd[begin],
                jacobian_perturbation = 1e-13
            )
        end

        @test orb.t ≈ orb_input.t
        @test orb.a ≈ orb_input.a
        @test orb.e ≈ orb_input.e
        @test orb.i ≈ orb_input.i
        @test orb.Ω ≈ orb_input.Ω
        @test orb.ω ≈ orb_input.ω
        @test orb.f ≈ orb_input.f atol = 1e-6

        orb, ~ = redirect_stdout(devnull) do
            Propagators.fit_mean_elements!(
                orbp,
                vjd,
                vr_i,
                vv_i;
                mean_elements_epoch = vjd[begin],
                jacobian_perturbation = 1e-13
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

    @testset "With Initial Guess" begin
        # Obtain the mean elements.
        orb, ~ = redirect_stdout(devnull) do
            Propagators.fit_mean_elements(
                Val(:J2),
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

        # Obtain the mean elements.
        orb, ~ = redirect_stdout(devnull) do
            Propagators.fit_mean_elements!(
                orbp,
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
        # Obtain the mean elements.
        orb, ~ = redirect_stdout(devnull) do
            Propagators.fit_mean_elements(
                Val(:J2),
                vjd,
                vr_i,
                vv_i;
                mean_elements_epoch = vjd[begin] + 1,
            )
        end

        @test orb.t ≈ orb_input.t + 1
        @test orb.a ≈ orb_input.a
        @test orb.e ≈ orb_input.e
        @test orb.i ≈ orb_input.i

        # The input orbit is the nominal orbit for the Amazonia-1 satellite, which is Sun
        # synchronous. Thus, the RAAN moves approximately 0.9856002605° per day.
        @test orb.Ω ≈ orb_input.Ω + deg2rad(0.9856002605) atol = 2e-5

        # Obtain the mean elements.
        orb, ~ = redirect_stdout(devnull) do
            Propagators.fit_mean_elements!(
                orbp,
                vjd,
                vr_i,
                vv_i;
                mean_elements_epoch = vjd[begin] + 1,
            )
        end

        @test orb.t ≈ orb_input.t + 1
        @test orb.a ≈ orb_input.a
        @test orb.e ≈ orb_input.e
        @test orb.i ≈ orb_input.i

        # The input orbit is the nominal orbit for the Amazonia-1 satellite, which is Sun
        # synchronous. Thus, the RAAN moves approximately 0.9856002605° per day.
        @test orb.Ω ≈ orb_input.Ω + deg2rad(0.9856002605) atol = 2e-5
    end

    @testset "Errors" begin
        # == Wrong dimensions in the input vectors =========================================

        @test_throws ArgumentError Propagators.fit_mean_elements(Val(:J2), vjd[1:end-1], vr_i, vv_i)
        @test_throws ArgumentError Propagators.fit_mean_elements(Val(:J2), vjd, vr_i[1:end-1], vv_i)
        @test_throws ArgumentError Propagators.fit_mean_elements(Val(:J2), vjd, vr_i, vv_i[1:end-1])

        # == Wrong dimensions in the weight vector =========================================

        @test_throws ArgumentError Propagators.fit_mean_elements(
            Val(:J2),
            vjd,
            vr_i,
            vv_i;
            weight_vector = [1, 2, 3, 4, 5]
        )
    end
end

@testset "Update J2 Mean Elements Epoch" begin
    orb_input = KeplerianElements(
        DateTime("2023-01-01") |> datetime2julian,
        7130.982e3,
        0.001111,
        98.405 |> deg2rad,
        90     |> deg2rad,
        200    |> deg2rad,
        45     |> deg2rad
    )

    orb = update_j2_mean_elements_epoch(orb_input, DateTime("2023-01-02"))

    @test orb.t ≈ orb_input.t + 1
    @test orb.a == orb_input.a
    @test orb.e == orb_input.e
    @test orb.i == orb_input.i

    # The input orbit is the nominal orbit for the Amazonia-1 satellite, which is Sun
    # synchronous. Thus, the RAAN moves approximately 0.9856002605° per day.
    @test orb.Ω ≈ orb_input.Ω + deg2rad(0.9856002605) atol = 2e-5
end

@testset "Copying Structure" verbose = true begin
    for (T, j2c) in ((Float64, j2c_egm2008), (Float32, j2c_egm2008_f32))
        @testset "$T" begin
            jd₀ = date_to_jd(2023, 1, 1, 0, 0, 0)

            orb = KeplerianElements(
                jd₀,
                T(8000e3),
                T(0.015),
                T(28.5) |> deg2rad,
                T(100)  |> deg2rad,
                T(200)  |> deg2rad,
                T(45)   |> deg2rad
            )

            orbp = Propagators.init(Val(:J2), orb; j2c = j2c)
            new_orbp = copy(orbp)

            for f in fieldnames(typeof(orbp.j2d))
                @test getfield(new_orbp.j2d, f) == getfield(orbp.j2d, f)
            end

            new_orbp.j2d.Δt = 1000
            @test new_orbp.j2d.Δt != orbp.j2d.Δt
        end
    end
end
