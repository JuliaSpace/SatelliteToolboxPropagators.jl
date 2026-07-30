## Description #############################################################################
#
#   Tests related to the J2 osculating orbit propagator.
#
############################################################################################

############################################################################################
#                                       Test Results                                       #
############################################################################################
#
# The following test results were obtained using a numerical propagator as shown in:
#
#   https://github.com/JuliaSpace/PropagatorTests/
#
# == Initial Mean Elements (TOD) ===========================================================
#
#            Epoch :    2.45995e6 (2023-01-01T00:00:00)
#  Semi-major axis : 7190.98     km
#     Eccentricity :    0.001111
#      Inclination :   98.405    °
#             RAAN :   90.0      °
#  Arg. of Perigee :  200.0      °
#     True Anomaly :   45.0      °
#
# == Numerical Propagation Results =========================================================
#
# Gravity model         : EGM-2008 (Degree 4, Order 0)
# Integration algorithm : AutoVern7(Rodas5())
#
#  Time [s] │ Pos. X (TOD)  Pos. Y (TOD)  Pos. Z (TOD)  Vel. X (TOD)  Vel. Y (TOD)  Vel. Z (TOD)
#           │           km            km            km        km / s        km / s        km / s
# ──────────┼────────────────────────────────────────────────────────────────────────────────────
#       0.0 │     -952.883     -3038.438     -6444.903        -0.460         6.745        -3.116
#     600.0 │    -1034.273      1320.078     -6994.339         0.197         7.315         1.342
#    1200.0 │     -731.130      5188.074     -4936.895         0.781         5.164         5.295
#    1800.0 │     -156.143      7126.675     -1039.300         1.074         1.090         7.277
#    2400.0 │      476.947      6413.653      3245.983         0.968        -3.389         6.546
#    3000.0 │      932.852      3316.764      6322.415         0.503        -6.601         3.379
#    3600.0 │     1042.661     -1010.960      7047.367        -0.149        -7.361        -1.041
#    4200.0 │      765.306     -4962.529      5148.982        -0.747        -5.385        -5.085
#    4800.0 │      202.796     -7063.290      1327.333        -1.068        -1.388        -7.242
#    5400.0 │     -435.494     -6521.701     -2991.520        -0.991         3.135        -6.686
#    6000.0 │     -910.991     -3540.650     -6189.318        -0.543         6.478        -3.629
#
############################################################################################

@testset "J2 Osculating Orbit Propagator" verbose = true begin
    # Create a matrix with the results.
    #
    # Notice that the results here are not exactly the same because in J2 osculating
    # propagator the long-term periodics are not taken into account. Hence, the tolerances
    # below (2 km in position, 3 m / s in velocity) reflect this modeling difference against
    # the numerical propagator, not the accuracy of the implementation. The short-period
    # correction itself is checked to a much tighter tolerance in the testset
    # "Osculating Radial Rate".
    results = [
           0.0  -952.883 -3038.438 -6444.903 -0.460  6.745 -3.116
         600.0 -1034.273  1320.078 -6994.339  0.197  7.315  1.342
        1200.0  -731.130  5188.074 -4936.895  0.781  5.164  5.295
        1800.0  -156.143  7126.675 -1039.300  1.074  1.090  7.277
        2400.0   476.947  6413.653  3245.983  0.968 -3.389  6.546
        3000.0   932.852  3316.764  6322.415  0.503 -6.601  3.379
        3600.0  1042.661 -1010.960  7047.367 -0.149 -7.361 -1.041
        4200.0   765.306 -4962.529  5148.982 -0.747 -5.385 -5.085
        4800.0   202.796 -7063.290  1327.333 -1.068 -1.388 -7.242
        5400.0  -435.494 -6521.701 -2991.520 -0.991  3.135 -6.686
        6000.0  -910.991 -3540.650 -6189.318 -0.543  6.478 -3.629
    ]

    # == Constructor =======================================================================

    @testset "Constructor" begin
        orb    = KeplerianElements(0.0, 8000.0e3, 0.0, 0.0, 0.0, 0.0, 0.0)
        j2d    = J2Propagator{Float64, Float64}(orb, orb, j2c_egm2008, 0, 0, 0, 0, 0, 0)
        j2oscd = J2OsculatingPropagator{Float64, Float64}(j2d, 0, orb)

        # Test some random fields.
        @test j2oscd.j2d == j2d
        @test j2oscd.Δt == 0
        @test j2d.orbk == orb

        orb = KeplerianElements(0.0, 8_000_000, 0, 0, 0, 0, 0)
        orbp = Propagators.init(Val(:J2osc), orb)
        @test orbp.j2oscd.j2d.orb₀ isa KeplerianElements{Float64, Float64}
    end

    # == General API Functions =============================================================

    @testset "General API Functions" begin
        orb = KeplerianElements(0.0, 8000.0e3, 0.0, 0.0, 0.0, 0.0, 0.0)
        orbp = Propagators.init(Val(:J2osc), orb)
        @test Propagators.name(orbp) == "J2 Osculating Orbit Propagator"
    end

    # == Float64 ===========================================================================

    @testset "Float64" begin
        T = Float64

        # Initialize the propagator.
        jd₀ = date_to_jd(2023, 1, 1, 0, 0, 0)
        orb = KeplerianElements(
            jd₀,
            T(7190.982e3),
            T(0.001111),
            T(98.405) |> deg2rad,
            T(90) |> deg2rad,
            T(200) |> deg2rad,
            T(45) |> deg2rad,
        )

        # Test all the results.
        orbp = Propagators.init(Val(:J2osc), orb; j2c = j2c_egm2008)

        for k in axes(results, 1)
            r, v = Propagators.propagate!(orbp, results[k, 1])

            @test results[k, 2] ≈ r[1] / 1000 atol = 2
            @test results[k, 3] ≈ r[2] / 1000 atol = 2
            @test results[k, 4] ≈ r[3] / 1000 atol = 2
            @test results[k, 5] ≈ v[1] / 1000 atol = 3e-3
            @test results[k, 6] ≈ v[2] / 1000 atol = 3e-3
            @test results[k, 7] ≈ v[3] / 1000 atol = 3e-3
            @test eltype(r) == T
            @test eltype(v) == T
        end

        # Re-initialize the propagator.
        orbp = Propagators.init(Val(:J2osc), orb; j2c = j2c_egm2008)

        orbk = Propagators.mean_elements(orbp)
        @test orbk isa KeplerianElements{Float64, Float64}
        @test orbk.t ≈ orb.t
        @test orbk.a ≈ orb.a
        @test orbk.e ≈ orb.e
        @test orbk.i ≈ orb.i
        @test orbk.Ω ≈ orb.Ω
        @test orbk.ω ≈ orb.ω
        @test orbk.f ≈ orb.f

        r, v = Propagators.step!(orbp, results[end, 1])

        @test results[end, 2] ≈ r[1] / 1000 atol = 2e-1
        @test results[end, 3] ≈ r[2] / 1000 atol = 2e-1
        @test results[end, 4] ≈ r[3] / 1000 atol = 2e-1
        @test results[end, 5] ≈ v[1] / 1000 atol = 1e-3
        @test results[end, 6] ≈ v[2] / 1000 atol = 1e-3
        @test results[end, 7] ≈ v[3] / 1000 atol = 1e-3
        @test eltype(r) == T
        @test eltype(v) == T

        r, v = Propagators.propagate_to_epoch!(orbp, jd₀ + results[end, 1] / 86400)

        @test results[end, 2] ≈ r[1] / 1000 atol = 2e-1
        @test results[end, 3] ≈ r[2] / 1000 atol = 2e-1
        @test results[end, 4] ≈ r[3] / 1000 atol = 2e-1
        @test results[end, 5] ≈ v[1] / 1000 atol = 1e-3
        @test results[end, 6] ≈ v[2] / 1000 atol = 1e-3
        @test results[end, 7] ≈ v[3] / 1000 atol = 1e-3
        @test eltype(r) == T
        @test eltype(v) == T

        # Test in-place initialization.
        orbp = OrbitPropagatorJ2Osculating(J2OsculatingPropagator{Float64, T}())
        j2d = J2Propagator{Float64, T}()
        j2d.j2c = j2c_egm2008
        orbp.j2oscd.j2d = j2d
        Propagators.init!(orbp, orb)

        r, v = Propagators.step!(orbp, results[end, 1])

        @test results[end, 2] ≈ r[1] / 1000 atol = 2e-1
        @test results[end, 3] ≈ r[2] / 1000 atol = 2e-1
        @test results[end, 4] ≈ r[3] / 1000 atol = 2e-1
        @test results[end, 5] ≈ v[1] / 1000 atol = 1e-3
        @test results[end, 6] ≈ v[2] / 1000 atol = 1e-3
        @test results[end, 7] ≈ v[3] / 1000 atol = 1e-3
        @test eltype(r) == T
        @test eltype(v) == T

        # Test simultaneous initialization and propagation.
        r, v, orbp = Propagators.propagate(Val(:J2osc), results[end, 1], orb)

        orbk = Propagators.mean_elements(orbp)
        @test orbk isa KeplerianElements{Float64, Float64}

        @test results[end, 2] ≈ r[1] / 1000 atol = 2e-1
        @test results[end, 3] ≈ r[2] / 1000 atol = 2e-1
        @test results[end, 4] ≈ r[3] / 1000 atol = 2e-1
        @test results[end, 5] ≈ v[1] / 1000 atol = 1e-3
        @test results[end, 6] ≈ v[2] / 1000 atol = 1e-3
        @test results[end, 7] ≈ v[3] / 1000 atol = 1e-3
        @test eltype(r) == T
        @test eltype(v) == T

        r, v, orbp = Propagators.propagate_to_epoch(
            Val(:J2osc), jd₀ + results[end, 1] / 86400, orb
        )

        orbk = Propagators.mean_elements(orbp)
        @test orbk isa KeplerianElements{Float64, Float64}

        @test results[end, 2] ≈ r[1] / 1000 atol = 2e-1
        @test results[end, 3] ≈ r[2] / 1000 atol = 2e-1
        @test results[end, 4] ≈ r[3] / 1000 atol = 2e-1
        @test results[end, 5] ≈ v[1] / 1000 atol = 1e-3
        @test results[end, 6] ≈ v[2] / 1000 atol = 1e-3
        @test results[end, 7] ≈ v[3] / 1000 atol = 1e-3
        @test eltype(r) == T
        @test eltype(v) == T

        r, v, j2oscd = j2osc(results[end, 1], orb)

        @test j2oscd isa J2OsculatingPropagator{Float64, Float64}

        @test results[end, 2] ≈ r[1] / 1000 atol = 2e-1
        @test results[end, 3] ≈ r[2] / 1000 atol = 2e-1
        @test results[end, 4] ≈ r[3] / 1000 atol = 2e-1
        @test results[end, 5] ≈ v[1] / 1000 atol = 1e-3
        @test results[end, 6] ≈ v[2] / 1000 atol = 1e-3
        @test results[end, 7] ≈ v[3] / 1000 atol = 1e-3
        @test eltype(r) == T
        @test eltype(v) == T
    end

    # == Float32 ===========================================================================

    @testset "Float32" begin
        T = Float32

        # Initialize the propagator.
        jd₀ = date_to_jd(2023, 1, 1, 0, 0, 0)
        orb = KeplerianElements(
            jd₀,
            T(7190.982e3),
            T(0.001111),
            T(98.405) |> deg2rad,
            T(90) |> deg2rad,
            T(200) |> deg2rad,
            T(45) |> deg2rad,
        )

        # Test all the results.
        orbp = Propagators.init(Val(:J2osc), orb; j2c = j2c_egm2008_f32)

        for k in axes(results, 1)
            r, v = Propagators.propagate!(orbp, results[k, 1])

            @test results[k, 2] ≈ r[1] / 1000 atol = 2
            @test results[k, 3] ≈ r[2] / 1000 atol = 2
            @test results[k, 4] ≈ r[3] / 1000 atol = 2
            @test results[k, 5] ≈ v[1] / 1000 atol = 3e-3
            @test results[k, 6] ≈ v[2] / 1000 atol = 3e-3
            @test results[k, 7] ≈ v[3] / 1000 atol = 3e-3
            @test eltype(r) == T
            @test eltype(v) == T
        end

        # Re-initialize the propagator.
        orbp = Propagators.init(Val(:J2osc), orb; j2c = j2c_egm2008_f32)

        orbk = Propagators.mean_elements(orbp)
        @test orbk isa KeplerianElements{Float64, Float32}
        @test orbk.t ≈ orb.t
        @test orbk.a ≈ orb.a
        @test orbk.e ≈ orb.e
        @test orbk.i ≈ orb.i
        @test orbk.Ω ≈ orb.Ω
        @test orbk.ω ≈ orb.ω
        @test orbk.f ≈ orb.f

        r, v = Propagators.step!(orbp, results[end, 1])

        @test results[end, 2] ≈ r[1] / 1000 atol = 2e-1
        @test results[end, 3] ≈ r[2] / 1000 atol = 2e-1
        @test results[end, 4] ≈ r[3] / 1000 atol = 2e-1
        @test results[end, 5] ≈ v[1] / 1000 atol = 1e-3
        @test results[end, 6] ≈ v[2] / 1000 atol = 1e-3
        @test results[end, 7] ≈ v[3] / 1000 atol = 1e-3
        @test eltype(r) == T
        @test eltype(v) == T

        r, v = Propagators.propagate_to_epoch!(orbp, jd₀ + results[end, 1] / 86400)

        @test results[end, 2] ≈ r[1] / 1000 atol = 2e-1
        @test results[end, 3] ≈ r[2] / 1000 atol = 2e-1
        @test results[end, 4] ≈ r[3] / 1000 atol = 2e-1
        @test results[end, 5] ≈ v[1] / 1000 atol = 1e-3
        @test results[end, 6] ≈ v[2] / 1000 atol = 1e-3
        @test results[end, 7] ≈ v[3] / 1000 atol = 1e-3
        @test eltype(r) == T
        @test eltype(v) == T

        # Test in-place initialization.
        orbp = OrbitPropagatorJ2Osculating(J2OsculatingPropagator{Float64, T}())
        j2d = J2Propagator{Float64, T}()
        j2d.j2c = j2c_egm2008_f32
        orbp.j2oscd.j2d = j2d
        Propagators.init!(orbp, orb)

        r, v = Propagators.step!(orbp, results[end, 1])

        @test results[end, 2] ≈ r[1] / 1000 atol = 2e-1
        @test results[end, 3] ≈ r[2] / 1000 atol = 2e-1
        @test results[end, 4] ≈ r[3] / 1000 atol = 2e-1
        @test results[end, 5] ≈ v[1] / 1000 atol = 1e-3
        @test results[end, 6] ≈ v[2] / 1000 atol = 1e-3
        @test results[end, 7] ≈ v[3] / 1000 atol = 1e-3
        @test eltype(r) == T
        @test eltype(v) == T

        # Test simultaneous initialization and propagation.
        r, v, orbp = Propagators.propagate(
            Val(:J2osc), results[end, 1], orb; j2c = j2c_egm2008_f32
        )

        orbk = Propagators.mean_elements(orbp)
        @test orbk isa KeplerianElements{Float64, Float32}

        @test results[end, 2] ≈ r[1] / 1000 atol = 2e-1
        @test results[end, 3] ≈ r[2] / 1000 atol = 2e-1
        @test results[end, 4] ≈ r[3] / 1000 atol = 2e-1
        @test results[end, 5] ≈ v[1] / 1000 atol = 1e-3
        @test results[end, 6] ≈ v[2] / 1000 atol = 1e-3
        @test results[end, 7] ≈ v[3] / 1000 atol = 1e-3
        @test eltype(r) == T
        @test eltype(v) == T

        r, v, orbp = Propagators.propagate_to_epoch(
            Val(:J2osc), jd₀ + results[end, 1] / 86400, orb; j2c = j2c_egm2008_f32
        )

        orbk = Propagators.mean_elements(orbp)
        @test orbk isa KeplerianElements{Float64, Float32}

        @test results[end, 2] ≈ r[1] / 1000 atol = 2e-1
        @test results[end, 3] ≈ r[2] / 1000 atol = 2e-1
        @test results[end, 4] ≈ r[3] / 1000 atol = 2e-1
        @test results[end, 5] ≈ v[1] / 1000 atol = 1e-3
        @test results[end, 6] ≈ v[2] / 1000 atol = 1e-3
        @test results[end, 7] ≈ v[3] / 1000 atol = 1e-3
        @test eltype(r) == T
        @test eltype(v) == T

        r, v, j2oscd = j2osc(results[end, 1], orb; j2c = j2c_egm2008_f32)

        @test j2oscd isa J2OsculatingPropagator{Float64, Float32}

        @test results[end, 2] ≈ r[1] / 1000 atol = 2e-1
        @test results[end, 3] ≈ r[2] / 1000 atol = 2e-1
        @test results[end, 4] ≈ r[3] / 1000 atol = 2e-1
        @test results[end, 5] ≈ v[1] / 1000 atol = 1e-3
        @test results[end, 6] ≈ v[2] / 1000 atol = 1e-3
        @test results[end, 7] ≈ v[3] / 1000 atol = 1e-3
        @test eltype(r) == T
        @test eltype(v) == T
    end
end

@testset "Fitting Mean Elements for the J2 Osculating Orbit Propagator" verbose = true begin
    # We just need to run the J2 osculating propagator, obtain the osculating elements,
    # convert them to mean Keplerian elements, and compare with the original mean elements.

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
    orbp = Propagators.init(Val(:J2osc), orb_input)
    ret  = Propagators.propagate!.(orbp, 0:10:12_000)
    vr_i = first.(ret)
    vv_i = last.(ret)
    vjd  = Propagators.epoch(orbp) .+ (0:10:12_000) ./ 86400

    @testset "Without Initial Guess" begin
        # Obtain the mean elements.
        orb, ~ = redirect_stdout(devnull) do
            Propagators.fit_mean_elements(
                Val(:J2osc), vjd, vr_i, vv_i; mean_elements_epoch = vjd[begin]
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
                orbp, vjd, vr_i, vv_i; mean_elements_epoch = vjd[begin]
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
                Val(:J2osc),
                vjd,
                vr_i,
                vv_i;
                mean_elements_epoch = vjd[begin],
                jacobian_perturbation = 1e-13,
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
                jacobian_perturbation = 1e-13,
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
                Val(:J2osc),
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
                Val(:J2osc), vjd, vr_i, vv_i; mean_elements_epoch = vjd[begin] + 1
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
                orbp, vjd, vr_i, vv_i; mean_elements_epoch = vjd[begin] + 1
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

        @test_throws ArgumentError Propagators.fit_mean_elements(
            Val(:J2osc), vjd[1:(end - 1)], vr_i, vv_i
        )
        @test_throws ArgumentError Propagators.fit_mean_elements(
            Val(:J2osc), vjd, vr_i[1:(end - 1)], vv_i
        )
        @test_throws ArgumentError Propagators.fit_mean_elements(
            Val(:J2osc), vjd, vr_i, vv_i[1:(end - 1)]
        )

        # == Wrong dimensions in the weight vector =========================================

        @test_throws ArgumentError Propagators.fit_mean_elements(
            Val(:J2osc), vjd, vr_i, vv_i; weight_vector = [1, 2, 3, 4, 5]
        )
    end
end

@testset "Fitting Mean Elements for the J2 Osculating Propagator (ForwardDiffJacobian)" verbose =
    true begin
    orb_input = KeplerianElements(
        DateTime("2023-01-01") |> datetime2julian,
        7130.982e3,
        0.001111,
        98.405 |> deg2rad,
        90 |> deg2rad,
        200 |> deg2rad,
        45 |> deg2rad,
    )

    orbp = Propagators.init(Val(:J2osc), orb_input)
    ret  = Propagators.propagate!.(orbp, 0:10:12_000)
    vr_i = first.(ret)
    vv_i = last.(ret)
    vjd  = Propagators.epoch(orbp) .+ (0:10:12_000) ./ 86400

    @testset "Without Initial Guess" begin
        orb, P = redirect_stdout(devnull) do
            Propagators.fit_mean_elements(
                Val(:J2osc),
                vjd,
                vr_i,
                vv_i;
                mean_elements_epoch = vjd[begin],
                jacobian_method = ForwardDiffJacobian(),
            )
        end

        @test orb.t ≈ orb_input.t
        @test orb.a ≈ orb_input.a atol = 100
        @test orb.e ≈ orb_input.e atol = 1e-3
        @test orb.i ≈ orb_input.i atol = 1e-5
        @test P isa SMatrix{6, 6, Float64}
    end

    @testset "With Initial Guess" begin
        orb, P = redirect_stdout(devnull) do
            Propagators.fit_mean_elements(
                Val(:J2osc),
                vjd,
                vr_i,
                vv_i;
                max_iterations = 3,
                mean_elements_epoch = vjd[begin],
                initial_guess = orb_input,
                jacobian_method = ForwardDiffJacobian(),
            )
        end

        @test orb.t ≈ orb_input.t
        @test orb.a ≈ orb_input.a atol = 100
        @test orb.e ≈ orb_input.e atol = 1e-3
        @test orb.i ≈ orb_input.i atol = 1e-5
    end
end

@testset "Update J2 Osculating Mean Elements Epoch" begin
    orb_input = KeplerianElements(
        DateTime("2023-01-01") |> datetime2julian,
        7130.982e3,
        0.001111,
        98.405 |> deg2rad,
        90 |> deg2rad,
        200 |> deg2rad,
        45 |> deg2rad,
    )

    orb = update_j2osc_mean_elements_epoch(orb_input, DateTime("2023-01-02"))

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
                T(100) |> deg2rad,
                T(200) |> deg2rad,
                T(25) |> deg2rad,
            )

            orbp = Propagators.init(Val(:J2osc), orb; j2c = j2c)
            new_orbp = copy(orbp)

            for f in fieldnames(typeof(orbp.j2oscd))
                f == :j2d && continue
                @test getfield(new_orbp.j2oscd, f) == getfield(orbp.j2oscd, f)
            end

            for f in fieldnames(typeof(orbp.j2oscd.j2d))
                @test getfield(new_orbp.j2oscd.j2d, f) == getfield(orbp.j2oscd.j2d, f)
            end

            new_orbp.j2oscd.Δt = 1000
            @test new_orbp.j2oscd.Δt != orbp.j2oscd.Δt

            new_orbp.j2oscd.j2d.Δt = 1000
            @test new_orbp.j2oscd.j2d.Δt != orbp.j2oscd.j2d.Δt
        end
    end
end

# == Osculating Radial Rate ================================================================
#
# The short-period correction to the radial rate is not visible in the reference tables
# above, because those come from a numerical propagator with a full gravity model and differ
# from these propagators by about 2 km in position and 2 m / s in velocity regardless.
#
# The reference values below isolate it. They were obtained by integrating the J2 equations
# of motion (two-body plus the J2 zonal acceleration, EGM-2008 constants) with a fixed-step
# RK4 integrator using 2000 steps per orbit, and extracting the mean elements from the
# resulting trajectory with a Hann-windowed average over six orbital periods, which is the
# Brouwer definition of mean elements to first order. Each row holds the extracted mean
# elements and the true osculating radial rate, r ⋅ v / |r|, at the centre of the window:
#
#     (a [m], e, i [rad], ω [rad], M [rad], ṙ [m / s])
#
# The inclination is 98.405°, the RAAN is 90°, and the eccentricities are 0.02, 0.08, and
# 0.15, each sampled at eight anomalies.
#
# The tolerance is the residue of the first-order theory. The expression previously
# implemented, which used (1 - e cos f)² instead of (1 + e cos f)² in the second term, is off
# by up to 1.33 m / s and fails this test.

@testset "Osculating Radial Rate" verbose = true begin
    reference = [
        (
            7183847.8491930151,
            0.0198000650,
            1.7175640808,
            3.4760929577,
            0.0139527782,
            0.0000000006,
        ),
        (
            7197161.6035945583,
            0.0211832948,
            1.7174314741,
            3.5221560837,
            0.7252183073,
            105.3117091008,
        ),
        (
            7197901.0608316371,
            0.0198150794,
            1.7174197442,
            3.5614404178,
            1.4606498191,
            148.9332472863,
        ),
        (
            7185324.8573491490,
            0.0198369298,
            1.7175489941,
            3.4773432427,
            2.3416455678,
            105.3117090999,
        ),
        (
            7184330.6082761073,
            0.0201996025,
            1.7175602961,
            3.5029620178,
            3.1287134373,
            -0.0000000007,
        ),
        (
            7196425.4519446185,
            0.0189399367,
            1.7174323769,
            3.4549929272,
            3.9905225353,
            -105.3117091007,
        ),
        (
            7197900.9734012149,
            0.0203497524,
            1.7174213290,
            3.4221276352,
            4.8214724651,
            -148.9332472862,
        ),
        (
            7185077.4157844847,
            0.0201873767,
            1.7175525769,
            3.5049500022,
            5.5121762531,
            -105.3117090998,
        ),
        (
            7182931.4258067394,
            0.0797237059,
            1.7175708190,
            3.4864859873,
            0.0034948928,
            0.0000000038,
        ),
        (
            7198506.5328443246,
            0.0812854673,
            1.7174294882,
            3.4991320139,
            0.6678027378,
            422.5168026141,
        ),
        (
            7198053.9156890949,
            0.0798197372,
            1.7174165804,
            3.5089504189,
            1.3934319931,
            597.5289925841,
        ),
        (
            7185599.9070517309,
            0.0798078064,
            1.7175443634,
            3.4881187713,
            2.2414460706,
            422.5168026099,
        ),
        (
            7184901.5037797149,
            0.0801397057,
            1.7175554962,
            3.4933659390,
            3.1383575555,
            -0.0000000028,
        ),
        (
            7195505.0691928314,
            0.0790072358,
            1.7174331431,
            3.4820872250,
            4.0528587968,
            -422.5168026140,
        ),
        (
            7198053.5585980229,
            0.0803584749,
            1.7174229967,
            3.4734442897,
            4.8898755736,
            -597.5289925840,
        ),
        (
            7184591.1173784742,
            0.0801547555,
            1.7175588688,
            3.4949929506,
            5.6025688184,
            -422.5168026087,
        ),
        (
            7181478.5102035254,
            0.1496306659,
            1.7175806949,
            3.4879659757,
            0.0019215108,
            0.0000000140,
        ),
        (
            7200571.5026989058,
            0.1514567721,
            1.7174260880,
            3.4954699797,
            0.5839455326,
            798.7165162242,
        ),
        (
            7198483.9469006751,
            0.1498757265,
            1.7174115022,
            3.5011015424,
            1.2624022277,
            1129.5557297251,
        ),
        (
            7185794.5462470576,
            0.1497891434,
            1.7175398742,
            3.4898013120,
            2.1280406405,
            798.7165162118,
        ),
        (
            7185378.2525120089,
            0.1500913695,
            1.7175510083,
            3.4918752282,
            3.1398949178,
            -0.0000000076,
        ),
        (
            7194644.1613620706,
            0.1490904054,
            1.7174331680,
            3.4859319271,
            4.1607463286,
            -798.7165162234,
        ),
        (
            7198483.2383694621,
            0.1504233653,
            1.7174239326,
            3.4812668302,
            5.0209647203,
            -1129.5557297248,
        ),
        (
            7183802.8619856136,
            0.1501249115,
            1.7175679770,
            3.4934212469,
            5.6916574363,
            -798.7165162062,
        ),
    ]

    @testset "$prop" for prop in (:J2osc,)
        for (a, e, i, ω, M, ṙ_ref) in reference
            orb = KeplerianElements(
                date_to_jd(2023, 1, 1, 0, 0, 0),
                a,
                e,
                i,
                90.0 |> deg2rad,
                ω,
                mean_to_true_anomaly(e, M),
            )

            orbp = Propagators.init(Val(prop), orb)
            r, v = Propagators.propagate!(orbp, 0.0)

            ṙ = sum(r .* v) / sqrt(sum(abs2, r))

            @test ṙ ≈ ṙ_ref atol = 0.35
        end
    end
end
