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
    # propagator the long-term periodics are not taken into account.
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
        j2d    = J2Propagator{Float64, Float64}(orb, orb, j2c_egm2008, 0, 0, 0, 0, 0)
        j2oscd = J2OsculatingPropagator{Float64, Float64}(j2d, 0, orb)

        # Test some random fields.
        @test j2oscd.j2d == j2d
        @test j2oscd.Δt  == 0
        @test j2d.orbk   == orb
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
            T(90)     |> deg2rad,
            T(200)    |> deg2rad,
            T(45)     |> deg2rad
        )

        # Test all the results.
        orbp = Propagators.init(Val(:J2osc), orb; j2c = j2c_egm2008)

        for k in size(results, 1)
            r, v = Propagators.propagate!(orbp, results[k, 1])

            @test results[k, 2] ≈ r[1] / 1000 atol = 2e-1
            @test results[k, 3] ≈ r[2] / 1000 atol = 2e-1
            @test results[k, 4] ≈ r[3] / 1000 atol = 2e-1
            @test results[k, 5] ≈ v[1] / 1000 atol = 1e-3
            @test results[k, 6] ≈ v[2] / 1000 atol = 1e-3
            @test results[k, 7] ≈ v[3] / 1000 atol = 1e-3
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
            Val(:J2osc),
            jd₀ + results[end, 1] / 86400,
            orb
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
            T(90)     |> deg2rad,
            T(200)    |> deg2rad,
            T(45)     |> deg2rad
        )

        # Test all the results.
        orbp = Propagators.init(Val(:J2osc), orb; j2c = j2c_egm2008_f32)

        for k in size(results, 1)
            r, v = Propagators.propagate!(orbp, results[k, 1])

            @test results[k, 2] ≈ r[1] / 1000 atol = 2e-1
            @test results[k, 3] ≈ r[2] / 1000 atol = 2e-1
            @test results[k, 4] ≈ r[3] / 1000 atol = 2e-1
            @test results[k, 5] ≈ v[1] / 1000 atol = 1e-3
            @test results[k, 6] ≈ v[2] / 1000 atol = 1e-3
            @test results[k, 7] ≈ v[3] / 1000 atol = 1e-3
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
            Val(:J2osc),
            results[end, 1],
            orb;
            j2c = j2c_egm2008_f32
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
            Val(:J2osc),
            jd₀ + results[end, 1] / 86400,
            orb;
            j2c = j2c_egm2008_f32
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
        90     |> deg2rad,
        200    |> deg2rad,
        45     |> deg2rad
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
                Val(:J2osc),
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
                Val(:J2osc),
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
                Val(:J2osc),
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

        @test_throws ArgumentError Propagators.fit_mean_elements(Val(:J2osc), vjd[1:end-1], vr_i, vv_i)
        @test_throws ArgumentError Propagators.fit_mean_elements(Val(:J2osc), vjd, vr_i[1:end-1], vv_i)
        @test_throws ArgumentError Propagators.fit_mean_elements(Val(:J2osc), vjd, vr_i, vv_i[1:end-1])

        # == Wrong dimensions in the weight vector =========================================

        @test_throws ArgumentError Propagators.fit_mean_elements(
            Val(:J2osc),
            vjd,
            vr_i,
            vv_i;
            weight_vector = [1, 2, 3, 4, 5]
        )
    end
end

@testset "Update J2 Osculating Mean Elements Epoch" begin
    orb_input = KeplerianElements(
        DateTime("2023-01-01") |> datetime2julian,
        7130.982e3,
        0.001111,
        98.405 |> deg2rad,
        90     |> deg2rad,
        200    |> deg2rad,
        45     |> deg2rad
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
                T(100)  |> deg2rad,
                T(200)  |> deg2rad,
                T(25)   |> deg2rad
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
