## Description #############################################################################
#
# Tests related to the printed representations of the propagators.
#
############################################################################################

# Propagator without a data structure to test the default printing of the API.
struct DummyShowPropagator{Tepoch, T} <: OrbitPropagator{Tepoch, T}
    epoch::Tepoch
    Δt::T
end

Propagators.epoch(orbp::DummyShowPropagator)        = orbp.epoch
Propagators.last_instant(orbp::DummyShowPropagator) = orbp.Δt

# == Reference Scenario ====================================================================

const SHOW_JD₀ = date_to_jd(2023, 1, 1, 0, 0, 0)

# The anomaly is the true anomaly (45°), which is converted to the mean anomaly
# (43.79419642°) by the propagators.
const SHOW_ORB = KeplerianElements(
    SHOW_JD₀,
    8000e3,
    0.015,
    28.5 |> deg2rad,
    100  |> deg2rad,
    200  |> deg2rad,
    45   |> deg2rad,
)

const SHOW_TLE = tle"""
    CBERS 2
    1 28057U 03049A   06177.78615833  .00000060  00000-0  35940-4 0  1836
    2 28057  98.4283 247.6961 0000884  88.1964 271.9322 14.35478080140550
    """

const SHOW_J2_ROWS = join(
    (
        "                   R₀ : 6378.14       km",
        "                   μm :    0.00123945 rad / s",
        "                   J₂ :    0.00108263",
        "      Semi-major axis : 8000.00000000 km",
        "         Eccentricity :    0.01500000",
        "          Inclination :   28.50000000 °",
        "                 RAAN :  100.00000000 °",
        "      Arg. of perigee :  200.00000000 °",
        "         Mean anomaly :   43.79419642 °",
        "          Mean motion :   12.14123799 rev / day",
        "            RAAN rate :   -3.96676915 ° / day",
        " Arg. of perigee rate :    6.45828175 ° / day",
        "     Last propagation :  100.0        s",
    ),
    '\n',
)

const SHOW_J4_ROWS = join(
    (
        "                   R₀ : 6378.14       km",
        "                   μm :    0.00123945 rad / s",
        "                   J₂ :    0.00108263",
        "                   J₄ :   -1.6199e-6",
        "      Semi-major axis : 8000.00000000 km",
        "         Eccentricity :    0.01500000",
        "          Inclination :   28.50000000 °",
        "                 RAAN :  100.00000000 °",
        "      Arg. of perigee :  200.00000000 °",
        "         Mean anomaly :   43.79419642 °",
        "          Mean motion :   12.14124728 rev / day",
        "            RAAN rate :   -3.97164272 ° / day",
        " Arg. of perigee rate :    6.46605836 ° / day",
        "     Last propagation :  100.0        s",
    ),
    '\n',
)

const SHOW_TWOBODY_ROWS = join(
    (
        "                μ :    3.986e14   m³ / s²",
        "  Semi-major axis : 8000.00000000 km",
        "     Eccentricity :    0.01500000",
        "      Inclination :   28.50000000 °",
        "             RAAN :  100.00000000 °",
        "  Arg. of perigee :  200.00000000 °",
        "     Mean anomaly :   43.79419642 °",
        "      Mean motion :   12.13298837 rev / day",
        " Last propagation :  100.0        s",
    ),
    '\n',
)

const SHOW_SGP4_ROWS = join(
    (
        "            Epoch :   2.45391e6 (2006-06-26T18:52:04.080)",
        "      Mean motion :  14.35478080 rev / day",
        "     Eccentricity :   0.00008840",
        "      Inclination :  98.42830000 °",
        "             RAAN : 247.69610000 °",
        "  Arg. of perigee :  88.19640000 °",
        "     Mean anomaly : 271.93220000 °",
        "               B* :   3.594e-5   1 / er",
        " Last propagation :   1.66667    min",
    ),
    '\n',
)

# The Keplerian propagators print the epoch with the same width of the other rows.
const SHOW_EPOCH_ROW = "Epoch :    2.45995e6 (2023-01-01T00:00:00)"

"""
    indent(str::String) -> String

Indent every line of `str` by two spaces, as the API wrappers print the wrapped structure.
"""
indent(str::String) = join("  " .* split(str, '\n'), '\n')

############################################################################################
#                                   Propagator Structures                                  #
############################################################################################

@testset "Propagator Structures" verbose = true begin
    scenarios = (
        (:J2,      "J2Propagator{Float64, Float64}",           SHOW_J2_ROWS),
        (:J2osc,   "J2OsculatingPropagator{Float64, Float64}", SHOW_J2_ROWS),
        (:J4,      "J4Propagator{Float64, Float64}",           SHOW_J4_ROWS),
        (:J4osc,   "J4OsculatingPropagator{Float64, Float64}", SHOW_J4_ROWS),
        (:TwoBody, "TwoBodyPropagator{Float64, Float64}",      SHOW_TWOBODY_ROWS),
    )

    for (tag, name, rows) in scenarios
        @testset "$name" begin
            orbp = Propagators.init(Val(tag), SHOW_ORB)
            Propagators.propagate!(orbp, 100)
            pd = Propagators.propagator_data(orbp)

            expected = "$name: Epoch = 2.45995e6 (2023-01-01T00:00:00)"
            @test sprint(show, pd) == expected

            # The epoch label is right-aligned with the widest label of the rows.
            label_width = maximum(length, strip.(first.(split.(split(rows, '\n'), " : "))))
            epoch_row   = lpad("", label_width - length("Epoch") + 1) * SHOW_EPOCH_ROW

            expected = "$name:\n" * epoch_row * "\n" * rows
            @test sprint(show, MIME("text/plain"), pd) == expected
        end
    end

    @testset "Float32" begin
        orbp = Propagators.init(Val(:J4), SHOW_ORB; j4c = J4C_EGM2008_F32)
        pd   = orbp.j4d

        expected = "J4Propagator{Float64, Float32}: Epoch = 2.45995e6 (2023-01-01T00:00:00)"
        @test sprint(show, pd) == expected

        result = sprint(show, MIME("text/plain"), pd)
        @test startswith(result, "J4Propagator{Float64, Float32}:\n")
        @test occursin("Semi-major axis : 8000.00000000 km", result)
        @test occursin("Last propagation :    0.0        s", result)
    end

    @testset "Uninitialized" begin
        for pd in (
            J2OsculatingPropagator{Float64, Float64}(),
            J4OsculatingPropagator{Float64, Float32}(),
        )
            name = string(nameof(typeof(pd)), "{", join(typeof(pd).parameters, ", "), "}")
            @test sprint(show, pd) == "$name (not initialized)"
            expected = "$name:\n Status : not initialized"
            @test sprint(show, MIME("text/plain"), pd) == expected
        end
    end

    @testset "Colors" begin
        orbp = Propagators.init(Val(:J2), SHOW_ORB)
        pd   = orbp.j2d

        plain   = sprint(show, MIME("text/plain"), pd)
        colored = sprint(show, MIME("text/plain"), pd; context = :color => true)

        @test !occursin("\e[", plain)
        @test occursin("\e[1m", colored)
        @test occursin("\e[1m                Epoch : \e[22m", colored)
        @test replace(colored, r"\e\[[0-9;]*m" => "") == plain
    end
end

############################################################################################
#                                           API                                            #
############################################################################################

@testset "Orbit Propagators" verbose = true begin
    scenarios = (
        (:J2,      "OrbitPropagatorJ2{Float64, Float64} (J2 Orbit Propagator)"),
        (:J2osc,   "OrbitPropagatorJ2Osculating{Float64, Float64} (J2 Osculating Orbit Propagator)"),
        (:J4,      "OrbitPropagatorJ4{Float64, Float64} (J4 Orbit Propagator)"),
        (:J4osc,   "OrbitPropagatorJ4Osculating{Float64, Float64} (J4 Osculating Orbit Propagator)"),
        (:TwoBody, "OrbitPropagatorTwoBody{Float64, Float64} (Two-Body Orbit Propagator)"),
    )

    for (tag, header) in scenarios
        @testset "$header" begin
            orbp = Propagators.init(Val(tag), SHOW_ORB)
            Propagators.propagate!(orbp, 100)
            pd = Propagators.propagator_data(orbp)

            @test Propagators.is_initialized(orbp)

            expected = "$header: Epoch = 2.45995e6 (2023-01-01T00:00:00)"
            @test sprint(show, orbp) == expected

            # The wrapper prints its header followed by the rich form of the structure.
            expected = "$header:\n" * indent(sprint(show, MIME("text/plain"), pd))
            @test sprint(show, MIME("text/plain"), orbp) == expected
        end
    end

    @testset "OrbitPropagatorJ2{Float64, Float64} (Exact Output)" begin
        orbp = Propagators.init(Val(:J2), SHOW_ORB)
        Propagators.propagate!(orbp, 100)

        expected = join(
            (
                "OrbitPropagatorJ2{Float64, Float64} (J2 Orbit Propagator):",
                "  J2Propagator{Float64, Float64}:",
                "                  Epoch :    2.45995e6 (2023-01-01T00:00:00)",
                indent(SHOW_J2_ROWS),
            ),
            '\n',
        )
        @test sprint(show, MIME("text/plain"), orbp) == expected
    end

    @testset "OrbitPropagatorSgp4{Float64, Float64}" begin
        orbp = Propagators.init(Val(:SGP4), SHOW_TLE)
        Propagators.propagate!(orbp, 100)

        @test Propagators.is_initialized(orbp)
        @test Propagators.propagator_data(orbp) === orbp.sgp4d

        header   = "OrbitPropagatorSgp4{Float64, Float64} (SGP4 Orbit Propagator)"
        expected = "$header: Epoch = 2.45391e6 (2006-06-26T18:52:04.080)"
        @test sprint(show, orbp) == expected

        expected = join(
            (
                "$header:",
                "  Sgp4Propagator{Float64, Float64} (SGP4):",
                indent(SHOW_SGP4_ROWS),
            ),
            '\n',
        )
        @test sprint(show, MIME("text/plain"), orbp) == expected
    end

    @testset "Uninitialized" begin
        scenarios = (
            (
                OrbitPropagatorJ2Osculating(J2OsculatingPropagator{Float64, Float64}()),
                "OrbitPropagatorJ2Osculating{Float64, Float64} (J2 Osculating Orbit Propagator)",
                "J2OsculatingPropagator{Float64, Float64}",
            ),
            (
                OrbitPropagatorJ4Osculating(J4OsculatingPropagator{Float64, Float64}()),
                "OrbitPropagatorJ4Osculating{Float64, Float64} (J4 Osculating Orbit Propagator)",
                "J4OsculatingPropagator{Float64, Float64}",
            ),
            (
                OrbitPropagatorSgp4(Sgp4Propagator{Float64}(SGP4C_WGS84)),
                "OrbitPropagatorSgp4{Float64, Float64} (SGP4 Orbit Propagator)",
                "Sgp4Propagator{Float64, Float64}",
            ),
        )

        for (orbp, header, name) in scenarios
            @test !Propagators.is_initialized(orbp)
            @test sprint(show, orbp) == "$header (not initialized)"

            expected = "$header:\n  $name:\n   Status : not initialized"
            @test sprint(show, MIME("text/plain"), orbp) == expected
        end
    end

    @testset "Default Printing" begin
        orbp = DummyShowPropagator(SHOW_JD₀, 10.0)

        @test Propagators.is_initialized(orbp)
        @test Propagators.propagator_data(orbp) === nothing

        # The default name is the type, which is not repeated in the header.
        header   = "DummyShowPropagator{Float64, Float64}"
        expected = "$header: Epoch = 2.45995e6 (2023-01-01T00:00:00)"
        @test sprint(show, orbp) == expected

        expected = join(
            (
                "DummyShowPropagator{Float64, Float64}:",
                "            Epoch :  2.45995e6 (2023-01-01T00:00:00)",
                " Last propagation : 10.0 s",
            ),
            '\n',
        )
        @test sprint(show, MIME("text/plain"), orbp) == expected
    end

    @testset "Colors" begin
        orbp = Propagators.init(Val(:J2), SHOW_ORB)

        plain   = sprint(show, MIME("text/plain"), orbp)
        colored = sprint(show, MIME("text/plain"), orbp; context = :color => true)

        @test !occursin("\e[", plain)
        @test occursin("\e[1m", colored)
        @test replace(colored, r"\e\[[0-9;]*m" => "") == plain
    end
end
