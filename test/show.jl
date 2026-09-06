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

# Bodies of the rich representations after propagating the reference scenario by 100 s.
const SHOW_MEAN_ELEMENTS = (
    "  ├─ Mean Elements",
    "  │    Epoch             : 2.45995e6 (2023-01-01T00:00:00)",
    "  │    Semi-Major Axis   : 8000.0 km",
    "  │    Eccentricity      : 0.015",
    "  │    Inclination       : 28.5°",
    "  │    RA of Asc. Node   : 100.0°",
    "  │    Arg. of Periapsis : 200.0°",
    "  │    Mean Anomaly      : 43.79419642°",
)

const SHOW_PROPAGATION = (
    "  └─ Propagation",
    "       Last Instant : 100.0 s",
)

const SHOW_J2_BODY = join(
    (
        SHOW_MEAN_ELEMENTS...,
        "  ├─ Secular Rates",
        "  │    Mean Motion            : 12.14123799 rev/day",
        "  │    RAAN Rate              : -3.966769151 °/day",
        "  │    Arg. of Periapsis Rate : 6.458281745 °/day",
        "  ├─ Constants",
        "  │    R₀ : 6378.137 km",
        "  │    μm : 0.001239447462 rad/s",
        "  │    J₂ : 0.001082626174",
        SHOW_PROPAGATION...,
    ),
    '\n',
)

const SHOW_J4_BODY = join(
    (
        SHOW_MEAN_ELEMENTS...,
        "  ├─ Secular Rates",
        "  │    Mean Motion            : 12.14124728 rev/day",
        "  │    RAAN Rate              : -3.971642723 °/day",
        "  │    Arg. of Periapsis Rate : 6.466058362 °/day",
        "  ├─ Constants",
        "  │    R₀ : 6378.137 km",
        "  │    μm : 0.001239447462 rad/s",
        "  │    J₂ : 0.001082626174",
        "  │    J₄ : -1.6198976e-6",
        SHOW_PROPAGATION...,
    ),
    '\n',
)

const SHOW_TWOBODY_BODY = join(
    (
        SHOW_MEAN_ELEMENTS...,
        "  ├─ Secular Rates",
        "  │    Mean Motion : 12.13298837 rev/day",
        "  ├─ Constants",
        "  │    μ : 3.986004415e14 m³/s²",
        SHOW_PROPAGATION...,
    ),
    '\n',
)

const SHOW_SGP4_BODY = join(
    (
        "  ├─ Mean Elements",
        "  │    Epoch             : 2.45391e6 (2006-06-26T18:52:04.080)",
        "  │    Semi-Major Axis   : 7151.615424 km",
        "  │    Mean Motion       : 14.3547808 rev/day",
        "  │    Eccentricity      : 8.84e-5",
        "  │    Inclination       : 98.4283°",
        "  │    RA of Asc. Node   : 247.6961°",
        "  │    Arg. of Periapsis : 88.1964°",
        "  │    Mean Anomaly      : 271.9322°",
        "  │    B*                : 3.594e-5 1/ER",
        "  ├─ Constants",
        "  │    R₀  : 6378.137 km",
        "  │    XKE : 0.07436685317 er^(3/2)/min",
        "  │    J₂  : 0.001082629989",
        "  │    J₃  : -2.53215306e-6",
        "  │    J₄  : -1.61098761e-6",
        "  └─ Propagation",
        "       Last Instant : 1.666666667 min",
    ),
    '\n',
)

############################################################################################
#                                   Propagator Structures                                  #
############################################################################################

@testset "Propagator Structures" verbose = true begin
    scenarios = (
        (:J2,      "J2Propagator{Float64, Float64}",           SHOW_J2_BODY),
        (:J2osc,   "J2OsculatingPropagator{Float64, Float64}", SHOW_J2_BODY),
        (:J4,      "J4Propagator{Float64, Float64}",           SHOW_J4_BODY),
        (:J4osc,   "J4OsculatingPropagator{Float64, Float64}", SHOW_J4_BODY),
        (:TwoBody, "TwoBodyPropagator{Float64, Float64}",      SHOW_TWOBODY_BODY),
    )

    for (tag, name, body) in scenarios
        @testset "$name" begin
            orbp = Propagators.init(Val(tag), SHOW_ORB)
            Propagators.propagate!(orbp, 100)
            pd = Propagators.propagator_data(orbp)

            expected = "$name: Epoch = 2.45995e6 (2023-01-01T00:00:00)"
            @test sprint(show, pd) == expected

            expected = "$name:\n" * body
            @test sprint(show, MIME("text/plain"), pd) == expected

            # The body can be printed under another header.
            @test sprint(SatelliteToolboxBase.print_tree_body, pd) == body
        end
    end

    @testset "Float32" begin
        orbp = Propagators.init(Val(:J4), SHOW_ORB; j4c = J4C_EGM2008_F32)
        pd   = orbp.j4d

        expected = "J4Propagator{Float64, Float32}: Epoch = 2.45995e6 (2023-01-01T00:00:00)"
        @test sprint(show, pd) == expected

        result = sprint(show, MIME("text/plain"), pd)
        @test startswith(result, "J4Propagator{Float64, Float32}:\n")
        @test occursin("  │    Semi-Major Axis   : 8000.0 km\n", result)
        @test occursin("  │    R₀ : 6378.14 km\n", result)
        @test endswith(result, "  └─ Propagation\n       Last Instant : 0.0 s")
    end

    @testset "Uninitialized" begin
        for pd in (
            J2Propagator{Float64, Float64}(),
            J2OsculatingPropagator{Float64, Float64}(),
            J4Propagator{Float64, Float32}(),
            J4OsculatingPropagator{Float64, Float32}(),
            TwoBodyPropagator{Float64, Float64}(),
        )
            # The empty constructors mark the structure with a `NaN` propagation instant.
            @test isnan(pd.Δt)

            name = string(nameof(typeof(pd)), "{", join(typeof(pd).parameters, ", "), "}")
            @test sprint(show, pd) == "$name (not initialized)"
            expected = "$name:\n  Status : not initialized"
            @test sprint(show, MIME("text/plain"), pd) == expected
        end
    end

    @testset "Colors" begin
        orbp = Propagators.init(Val(:J2), SHOW_ORB)
        pd   = orbp.j2d

        plain   = sprint(show, MIME("text/plain"), pd)
        colored = sprint(show, MIME("text/plain"), pd; context = :color => true)

        @test !occursin("\e[", plain)
        @test startswith(colored, "\e[1mJ2Propagator{Float64, Float64}:\e[22m\n")
        @test occursin("\e[90mkm\e[39m", colored)
        @test replace(colored, r"\e\[[0-9;]*m" => "") == plain
    end
end

############################################################################################
#                                           API                                            #
############################################################################################

@testset "Orbit Propagators" verbose = true begin
    scenarios = (
        (
            :J2,
            "OrbitPropagatorJ2{Float64, Float64} (J2 Orbit Propagator)",
            SHOW_J2_BODY,
        ),
        (
            :J2osc,
            "OrbitPropagatorJ2Osculating{Float64, Float64} (J2 Osculating Orbit Propagator)",
            SHOW_J2_BODY,
        ),
        (
            :J4,
            "OrbitPropagatorJ4{Float64, Float64} (J4 Orbit Propagator)",
            SHOW_J4_BODY,
        ),
        (
            :J4osc,
            "OrbitPropagatorJ4Osculating{Float64, Float64} (J4 Osculating Orbit Propagator)",
            SHOW_J4_BODY,
        ),
        (
            :TwoBody,
            "OrbitPropagatorTwoBody{Float64, Float64} (Two-Body Orbit Propagator)",
            SHOW_TWOBODY_BODY,
        ),
    )

    for (tag, header, body) in scenarios
        @testset "$header" begin
            orbp = Propagators.init(Val(tag), SHOW_ORB)
            Propagators.propagate!(orbp, 100)

            @test Propagators.is_initialized(orbp)

            expected = "$header: Epoch = 2.45995e6 (2023-01-01T00:00:00)"
            @test sprint(show, orbp) == expected

            # The wrapper prints its header followed by the body of the structure.
            expected = "$header:\n" * body
            @test sprint(show, MIME("text/plain"), orbp) == expected
        end
    end

    @testset "OrbitPropagatorSgp4{Float64, Float64}" begin
        orbp = Propagators.init(Val(:SGP4), SHOW_TLE)
        Propagators.propagate!(orbp, 100)

        @test Propagators.is_initialized(orbp)
        @test Propagators.propagator_data(orbp) === orbp.sgp4d

        header   = "OrbitPropagatorSgp4{Float64, Float64} (SGP4 Orbit Propagator)"
        expected = "$header: Epoch = 2.45391e6 (2006-06-26T18:52:04.080)"
        @test sprint(show, orbp) == expected

        expected = "$header:\n" * SHOW_SGP4_BODY
        @test sprint(show, MIME("text/plain"), orbp) == expected
    end

    @testset "Uninitialized" begin
        scenarios = (
            (
                OrbitPropagatorJ2(J2Propagator{Float64, Float64}()),
                "OrbitPropagatorJ2{Float64, Float64} (J2 Orbit Propagator)",
            ),
            (
                OrbitPropagatorJ2Osculating(J2OsculatingPropagator{Float64, Float64}()),
                "OrbitPropagatorJ2Osculating{Float64, Float64} (J2 Osculating Orbit Propagator)",
            ),
            (
                OrbitPropagatorJ4(J4Propagator{Float64, Float64}()),
                "OrbitPropagatorJ4{Float64, Float64} (J4 Orbit Propagator)",
            ),
            (
                OrbitPropagatorJ4Osculating(J4OsculatingPropagator{Float64, Float64}()),
                "OrbitPropagatorJ4Osculating{Float64, Float64} (J4 Osculating Orbit Propagator)",
            ),
            (
                OrbitPropagatorTwoBody(TwoBodyPropagator{Float64, Float64}()),
                "OrbitPropagatorTwoBody{Float64, Float64} (Two-Body Orbit Propagator)",
            ),
            (
                OrbitPropagatorSgp4(Sgp4Propagator{Float64}(SGP4C_WGS84)),
                "OrbitPropagatorSgp4{Float64, Float64} (SGP4 Orbit Propagator)",
            ),
        )

        for (orbp, header) in scenarios
            @test !Propagators.is_initialized(orbp)

            # The propagators of this package mark the structure with a `NaN` instant.
            orbp isa OrbitPropagatorSgp4 || @test isnan(Propagators.last_instant(orbp))
            @test sprint(show, orbp) == "$header (not initialized)"

            expected = "$header:\n  Status : not initialized"
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
                "$header:",
                "  └─ Propagation",
                "       Epoch        : 2.45995e6 (2023-01-01T00:00:00)",
                "       Last Instant : 10.0 s",
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
