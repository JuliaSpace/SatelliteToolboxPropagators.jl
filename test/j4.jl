# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#
# Description
# ==========================================================================================
#
#   Tests of the J4 orbit propagator.
#
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#
# References
# ==========================================================================================
#
#   [1] Vallado, D. A (2013). Fundamentals of Astrodynamics and Applications.
#       Microcosm Press, Hawthorn, CA, USA.
#
#   [2] Hoots, F. R., Roehrich, R. L (1980). Models for Propagation of NORAD
#       Elements Set. Spacetrack Report No. 3.
#
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

############################################################################################
#                                       Test Results
############################################################################################
#
# Scenario 01
# ==========================================================================================
#
# Source: https://github.com/JuliaSpace/SatelliteToolbox.jl/issues/91
#
# The STK provided the following result for the J4Perturbation orbit propagator.
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
# Classical Keplerian elements on 2023-01-05T00:00:00.000
#
#   a = 7000 km
#   e = 0.015
#   i = 28.5°
#   Ω = 84.158846°
#   ω = 225.864212°
#   M = 247.190278°
#   f = 245.617459°
#
# NOTE: We cannot test RAAN value.
#
# We found a problem when implementing the analytical J4 orbit propagation when comparing
# with the STK's results. We needed to flip the sign of the J4 perturbation term to match
# the results provided by STK. However, this operation does not seem right given the
# equations in [1] and the secular perturbation term of SGP4.
#
############################################################################################

@testset "J4 orbit propagator" begin
    jd₀ = date_to_jd(2023, 1, 1, 0, 0, 0)
    jd₁ = date_to_jd(2023, 1, 5, 0, 0, 0)

    # General API Functions
    # ======================================================================================

    let
        orb = KeplerianElements(0.0, 8000.0e3, 0.0, 0.0, 0.0, 0.0, 0.0)
        orbp = Propagators.init(Val(:J4), orb)
        @test Propagators.name(orbp) == "J4 Orbit Propagator"
    end

    # Float64
    # ======================================================================================

    let T = Float64
        orb = KeplerianElements(
            jd₀,
            T(8000e3),
            T(0.015),
            T(28.5) |>deg2rad,
            T(100) |>deg2rad,
            T(200) |>deg2rad,
            T(45) |> deg2rad
        )

        orbp = Propagators.init(Val(:J4), orb; j4c = j4c_egm08)
        r, v = Propagators.propagate_to_epoch!(orbp, jd₁)

        @test eltype(r) == T
        @test eltype(v) == T

        orbk = Propagators.mean_elements(orbp)
        M_k  = true_to_mean_anomaly(orbk.e, orbk.f)

        @test orbk.a            ≈ 8000.0e3    (atol = 1e-3)
        @test orbk.e            ≈    0.015    (atol = 1e-6)
        @test orbk.i |> rad2deg ≈   28.5      (atol = 1e-6)
        @test orbk.ω |> rad2deg ≈  225.864212 (atol = 4e-3)
        @test orbk.f |> rad2deg ≈  245.617459 (atol = 4e-3)

        @test_broken orbk.Ω |> rad2deg ≈ 84.158846 (atol = 4e-3)
    end

    let T = Float32
        orb = KeplerianElements(
            jd₀,
            T(8000e3),
            T(0.015),
            T(28.5) |>deg2rad,
            T(100) |>deg2rad,
            T(200) |>deg2rad,
            T(45) |> deg2rad
        )

        orbp = Propagators.init(Val(:J4), orb; j4c = j4c_egm08_f32)
        r, v = Propagators.propagate_to_epoch!(orbp, jd₁)

        @test eltype(r) == T
        @test eltype(v) == T

        orbk = Propagators.mean_elements(orbp)
        M_k  = true_to_mean_anomaly(orbk.e, orbk.f)

        @test orbk.a            ≈ 8000.0e3    (atol = 1e-1)
        @test orbk.e            ≈    0.015    (atol = 1e-6)
        @test orbk.i |> rad2deg ≈   28.5      (atol = 2e-6)
        @test orbk.ω |> rad2deg ≈  225.864212 (atol = 4e-3)
        @test orbk.f |> rad2deg ≈  245.617459 (atol = 4e-3)

        @test_broken orbk.Ω |> rad2deg ≈ 84.158846 (atol = 4e-3)
    end
end
