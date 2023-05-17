# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#
# Description
# ==========================================================================================
#
#   Tests of the SGP4 orbit propagator.
#
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

############################################################################################
#                                       Test Results
############################################################################################
#
# Scenario 01
# ==========================================================================================
#
# The SGP4 algorithm is highly tested in the package SatelliteToolboxSgp4.jl. Hence, we are
# only performing one of the test sets here to check the interface using the API.
#
# Using the TLE:
#
#    CBERS 2
#    1 28057U 03049A   06177.78615833  .00000060  00000-0  35940-4 0  1836
#    2 28057  98.4283 247.6961 0000884  88.1964 271.9322 14.35478080140550
#
# We must obtain the following data after propagation:
#
#      t [min]        X [km]         Y [km]         Z [km]      VX [km/s]    VY [km/s]    VZ [km/s]
#      0.00000000 -2715.28237486 -6619.26436889    -0.01341443 -1.008587273  0.422782003  7.385272942
#    120.00000000 -1816.87920942 -1835.78762132  6661.07926465  2.325140071  6.655669329  2.463394512
#    240.00000000  1483.17364291  5395.21248786  4448.65907172  2.560540387  4.039025766 -5.736648561
#    360.00000000  2801.25607157  5455.03931333 -3692.12865695 -0.595095864 -3.951923117 -6.298799125
#    480.00000000   411.09332812 -1728.99769152 -6935.45548810 -2.935970964 -6.684085058  1.492800886
#    600.00000000 -2506.52558454 -6628.98655094  -988.07784497 -1.390577189 -0.556164143  7.312736468
#    720.00000000 -2090.79884266 -2723.22832193  6266.13356576  1.992640665  6.337529519  3.411803080
#    840.00000000  1091.80560222  4809.88229503  5172.42897894  2.717483546  4.805518977 -5.030019896
#    960.00000000  2811.14062300  5950.65707171 -2813.23705389 -0.159662742 -3.121215491 -6.775341949
#   1080.00000000   805.72698304  -812.16627907 -7067.58483968 -2.798936020 -6.889265977  0.472770873
#   1200.00000000 -2249.59837532 -6505.84890714 -1956.72365062 -1.731234729 -1.528750230  7.096660885
#   1320.00000000 -2311.57375797 -3560.99112891  5748.16749600  1.626569751  5.890482233  4.293545048
#   1440.00000000   688.16056594  4124.87618964  5794.55994449  2.810973665  5.479585563 -4.224866316
#   1560.00000000  2759.94088230  6329.87271798 -1879.19518331  0.266930672 -2.222670878 -7.119390567
#   1680.00000000  1171.50677137   125.82053748 -7061.96626202 -2.605687852 -6.958489749 -0.556333225
#   1800.00000000 -1951.43708472 -6251.71945820 -2886.95472355 -2.024131483 -2.475214272  6.741537478
#   1920.00000000 -2475.70722288 -4331.90569958  5117.31234924  1.235823539  5.322743371  5.091281211
#   2040.00000000   281.46097847  3353.51057102  6302.87900650  2.840647273  6.047222485 -3.337085992
#   2160.00000000  2650.33118860  6584.33434851  -908.29027134  0.675457235 -1.274044972 -7.323921567
#   2280.00000000  1501.17226597  1066.31132756 -6918.71472952 -2.361891904 -6.889669974 -1.574718619
#   2400.00000000 -1619.73468334 -5871.14051991 -3760.56587071 -2.264093975 -3.376316601  6.254622256
#   2520.00000000 -2581.04202505 -5020.05572531  4385.92329047  0.829668458  4.645048038  5.789262667
#   2640.00000000  -119.22080628  2510.90620488  6687.45615459  2.807575712  6.496549689 -2.384136661
#   2760.00000000  2486.23806726  6708.18210028    80.43349581  1.057274905 -0.294294027 -7.384689123
#   2880.00000000  1788.42334580  1990.50530957 -6640.59337725 -2.074169091 -6.683381288 -2.562777776
#
############################################################################################

@testset "SGP4 orbit propagator" verbose = true begin
    expected_results = [
           0.00000000 -2715.28237486 -6619.26436889    -0.01341443 -1.008587273  0.422782003  7.385272942
         120.00000000 -1816.87920942 -1835.78762132  6661.07926465  2.325140071  6.655669329  2.463394512
         240.00000000  1483.17364291  5395.21248786  4448.65907172  2.560540387  4.039025766 -5.736648561
         360.00000000  2801.25607157  5455.03931333 -3692.12865695 -0.595095864 -3.951923117 -6.298799125
         480.00000000   411.09332812 -1728.99769152 -6935.45548810 -2.935970964 -6.684085058  1.492800886
         600.00000000 -2506.52558454 -6628.98655094  -988.07784497 -1.390577189 -0.556164143  7.312736468
         720.00000000 -2090.79884266 -2723.22832193  6266.13356576  1.992640665  6.337529519  3.411803080
         840.00000000  1091.80560222  4809.88229503  5172.42897894  2.717483546  4.805518977 -5.030019896
         960.00000000  2811.14062300  5950.65707171 -2813.23705389 -0.159662742 -3.121215491 -6.775341949
        1080.00000000   805.72698304  -812.16627907 -7067.58483968 -2.798936020 -6.889265977  0.472770873
        1200.00000000 -2249.59837532 -6505.84890714 -1956.72365062 -1.731234729 -1.528750230  7.096660885
        1320.00000000 -2311.57375797 -3560.99112891  5748.16749600  1.626569751  5.890482233  4.293545048
        1440.00000000   688.16056594  4124.87618964  5794.55994449  2.810973665  5.479585563 -4.224866316
        1560.00000000  2759.94088230  6329.87271798 -1879.19518331  0.266930672 -2.222670878 -7.119390567
        1680.00000000  1171.50677137   125.82053748 -7061.96626202 -2.605687852 -6.958489749 -0.556333225
        1800.00000000 -1951.43708472 -6251.71945820 -2886.95472355 -2.024131483 -2.475214272  6.741537478
        1920.00000000 -2475.70722288 -4331.90569958  5117.31234924  1.235823539  5.322743371  5.091281211
        2040.00000000   281.46097847  3353.51057102  6302.87900650  2.840647273  6.047222485 -3.337085992
        2160.00000000  2650.33118860  6584.33434851  -908.29027134  0.675457235 -1.274044972 -7.323921567
        2280.00000000  1501.17226597  1066.31132756 -6918.71472952 -2.361891904 -6.889669974 -1.574718619
        2400.00000000 -1619.73468334 -5871.14051991 -3760.56587071 -2.264093975 -3.376316601  6.254622256
        2520.00000000 -2581.04202505 -5020.05572531  4385.92329047  0.829668458  4.645048038  5.789262667
        2640.00000000  -119.22080628  2510.90620488  6687.45615459  2.807575712  6.496549689 -2.384136661
        2760.00000000  2486.23806726  6708.18210028    80.43349581  1.057274905 -0.294294027 -7.384689123
        2880.00000000  1788.42334580  1990.50530957 -6640.59337725 -2.074169091 -6.683381288 -2.562777776
    ]

    tle = tle"""
        CBERS 2
        1 28057U 03049A   06177.78615833  .00000060  00000-0  35940-4 0  1836
        2 28057  98.4283 247.6961 0000884  88.1964 271.9322 14.35478080140550
        """

    jd₀ = tle_epoch(tle)

    # General API Functions
    # ======================================================================================

    @testset "General API Functions" begin
        orbp = Propagators.init(Val(:SGP4), tle)
        @test Propagators.name(orbp) == "SGP4 Orbit Propagator"

        orbk = Propagators.mean_elements(orbp)
        @test orbk.t == Propagators.epoch(orbp)
        @test orbk.a == orbp.sgp4d.a_k * orbp.sgp4d.sgp4c.R0
        @test orbk.e == orbp.sgp4d.e_k
        @test orbk.i == orbp.sgp4d.i_k
        @test orbk.Ω == orbp.sgp4d.Ω_k
        @test orbk.ω == orbp.sgp4d.ω_k
        @test orbk.f == mean_to_true_anomaly(orbp.sgp4d.e_k, orbp.sgp4d.M_k)
    end

    # Float64
    # ======================================================================================

    @testset "Float64" begin
        T = Float64

        # Initialization Using TLE
        # ----------------------------------------------------------------------------------

        orbp = Propagators.init(Val(:SGP4), tle; sgp4c = sgp4c_wgs72)

        for k in size(expected_results)[1]
            r_teme, v_teme = Propagators.propagate!(orbp, 60 * expected_results[k, 1])

            @test Propagators.last_instant(orbp) == 60 * expected_results[k, 1]

            @test eltype(r_teme) == T
            @test eltype(v_teme) == T

            @test r_teme[1] ≈ 1000 * expected_results[k, 2] atol = 1e-5
            @test r_teme[2] ≈ 1000 * expected_results[k, 3] atol = 1e-5
            @test r_teme[3] ≈ 1000 * expected_results[k, 4] atol = 1e-5
            @test v_teme[1] ≈ 1000 * expected_results[k, 5] atol = 1e-6
            @test v_teme[2] ≈ 1000 * expected_results[k, 6] atol = 1e-6
            @test v_teme[3] ≈ 1000 * expected_results[k, 7] atol = 1e-6

            r_teme, v_teme = Propagators.propagate_to_epoch!(
                orbp,
                jd₀ + expected_results[k, 1] / 1440
            )

            @test Propagators.last_instant(orbp) == 60 * expected_results[k, 1]

            @test eltype(r_teme) == T
            @test eltype(v_teme) == T

            @test r_teme[1] ≈ 1000 * expected_results[k, 2] atol = 1e-5
            @test r_teme[2] ≈ 1000 * expected_results[k, 3] atol = 1e-5
            @test r_teme[3] ≈ 1000 * expected_results[k, 4] atol = 1e-5
            @test v_teme[1] ≈ 1000 * expected_results[k, 5] atol = 1e-6
            @test v_teme[2] ≈ 1000 * expected_results[k, 6] atol = 1e-6
            @test v_teme[3] ≈ 1000 * expected_results[k, 7] atol = 1e-6
        end

        # Test in-place initialization.
        orbp = OrbitPropagatorSgp4(Sgp4Propagator{Float64, T}())
        orbp.sgp4d.sgp4c = sgp4c_wgs72
        orbp.sgp4d.sgp4ds = SatelliteToolboxSgp4.Sgp4DeepSpace{T}()
        Propagators.init!(orbp, tle)

        for k in size(expected_results)[1]
            r_teme, v_teme = Propagators.propagate!(orbp, 60 * expected_results[k, 1])

            @test Propagators.last_instant(orbp) == 60 * expected_results[k, 1]

            @test r_teme[1] ≈ 1000 * expected_results[k, 2] atol = 1e-5
            @test r_teme[2] ≈ 1000 * expected_results[k, 3] atol = 1e-5
            @test r_teme[3] ≈ 1000 * expected_results[k, 4] atol = 1e-5
            @test v_teme[1] ≈ 1000 * expected_results[k, 5] atol = 1e-6
            @test v_teme[2] ≈ 1000 * expected_results[k, 6] atol = 1e-6
            @test v_teme[3] ≈ 1000 * expected_results[k, 7] atol = 1e-6
        end

        # Test simultaneous initialization and propagation.
        r_teme, v_teme, orbp = Propagators.propagate(
            Val(:SGP4),
            60 * expected_results[end, 1],
            tle;
            sgp4c = sgp4c_wgs72
        )

        @test Propagators.last_instant(orbp) == 60 * expected_results[end, 1]

        @test r_teme[1] ≈ 1000 * expected_results[end, 2] atol = 1e-5
        @test r_teme[2] ≈ 1000 * expected_results[end, 3] atol = 1e-5
        @test r_teme[3] ≈ 1000 * expected_results[end, 4] atol = 1e-5
        @test v_teme[1] ≈ 1000 * expected_results[end, 5] atol = 1e-6
        @test v_teme[2] ≈ 1000 * expected_results[end, 6] atol = 1e-6
        @test v_teme[3] ≈ 1000 * expected_results[end, 7] atol = 1e-6

        r_teme, v_teme, orbp = Propagators.propagate_to_epoch(
            Val(:SGP4),
            jd₀ + expected_results[end, 1] / 1440,
            tle;
            sgp4c = sgp4c_wgs72
        )

        @test Propagators.last_instant(orbp) == 60 * expected_results[end, 1]

        @test r_teme[1] ≈ 1000 * expected_results[end, 2] atol = 1e-5
        @test r_teme[2] ≈ 1000 * expected_results[end, 3] atol = 1e-5
        @test r_teme[3] ≈ 1000 * expected_results[end, 4] atol = 1e-5
        @test v_teme[1] ≈ 1000 * expected_results[end, 5] atol = 1e-6
        @test v_teme[2] ≈ 1000 * expected_results[end, 6] atol = 1e-6
        @test v_teme[3] ≈ 1000 * expected_results[end, 7] atol = 1e-6

        # Initialization Using Individual Elements
        # ----------------------------------------------------------------------------------

        orbp = Propagators.init(
            Val(:SGP4),
            tle_epoch(tle),
            tle.mean_motion * 2π / 86400,
            tle.eccentricity,
            tle.inclination         |> deg2rad,
            tle.raan                |> deg2rad,
            tle.argument_of_perigee |> deg2rad,
            tle.mean_anomaly        |> deg2rad,
            tle.bstar;
            sgp4c = sgp4c_wgs72
        )

        for k in size(expected_results)[1]
            r_teme, v_teme = Propagators.propagate!(orbp, 60 * expected_results[k, 1])

            @test Propagators.last_instant(orbp) == 60 * expected_results[k, 1]

            @test eltype(r_teme) == T
            @test eltype(v_teme) == T

            @test r_teme[1] ≈ 1000 * expected_results[k, 2] atol = 1e-5
            @test r_teme[2] ≈ 1000 * expected_results[k, 3] atol = 1e-5
            @test r_teme[3] ≈ 1000 * expected_results[k, 4] atol = 1e-5
            @test v_teme[1] ≈ 1000 * expected_results[k, 5] atol = 1e-6
            @test v_teme[2] ≈ 1000 * expected_results[k, 6] atol = 1e-6
            @test v_teme[3] ≈ 1000 * expected_results[k, 7] atol = 1e-6

            r_teme, v_teme = Propagators.propagate_to_epoch!(
                orbp,
                jd₀ + expected_results[k, 1] / 1440
            )

            @test Propagators.last_instant(orbp) == 60 * expected_results[k, 1]

            @test eltype(r_teme) == T
            @test eltype(v_teme) == T

            @test r_teme[1] ≈ 1000 * expected_results[k, 2] atol = 1e-5
            @test r_teme[2] ≈ 1000 * expected_results[k, 3] atol = 1e-5
            @test r_teme[3] ≈ 1000 * expected_results[k, 4] atol = 1e-5
            @test v_teme[1] ≈ 1000 * expected_results[k, 5] atol = 1e-6
            @test v_teme[2] ≈ 1000 * expected_results[k, 6] atol = 1e-6
            @test v_teme[3] ≈ 1000 * expected_results[k, 7] atol = 1e-6
        end

        # Test in-place initialization.
        orbp = OrbitPropagatorSgp4(Sgp4Propagator{Float64, T}())
        orbp.sgp4d.sgp4c = sgp4c_wgs72
        orbp.sgp4d.sgp4ds = SatelliteToolboxSgp4.Sgp4DeepSpace{T}()
        Propagators.init!(
            orbp,
            tle_epoch(tle),
            tle.mean_motion * 2π / 86400,
            tle.eccentricity,
            tle.inclination         |> deg2rad,
            tle.raan                |> deg2rad,
            tle.argument_of_perigee |> deg2rad,
            tle.mean_anomaly        |> deg2rad,
            tle.bstar
        )

        for k in size(expected_results)[1]
            r_teme, v_teme = Propagators.propagate!(orbp, 60 * expected_results[k, 1])

            @test Propagators.last_instant(orbp) == 60 * expected_results[k, 1]

            @test r_teme[1] ≈ 1000 * expected_results[k, 2] atol = 1e-5
            @test r_teme[2] ≈ 1000 * expected_results[k, 3] atol = 1e-5
            @test r_teme[3] ≈ 1000 * expected_results[k, 4] atol = 1e-5
            @test v_teme[1] ≈ 1000 * expected_results[k, 5] atol = 1e-6
            @test v_teme[2] ≈ 1000 * expected_results[k, 6] atol = 1e-6
            @test v_teme[3] ≈ 1000 * expected_results[k, 7] atol = 1e-6
        end

        # Test simultaneous initialization and propagation.
        r_teme, v_teme, orbp = Propagators.propagate(
            Val(:SGP4),
            60 * expected_results[end, 1],
            tle_epoch(tle),
            tle.mean_motion * 2π / 86400,
            tle.eccentricity,
            tle.inclination         |> deg2rad,
            tle.raan                |> deg2rad,
            tle.argument_of_perigee |> deg2rad,
            tle.mean_anomaly        |> deg2rad,
            tle.bstar;
            sgp4c = sgp4c_wgs72
        )

        @test Propagators.last_instant(orbp) == 60 * expected_results[end, 1]

        @test r_teme[1] ≈ 1000 * expected_results[end, 2] atol = 1e-5
        @test r_teme[2] ≈ 1000 * expected_results[end, 3] atol = 1e-5
        @test r_teme[3] ≈ 1000 * expected_results[end, 4] atol = 1e-5
        @test v_teme[1] ≈ 1000 * expected_results[end, 5] atol = 1e-6
        @test v_teme[2] ≈ 1000 * expected_results[end, 6] atol = 1e-6
        @test v_teme[3] ≈ 1000 * expected_results[end, 7] atol = 1e-6

        r_teme, v_teme, orbp = Propagators.propagate_to_epoch(
            Val(:SGP4),
            jd₀ + expected_results[end, 1] / 1440,
            tle_epoch(tle),
            tle.mean_motion * 2π / 86400,
            tle.eccentricity,
            tle.inclination         |> deg2rad,
            tle.raan                |> deg2rad,
            tle.argument_of_perigee |> deg2rad,
            tle.mean_anomaly        |> deg2rad,
            tle.bstar;
            sgp4c = sgp4c_wgs72
        )

        @test Propagators.last_instant(orbp) == 60 * expected_results[end, 1]

        @test r_teme[1] ≈ 1000 * expected_results[end, 2] atol = 1e-5
        @test r_teme[2] ≈ 1000 * expected_results[end, 3] atol = 1e-5
        @test r_teme[3] ≈ 1000 * expected_results[end, 4] atol = 1e-5
        @test v_teme[1] ≈ 1000 * expected_results[end, 5] atol = 1e-6
        @test v_teme[2] ≈ 1000 * expected_results[end, 6] atol = 1e-6
        @test v_teme[3] ≈ 1000 * expected_results[end, 7] atol = 1e-6
    end

    # Float32
    # ======================================================================================

    @testset "Float32" begin
        T = Float32

        # Initialization Using TLE
        # ----------------------------------------------------------------------------------

        orbp = Propagators.init(Val(:SGP4), tle; sgp4c = sgp4c_wgs72_f32)

        for k in size(expected_results)[1]
            r_teme, v_teme = Propagators.propagate!(orbp, 60 * expected_results[k, 1])

            @test Propagators.last_instant(orbp) == 60 * expected_results[k, 1]

            @test eltype(r_teme) == T
            @test eltype(v_teme) == T

            @test r_teme[1] ≈ 1000 * expected_results[k, 2] atol = 60
            @test r_teme[2] ≈ 1000 * expected_results[k, 3] atol = 60
            @test r_teme[3] ≈ 1000 * expected_results[k, 4] atol = 60
            @test v_teme[1] ≈ 1000 * expected_results[k, 5] atol = 8e-1
            @test v_teme[2] ≈ 1000 * expected_results[k, 6] atol = 8e-1
            @test v_teme[3] ≈ 1000 * expected_results[k, 7] atol = 8e-1

            r_teme, v_teme = Propagators.propagate_to_epoch!(
                orbp,
                jd₀ + expected_results[k, 1] / 1440
            )

            @test Propagators.last_instant(orbp) == 60 * expected_results[k, 1]

            @test eltype(r_teme) == T
            @test eltype(v_teme) == T

            @test r_teme[1] ≈ 1000 * expected_results[k, 2] atol = 60
            @test r_teme[2] ≈ 1000 * expected_results[k, 3] atol = 60
            @test r_teme[3] ≈ 1000 * expected_results[k, 4] atol = 60
            @test v_teme[1] ≈ 1000 * expected_results[k, 5] atol = 8e-1
            @test v_teme[2] ≈ 1000 * expected_results[k, 6] atol = 8e-1
            @test v_teme[3] ≈ 1000 * expected_results[k, 7] atol = 8e-1
        end

        # Test in-place initialization.
        orbp = OrbitPropagatorSgp4(Sgp4Propagator{Float64, T}())
        orbp.sgp4d.sgp4c = sgp4c_wgs72_f32
        orbp.sgp4d.sgp4ds = SatelliteToolboxSgp4.Sgp4DeepSpace{T}()
        Propagators.init!(orbp, tle)

        for k in size(expected_results)[1]
            r_teme, v_teme = Propagators.propagate!(orbp, 60 * expected_results[k, 1])

            @test Propagators.last_instant(orbp) == 60 * expected_results[k, 1]

            @test r_teme[1] ≈ 1000 * expected_results[k, 2] atol = 60
            @test r_teme[2] ≈ 1000 * expected_results[k, 3] atol = 60
            @test r_teme[3] ≈ 1000 * expected_results[k, 4] atol = 60
            @test v_teme[1] ≈ 1000 * expected_results[k, 5] atol = 8e-1
            @test v_teme[2] ≈ 1000 * expected_results[k, 6] atol = 8e-1
            @test v_teme[3] ≈ 1000 * expected_results[k, 7] atol = 8e-1
        end

        # Test simultaneous initialization and propagation.
        r_teme, v_teme, orbp = Propagators.propagate(
            Val(:SGP4),
            60 * expected_results[end, 1],
            tle;
            sgp4c = sgp4c_wgs72_f32
        )

        @test Propagators.last_instant(orbp) == 60 * expected_results[end, 1]

        @test r_teme[1] ≈ 1000 * expected_results[end, 2] atol = 60
        @test r_teme[2] ≈ 1000 * expected_results[end, 3] atol = 60
        @test r_teme[3] ≈ 1000 * expected_results[end, 4] atol = 60
        @test v_teme[1] ≈ 1000 * expected_results[end, 5] atol = 8e-1
        @test v_teme[2] ≈ 1000 * expected_results[end, 6] atol = 8e-1
        @test v_teme[3] ≈ 1000 * expected_results[end, 7] atol = 8e-1

        r_teme, v_teme, orbp = Propagators.propagate_to_epoch(
            Val(:SGP4),
            jd₀ + expected_results[end, 1] / 1440,
            tle;
            sgp4c = sgp4c_wgs72_f32
        )

        @test Propagators.last_instant(orbp) == 60 * expected_results[end, 1]

        @test r_teme[1] ≈ 1000 * expected_results[end, 2] atol = 60
        @test r_teme[2] ≈ 1000 * expected_results[end, 3] atol = 60
        @test r_teme[3] ≈ 1000 * expected_results[end, 4] atol = 60
        @test v_teme[1] ≈ 1000 * expected_results[end, 5] atol = 8e-1
        @test v_teme[2] ≈ 1000 * expected_results[end, 6] atol = 8e-1
        @test v_teme[3] ≈ 1000 * expected_results[end, 7] atol = 8e-1

        # Initialization Using Individual Elements
        # ----------------------------------------------------------------------------------

        orbp = Propagators.init(
            Val(:SGP4),
            tle_epoch(tle),
            tle.mean_motion * 2π / 86400,
            tle.eccentricity,
            tle.inclination         |> deg2rad,
            tle.raan                |> deg2rad,
            tle.argument_of_perigee |> deg2rad,
            tle.mean_anomaly        |> deg2rad,
            tle.bstar;
            sgp4c = sgp4c_wgs72_f32
        )

        for k in size(expected_results)[1]
            r_teme, v_teme = Propagators.propagate!(orbp, 60 * expected_results[k, 1])

            @test Propagators.last_instant(orbp) == 60 * expected_results[k, 1]

            @test eltype(r_teme) == T
            @test eltype(v_teme) == T

            @test r_teme[1] ≈ 1000 * expected_results[k, 2] atol = 60
            @test r_teme[2] ≈ 1000 * expected_results[k, 3] atol = 60
            @test r_teme[3] ≈ 1000 * expected_results[k, 4] atol = 60
            @test v_teme[1] ≈ 1000 * expected_results[k, 5] atol = 8e-1
            @test v_teme[2] ≈ 1000 * expected_results[k, 6] atol = 8e-1
            @test v_teme[3] ≈ 1000 * expected_results[k, 7] atol = 8e-1

            r_teme, v_teme = Propagators.propagate_to_epoch!(
                orbp,
                jd₀ + expected_results[k, 1] / 1440
            )

            @test Propagators.last_instant(orbp) == 60 * expected_results[k, 1]

            @test eltype(r_teme) == T
            @test eltype(v_teme) == T

            @test r_teme[1] ≈ 1000 * expected_results[k, 2] atol = 60
            @test r_teme[2] ≈ 1000 * expected_results[k, 3] atol = 60
            @test r_teme[3] ≈ 1000 * expected_results[k, 4] atol = 60
            @test v_teme[1] ≈ 1000 * expected_results[k, 5] atol = 8e-1
            @test v_teme[2] ≈ 1000 * expected_results[k, 6] atol = 8e-1
            @test v_teme[3] ≈ 1000 * expected_results[k, 7] atol = 8e-1
        end

        # Test in-place initialization.
        orbp = OrbitPropagatorSgp4(Sgp4Propagator{Float64, T}())
        orbp.sgp4d.sgp4c = sgp4c_wgs72_f32
        orbp.sgp4d.sgp4ds = SatelliteToolboxSgp4.Sgp4DeepSpace{T}()
        Propagators.init!(
            orbp,
            tle_epoch(tle),
            tle.mean_motion * 2π / 86400,
            tle.eccentricity,
            tle.inclination         |> deg2rad,
            tle.raan                |> deg2rad,
            tle.argument_of_perigee |> deg2rad,
            tle.mean_anomaly        |> deg2rad,
            tle.bstar
        )

        for k in size(expected_results)[1]
            r_teme, v_teme = Propagators.propagate!(orbp, 60 * expected_results[k, 1])

            @test Propagators.last_instant(orbp) == 60 * expected_results[k, 1]

            @test r_teme[1] ≈ 1000 * expected_results[k, 2] atol = 60
            @test r_teme[2] ≈ 1000 * expected_results[k, 3] atol = 60
            @test r_teme[3] ≈ 1000 * expected_results[k, 4] atol = 60
            @test v_teme[1] ≈ 1000 * expected_results[k, 5] atol = 8e-1
            @test v_teme[2] ≈ 1000 * expected_results[k, 6] atol = 8e-1
            @test v_teme[3] ≈ 1000 * expected_results[k, 7] atol = 8e-1
        end

        # Test simultaneous initialization and propagation.
        r_teme, v_teme, orbp = Propagators.propagate(
            Val(:SGP4),
            60 * expected_results[end, 1],
            tle_epoch(tle),
            tle.mean_motion * 2π / 86400,
            tle.eccentricity,
            tle.inclination         |> deg2rad,
            tle.raan                |> deg2rad,
            tle.argument_of_perigee |> deg2rad,
            tle.mean_anomaly        |> deg2rad,
            tle.bstar;
            sgp4c = sgp4c_wgs72_f32
        )

        @test Propagators.last_instant(orbp) == 60 * expected_results[end, 1]

        @test r_teme[1] ≈ 1000 * expected_results[end, 2] atol = 60
        @test r_teme[2] ≈ 1000 * expected_results[end, 3] atol = 60
        @test r_teme[3] ≈ 1000 * expected_results[end, 4] atol = 60
        @test v_teme[1] ≈ 1000 * expected_results[end, 5] atol = 8e-1
        @test v_teme[2] ≈ 1000 * expected_results[end, 6] atol = 8e-1
        @test v_teme[3] ≈ 1000 * expected_results[end, 7] atol = 8e-1

        r_teme, v_teme, orbp = Propagators.propagate_to_epoch(
            Val(:SGP4),
            jd₀ + expected_results[end, 1] / 1440,
            tle_epoch(tle),
            tle.mean_motion * 2π / 86400,
            tle.eccentricity,
            tle.inclination         |> deg2rad,
            tle.raan                |> deg2rad,
            tle.argument_of_perigee |> deg2rad,
            tle.mean_anomaly        |> deg2rad,
            tle.bstar;
            sgp4c = sgp4c_wgs72_f32
        )

        @test Propagators.last_instant(orbp) == 60 * expected_results[end, 1]

        @test r_teme[1] ≈ 1000 * expected_results[end, 2] atol = 60
        @test r_teme[2] ≈ 1000 * expected_results[end, 3] atol = 60
        @test r_teme[3] ≈ 1000 * expected_results[end, 4] atol = 60
        @test v_teme[1] ≈ 1000 * expected_results[end, 5] atol = 8e-1
        @test v_teme[2] ≈ 1000 * expected_results[end, 6] atol = 8e-1
        @test v_teme[3] ≈ 1000 * expected_results[end, 7] atol = 8e-1
    end
end
