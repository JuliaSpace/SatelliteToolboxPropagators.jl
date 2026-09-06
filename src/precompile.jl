## Description #############################################################################
#
#   Precompilation.
#
############################################################################################

PrecompileTools.@setup_workload begin
    orb = KeplerianElements(
        DateTime("2023-01-01") |> datetime2julian,
        Float64(8000e3),
        Float64(0.015),
        Float64(28.5) |> deg2rad,
        Float64(100) |> deg2rad,
        Float64(200) |> deg2rad,
        Float64(45) |> deg2rad,
    )

    tle = tle"""
          AMAZONIA 1
          1 47699U 21015A   23083.68657856  .00000000  00000-8  43000-3 0  9999
          2 47699  98.4304 162.1097 0001247 136.2017 223.9283 14.40814394108652
          """

    vr_i = [
        @SVector([-6792.402703741442, 2192.6458461287293, 0.18851758695295118]) .* 1000,
        @SVector([-6357.88873265975, 2391.9476768911686, 2181.838771262736]) .* 1000,
    ]

    vv_i = [
        @SVector([0.3445760107690598, 1.0395135806993514, 7.393686131436984]) .* 1000,
        @SVector([2.5285015912807003, 0.27812476784300005, 7.030323100703928]) .* 1000,
    ]

    vjd = [2.46002818657856e6, 2.460028190050782e6]

    omm = parse_omm(
        """
        <?xml version="1.0" encoding="utf-8"?>
        <ndm><omm id="CCSDS_OMM_VERS" version="3.0">
        <header><CREATION_DATE>2025-12-30T23:36:37</CREATION_DATE><ORIGINATOR>18 SPCS</ORIGINATOR></header>
        <body><segment>
        <metadata><OBJECT_NAME>AMAZONIA 1</OBJECT_NAME><OBJECT_ID>2021-015A</OBJECT_ID><CENTER_NAME>EARTH</CENTER_NAME><REF_FRAME>TEME</REF_FRAME><TIME_SYSTEM>UTC</TIME_SYSTEM><MEAN_ELEMENT_THEORY>SGP4</MEAN_ELEMENT_THEORY></metadata>
        <data>
        <meanElements><EPOCH>2025-12-30T18:12:04.533984</EPOCH><MEAN_MOTION>14.40772474</MEAN_MOTION><ECCENTRICITY>0.00011240</ECCENTRICITY><INCLINATION>98.3721</INCLINATION><RA_OF_ASC_NODE>75.0877</RA_OF_ASC_NODE><ARG_OF_PERICENTER>97.3772</ARG_OF_PERICENTER><MEAN_ANOMALY>262.7545</MEAN_ANOMALY></meanElements>
        <tleParameters><EPHEMERIS_TYPE>0</EPHEMERIS_TYPE><CLASSIFICATION_TYPE>U</CLASSIFICATION_TYPE><NORAD_CAT_ID>47699</NORAD_CAT_ID><ELEMENT_SET_NO>999</ELEMENT_SET_NO><REV_AT_EPOCH>25439</REV_AT_EPOCH><BSTAR>0.00015330000000</BSTAR><MEAN_MOTION_DOT>0.00000447</MEAN_MOTION_DOT><MEAN_MOTION_DDOT>0.0000000000000</MEAN_MOTION_DDOT></tleParameters>
        </data>
        </segment></body>
        </omm></ndm>
        """,
    )

    # The fitting functions print their progress, which we silence here.
    redirect_stdout(devnull) do
        PrecompileTools.@compile_workload begin
            # Exercise, for every propagator and for both `Float64` and `Float32`, the
            # initialization, the propagation with every time representation and sink, the
            # epoch-based propagation, the stepping, and the mean elements fitting.
            for (prop, f32_kwargs) in (
                (:J2, (; j2c = J2C_EGM2008_F32)),
                (:J2osc, (; j2c = J2C_EGM2008_F32)),
                (:J4, (; j4c = J4C_EGM2008_F32)),
                (:J4osc, (; j4c = J4C_EGM2008_F32)),
                (:SGP4, (; sgp4c = sgp4c_wgs84_f32)),
                (:TwoBody, (; m0 = TBC_M0_F32)),
            )
                mean_elements = prop != :SGP4 ? orb : tle

                # == Float64 ===============================================================

                orbp = Propagators.init(Val(prop), mean_elements)

                Propagators.propagate!(orbp, 0.0)
                Propagators.propagate!(orbp, [0.0, 1.0])
                Propagators.propagate!(orbp, [0.0, 1.0, 2.0])

                Propagators.propagate!(orbp, 0.0, Tuple)
                Propagators.propagate!(orbp, [0.0, 1.0], Tuple)
                Propagators.propagate!(orbp, [0.0, 1.0, 2.0], Tuple)

                Propagators.propagate!(orbp, 0.0, OrbitStateVector)
                Propagators.propagate!(orbp, [0.0, 1.0], OrbitStateVector)
                Propagators.propagate!(orbp, [0.0, 1.0, 2.0], OrbitStateVector)

                Propagators.propagate!(orbp, Dates.Second(1))
                Propagators.propagate!(orbp, Dates.Second(1) + Dates.Minute(1))
                Propagators.propagate!(orbp, [Dates.Second(1) for _ in 1:2])
                Propagators.propagate!(
                    orbp, [Dates.Second(i) + Dates.Minute(1) for i in 1:2]
                )

                Propagators.propagate!(orbp, Dates.Second(1), Tuple)
                Propagators.propagate!(orbp, Dates.Second(1) + Dates.Minute(1), Tuple)
                Propagators.propagate!(orbp, [Dates.Second(1) for _ in 1:2], Tuple)
                Propagators.propagate!(
                    orbp, [Dates.Second(i) + Dates.Minute(1) for i in 1:2], Tuple
                )

                Propagators.propagate!(orbp, Dates.Second(1), OrbitStateVector)
                Propagators.propagate!(
                    orbp, Dates.Second(1) + Dates.Minute(1), OrbitStateVector
                )
                Propagators.propagate!(
                    orbp, [Dates.Second(1) for _ in 1:2], OrbitStateVector
                )
                Propagators.propagate!(
                    orbp, [Dates.Second(i) + Dates.Minute(1) for i in 1:2], OrbitStateVector
                )

                Propagators.propagate_to_epoch!(orbp, JD_J2000)
                Propagators.propagate_to_epoch!(orbp, [JD_J2000, JD_J2000])

                Propagators.propagate_to_epoch!(orbp, JD_J2000, Tuple)
                Propagators.propagate_to_epoch!(orbp, [JD_J2000, JD_J2000], Tuple)

                Propagators.propagate_to_epoch!(orbp, JD_J2000, OrbitStateVector)
                Propagators.propagate_to_epoch!(
                    orbp, [JD_J2000, JD_J2000], OrbitStateVector
                )

                Propagators.propagate_to_epoch!(orbp, DateTime(2024, 1, 1))
                Propagators.propagate_to_epoch!(orbp, [DateTime(2024, 1, i) for i in 1:2])

                Propagators.propagate_to_epoch!(orbp, DateTime(2024, 1, 1), Tuple)
                Propagators.propagate_to_epoch!(
                    orbp, [DateTime(2024, 1, i) for i in 1:2], Tuple
                )

                Propagators.propagate_to_epoch!(
                    orbp, DateTime(2024, 1, 1), OrbitStateVector
                )
                Propagators.propagate_to_epoch!(
                    orbp, [DateTime(2024, 1, i) for i in 1:2], OrbitStateVector
                )

                Propagators.step!(orbp, 1.0)

                Propagators.step!(orbp, 1.0, Tuple)

                Propagators.step!(orbp, 1.0, OrbitStateVector)

                Propagators.step!(orbp, Dates.Second(1))
                Propagators.step!(orbp, Dates.Second(1) + Dates.Minute(1))

                Propagators.step!(orbp, Dates.Second(1), Tuple)
                Propagators.step!(orbp, Dates.Second(1) + Dates.Minute(1), Tuple)

                Propagators.step!(orbp, Dates.Second(1), OrbitStateVector)
                Propagators.step!(orbp, Dates.Second(1) + Dates.Minute(1), OrbitStateVector)

                Propagators.fit_mean_elements!(orbp, vjd, vr_i, vv_i)
                Propagators.fit_mean_elements(Val(prop), vjd, vr_i, vv_i)

                # == Float32 ===============================================================

                orbp = Propagators.init(Val(prop), mean_elements; f32_kwargs...)

                Propagators.propagate!(orbp, 0.0f0)
                Propagators.propagate!(orbp, [0.0f0, 1.0f0])
                Propagators.propagate!(orbp, [0.0f0, 1.0f0, 2.0f0])

                Propagators.propagate!(orbp, 0.0f0, Tuple)
                Propagators.propagate!(orbp, [0.0f0, 1.0f0], Tuple)
                Propagators.propagate!(orbp, [0.0f0, 1.0f0, 2.0f0], Tuple)

                Propagators.propagate!(orbp, 0.0f0, OrbitStateVector)
                Propagators.propagate!(orbp, [0.0f0, 1.0f0], OrbitStateVector)
                Propagators.propagate!(orbp, [0.0f0, 1.0f0, 2.0f0], OrbitStateVector)

                Propagators.propagate!(orbp, Dates.Second(1))
                Propagators.propagate!(orbp, Dates.Second(1) + Dates.Minute(1))
                Propagators.propagate!(orbp, [Dates.Second(1) for _ in 1:2])
                Propagators.propagate!(
                    orbp, [Dates.Second(1) + Dates.Minute(1) for _ in 1:2]
                )

                Propagators.propagate!(orbp, Dates.Second(1), Tuple)
                Propagators.propagate!(orbp, Dates.Second(1) + Dates.Minute(1), Tuple)
                Propagators.propagate!(orbp, [Dates.Second(1) for _ in 1:2], Tuple)
                Propagators.propagate!(
                    orbp, [Dates.Second(1) + Dates.Minute(1) for _ in 1:2], Tuple
                )

                Propagators.propagate!(orbp, Dates.Second(1), OrbitStateVector)
                Propagators.propagate!(
                    orbp, Dates.Second(1) + Dates.Minute(1), OrbitStateVector
                )
                Propagators.propagate!(
                    orbp, [Dates.Second(1) for _ in 1:2], OrbitStateVector
                )
                Propagators.propagate!(
                    orbp, [Dates.Second(1) + Dates.Minute(1) for _ in 1:2], OrbitStateVector
                )

                Propagators.propagate_to_epoch!(orbp, JD_J2000)
                Propagators.propagate_to_epoch!(orbp, [JD_J2000, JD_J2000])

                Propagators.propagate_to_epoch!(orbp, JD_J2000, Tuple)
                Propagators.propagate_to_epoch!(orbp, [JD_J2000, JD_J2000], Tuple)

                Propagators.propagate_to_epoch!(orbp, JD_J2000, OrbitStateVector)
                Propagators.propagate_to_epoch!(
                    orbp, [JD_J2000, JD_J2000], OrbitStateVector
                )

                Propagators.propagate_to_epoch!(orbp, DateTime(2024, 1, 1))
                Propagators.propagate_to_epoch!(orbp, [DateTime(2024, 1, i) for i in 1:2])

                Propagators.propagate_to_epoch!(orbp, DateTime(2024, 1, 1), Tuple)
                Propagators.propagate_to_epoch!(
                    orbp, [DateTime(2024, 1, i) for i in 1:2], Tuple
                )

                Propagators.propagate_to_epoch!(
                    orbp, DateTime(2024, 1, 1), OrbitStateVector
                )
                Propagators.propagate_to_epoch!(
                    orbp, [DateTime(2024, 1, i) for i in 1:2], OrbitStateVector
                )

                Propagators.step!(orbp, 1.0f0)

                Propagators.step!(orbp, 1.0f0, Tuple)

                Propagators.step!(orbp, 1.0f0, OrbitStateVector)

                Propagators.step!(orbp, Dates.Second(1))
                Propagators.step!(orbp, Dates.Second(1) + Dates.Minute(1))

                Propagators.step!(orbp, Dates.Second(1), Tuple)
                Propagators.step!(orbp, Dates.Second(1) + Dates.Minute(1), Tuple)

                Propagators.step!(orbp, Dates.Second(1), OrbitStateVector)
                Propagators.step!(orbp, Dates.Second(1) + Dates.Minute(1), OrbitStateVector)
            end

            # Exercise the SGP4 initialization from an Orbit Mean-Elements Message.
            orbp = Propagators.init(Val(:SGP4), omm)
            Propagators.propagate!(orbp, 0.0)
            Propagators.init!(orbp, omm)

            orbp = Propagators.init(Val(:SGP4), omm; sgp4c = sgp4c_wgs84_f32)
            Propagators.propagate!(orbp, 0.0f0)
        end
    end
end
