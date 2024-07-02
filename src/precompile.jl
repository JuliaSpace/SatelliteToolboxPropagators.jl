## Description #############################################################################
#
#   Precompilation.
#
############################################################################################

import PrecompileTools

PrecompileTools.@setup_workload begin
    orb = KeplerianElements(
        DateTime("2023-01-01") |> datetime2julian,
        Float64(8000e3),
        Float64(0.015),
        Float64(28.5) |> deg2rad,
        Float64(100)  |> deg2rad,
        Float64(200)  |> deg2rad,
        Float64(45)   |> deg2rad
    )

    tle = tle"""
          AMAZONIA 1
          1 47699U 21015A   23083.68657856  .00000000  00000-8  43000-3 0  9999
          2 47699  98.4304 162.1097 0001247 136.2017 223.9283 14.40814394108652
          """

    vr_i = [
        @SVector([-6792.402703741442, 2192.6458461287293, 0.18851758695295118]) .* 1000,
        @SVector([-6357.88873265975, 2391.9476768911686, 2181.838771262736]) .* 1000
    ]

    vv_i = [
        @SVector([0.3445760107690598, 1.0395135806993514, 7.393686131436984]) .* 1000,
        @SVector([2.5285015912807003, 0.27812476784300005, 7.030323100703928]) .* 1000
    ]

    vjd = [
        2.46002818657856e6,
        2.460028190050782e6
    ]

    redirect_stdout(devnull) do
        PrecompileTools.@compile_workload begin
            for (prop, f32_kwargs) in (
                (:J2,      (; j2c = j2c_egm2008_f32)),
                (:J2osc,   (; j2c = j2c_egm2008_f32)),
                (:J4,      (; j4c = j4c_egm2008_f32)),
                (:J4osc,   (; j4c = j4c_egm2008_f32)),
                (:SGP4,    (; sgp4c = sgp4c_wgs84_f32)),
                (:TwoBody, (; m0 = tbc_m0_f32))
            )
                mean_elements = prop != :SGP4 ? orb : tle
                orbp = Propagators.init(Val(prop), mean_elements)

                Propagators.propagate!(orbp, 0.0)
                Propagators.propagate!(orbp, [0.0, 1.0])
                Propagators.propagate!(orbp, [0.0, 1.0, 2.0])

                if prop != :TwoBody
                    Propagators.fit_mean_elements!(orbp, vjd, vr_i, vv_i)
                    Propagators.fit_mean_elements(Val(prop), vjd, vr_i, vv_i)
                end

                orbp = Propagators.init(Val(prop), mean_elements; f32_kwargs...)
                Propagators.propagate!(orbp, 0.0f0)
                Propagators.propagate!(orbp, [0.0f0, 1.0f0])
                Propagators.propagate!(orbp, [0.0f0, 1.0f0, 2.0f0])
            end
        end
    end
end
