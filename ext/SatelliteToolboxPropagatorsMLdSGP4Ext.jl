module SatelliteToolboxPropagatorsMLdSGP4Ext

using SatelliteToolboxPropagators
using MLdSGP4

import SatelliteToolboxPropagators: Propagators, OrbitPropagatorMLdSGP4

# ==========================================================================================
#                                   Propagators API
# ==========================================================================================

"""
    Propagators.init(Val(:MLdSGP4), tle::TLE; kwargs...) -> OrbitPropagatorMLdSGP4

Create and initialize the ML-∂SGP4 orbit propagator structure using a TLE.

# Keywords

- `sgp4c::Sgp4Constants`: SGP4 orbit propagator constants (see `Sgp4Constants`).
    (**Default** = `sgp4c_wgs84`)
- `model`: Trained `MLdSGP4Model` providing the neural-network correction weights.
    If omitted, a zero-correction model is used (equivalent to plain SGP4).
"""
function Propagators.init(::Val{:MLdSGP4}, tle::TLE; kwargs...)
    mlsgp4d = ml_dsgp4_init(tle; kwargs...)
    sgp4d   = mlsgp4d.sgp4d
    Tepoch  = typeof(sgp4d).parameters[1]
    T       = typeof(sgp4d).parameters[2]
    return OrbitPropagatorMLdSGP4{Tepoch, T, typeof(mlsgp4d)}(mlsgp4d)
end

Propagators.epoch(orbp::OrbitPropagatorMLdSGP4)        = orbp.mlsgp4d.epoch
Propagators.last_instant(orbp::OrbitPropagatorMLdSGP4) = orbp.mlsgp4d.sgp4d.Δt * 60
Propagators.name(orbp::OrbitPropagatorMLdSGP4)         = "ML-∂SGP4 Orbit Propagator"

function Propagators.propagate!(orbp::OrbitPropagatorMLdSGP4, t::Number)
    r_i, v_i = ml_dsgp4!(orbp.mlsgp4d, t / 60)
    return 1000r_i, 1000v_i
end

end # module SatelliteToolboxPropagatorsMLdSGP4Ext
