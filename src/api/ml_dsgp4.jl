## Description #############################################################################
#
#  API stubs for ML-∂SGP4 orbit propagator.
#
#  The real implementations are provided by SatelliteToolboxPropagatorsMLdSGP4Ext, which is
#  loaded automatically when MLdSGP4.jl is available.
#
############################################################################################

const _MLDSGP4_EXT_ERR = "ML-∂SGP4 requires MLdSGP4.jl. Load it with `using MLdSGP4`."

"""
    Propagators.init(Val(:MLdSGP4), tle::TLE; kwargs...) -> OrbitPropagatorMLdSGP4

Create and initialize the ML-∂SGP4 orbit propagator structure using a TLE.

!!! note

    This function requires the MLdSGP4.jl package to be loaded.  Make sure to run
    `using MLdSGP4` before calling this function.

# Keywords

- `model`: Trained `MLdSGP4Model` providing the neural-network correction weights.
    If omitted, a zero-correction model is used (equivalent to plain SGP4).
"""
function Propagators.init(::Val{:MLdSGP4}, args...; kwargs...)
    error(_MLDSGP4_EXT_ERR)
end
