## Description #############################################################################
#
#  API stubs for ML-∂SGP4 orbit propagator.
#
#  The real implementations are provided by SatelliteToolboxPropagatorsMLdSGP4Ext, which is
#  loaded automatically when Lux.jl, Optimisers.jl, and Zygote.jl are available.
#
############################################################################################

const _MLDSGP4_EXT_ERR = "ML-∂SGP4 requires Lux.jl, Optimisers.jl, and Zygote.jl. " *
                          "Load them with `using Lux, Optimisers, Zygote`."

"""
    Propagators.init(Val(:MLdSGP4), tle::TLE; kwargs...) -> OrbitPropagatorMLdSGP4

Create and initialize the ML-∂SGP4 orbit propagator structure using a TLE.

!!! note

    This function requires the Lux.jl extension to be loaded.  Make sure to run
    `using Lux, Optimisers, Zygote` before calling this function.

# Keywords

- `model`: Trained `MLdSGP4Model` providing the neural-network correction weights.
    If omitted, a zero-correction model is used (equivalent to plain SGP4).
"""
function Propagators.init(::Val{:MLdSGP4}, args...; kwargs...)
    error(_MLDSGP4_EXT_ERR)
end
