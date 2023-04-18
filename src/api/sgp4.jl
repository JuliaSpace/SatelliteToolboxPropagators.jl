# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#
# Description
# ==========================================================================================
#
#    API implementation for SGP4 orbit propagator.
#
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

Propagators.epoch(orbp::OrbitPropagatorSgp4)         = orbp.sgp4d.epoch
Propagators.last_instant(orbp::OrbitPropagatorSgp4)  = orbp.sgp4d.Δt * 60
Propagators.name(orbp::OrbitPropagatorSgp4)          = "SGP4 Orbit Propagator"

function Propagators.mean_elements(orbp::OrbitPropagatorSgp4)
    sgp4d = orbp.sgp4d

    return KeplerianElements(
        sgp4d.epoch + sgp4d.Δt / 86400,
        sgp4d.a_k * sgp4d.sgp4c.R0,
        sgp4d.e_k,
        sgp4d.i_k,
        sgp4d.Ω_k,
        sgp4d.ω_k,
        mean_to_true_anomaly(sgp4d.e_k, sgp4d.M_k)
    )
end

"""
    Propagators.init(Val(:SGP4), epoch::Number, n₀::Number, e₀::Number, i₀::Number, Ω₀::Number, ω₀::Number, M₀::Number, bstar::Number; kwargs...) -> OrbitPropagatorSgp4
    Propagators.init(Val(:SGP4), tle::TLE; kwargs...) -> OrbitPropagatorSgp4

Create and initialize the SGP4 orbit propagator structure using the initial orbit specified
by the arguments.

!!! note
    The type used in the propagation will be the same as used to define the constants in the
    structure `sgp4c`.

# Arguments

- `epoch::Number`: Epoch of the orbital elements [Julian Day].
- `n₀::Number`: SGP type "mean" mean motion at epoch [rad/s].
- `e₀::Number`: "Mean" eccentricity at epoch.
- `i₀::Number`: "Mean" inclination at epoch [rad].
- `Ω₀::Number`: "Mean" longitude of the ascending node at epoch [rad].
- `ω₀::Number`: "Mean" argument of perigee at epoch [rad].
- `M₀::Number`: "Mean" mean anomaly at epoch [rad].
- `bstar::Number`: Drag parameter (B*).
- `tle::TLE`: Two-line elements used for the initialization.

# Keywords

- `sgp4c::Sgp4Constants`: SGP4 orbit propagator constants (see [`Sgp4Constants`](@ref)).
    (**Default** = `sgp4c_wgs84`)
"""
function Propagators.init(::Val{:SGP4}, tle::TLE; sgp4c::Sgp4Constants = sgp4c_wgs84)
    sgp4d = sgp4_init(tle; sgp4c = sgp4c)
    return OrbitPropagatorSgp4(sgp4d)
end

function Propagators.init(
    ::Val{:SGP4},
    epoch::Number,
    n₀::Number,
    e₀::Number,
    i₀::Number,
    Ω₀::Number,
    ω₀::Number,
    M₀::Number,
    bstar::Number;
    sgp4c::Sgp4Constants = sgp4c_wgs84
)
    sgp4d = sgp4_init(epoch, 60n₀, e₀, i₀, Ω₀, ω₀, M₀, bstar; sgp4c = sgp4c)
    return OrbitPropagatorSgp4(sgp4d)
end

"""
    Propagators.init!(orbp::OrbitPropagatorSgp4, epoch::Number, n₀::Number, e₀::Number, i₀::Number, Ω₀::Number, ω₀::Number, M₀::Number, bstar::Number; kwargs...) -> Nothing
    Propagators.init!(orbp::OrbitPropagatorSgp4, tle::TLE; kwargs...) -> Nothing

Initialize the SGP4 orbit propagator structure `orbp` using the initial orbit specified
by the arguments.

!!! warning
    The propagation constants `sgp4c::Sgp4Constants` in `orbp.sgp4d` will not be changed.
    Hence, they must be initialized.

# Arguments

- `epoch::Number`: Epoch of the orbital elements [Julian Day].
- `n₀::Number`: SGP type "mean" mean motion at epoch [rad/s].
- `e₀::Number`: "Mean" eccentricity at epoch.
- `i₀::Number`: "Mean" inclination at epoch [rad].
- `Ω₀::Number`: "Mean" longitude of the ascending node at epoch [rad].
- `ω₀::Number`: "Mean" argument of perigee at epoch [rad].
- `M₀::Number`: "Mean" mean anomaly at epoch [rad].
- `bstar::Number`: Drag parameter (B*).
- `tle::TLE`: Two-line elements used for the initialization.
"""
function Propagators.init!(orbp::OrbitPropagatorSgp4, tle::TLE)
    sgp4_init!(orbp.sgp4d, tle)
    return nothing
end

function Propagators.init!(
    orbp::OrbitPropagatorSgp4,
    epoch::Number,
    n₀::Number,
    e₀::Number,
    i₀::Number,
    Ω₀::Number,
    ω₀::Number,
    M₀::Number,
    bstar::Number
)
    sgp4_init!(orbp.sgp4d, epoch, 60n₀, e₀, i₀, Ω₀, ω₀, M₀, bstar)
    return nothing
end

"""
    Propagators.propagate(Val(:SGP4), Δt::Number, epoch::Number, n₀::Number, e₀::Number, i₀::Number, Ω₀::Number, ω₀::Number, M₀::Number, bstar::Number; kwargs...)
    Propagators.propagate(Val(:SGP4), Δt::Number, tle::TLE; kwargs...)

Initialize the SGP4 propagator structure using the input arguments and propagate the orbit
until the time Δt [s].

!!! note
    The type used in the propagation will be the same as used to define the constants in the
    structure `sgp4c`.

# Arguments

- `epoch::Number`: Epoch of the orbital elements [Julian Day].
- `n₀::Number`: SGP type "mean" mean motion at epoch [rad/s].
- `e₀::Number`: "Mean" eccentricity at epoch.
- `i₀::Number`: "Mean" inclination at epoch [rad].
- `Ω₀::Number`: "Mean" longitude of the ascending node at epoch [rad].
- `ω₀::Number`: "Mean" argument of perigee at epoch [rad].
- `M₀::Number`: "Mean" mean anomaly at epoch [rad].
- `bstar::Number`: Drag parameter (B*).
- `tle::TLE`: Two-line elements used for the initialization.

# Keywords

- `sgp4c::Sgp4Constants{T}`: SGP4 orbit propagator constants (see [`Sgp4Constants`](@ref)).
    (**Default** = `sgp4c_wgs84`)

# Returns

- `SVector{3, T}`: Position vector [m] represented in the inertial frame at propagation
    instant.
- `SVector{3, T}`: Velocity vector [m / s] represented in the inertial frame at propagation
    instant.
- [`OrbitPropagatorSgp4`](@ref): Structure with the initialized parameters.
"""
function Propagators.propagate(
    ::Val{:SGP4},
    Δt::Number,
    tle::TLE;
    sgp4c::Sgp4Constants = sgp4c_wgs84
)
    r_i, v_i, sgp4d = sgp4(Δt / 60, tle; sgp4c = sgp4c)
    return 1000r_i, 1000v_i, OrbitPropagatorSgp4(sgp4d)
end

function Propagators.propagate(
    ::Val{:SGP4},
    Δt::Number,
    epoch::Number,
    n₀::Number,
    e₀::Number,
    i₀::Number,
    Ω₀::Number,
    ω₀::Number,
    M₀::Number,
    bstar::Number;
    sgp4c::Sgp4Constants = sgp4c_wgs84
)
    r_i, v_i, sgp4d = sgp4(
        Δt / 60,
        epoch,
        60n₀,
        e₀,
        i₀,
        Ω₀,
        ω₀,
        M₀,
        bstar;
        sgp4c = sgp4c
    )
    return 1000r_i, 1000v_i, OrbitPropagatorSgp4(sgp4d)
end

function Propagators.propagate!(orbp::OrbitPropagatorSgp4, t::Number)
    # Auxiliary variables.
    sgp4d = orbp.sgp4d

    # Propagate the orbit.
    r_i, v_i = sgp4!(sgp4d, t / 60)

    return 1000r_i, 1000v_i
end
