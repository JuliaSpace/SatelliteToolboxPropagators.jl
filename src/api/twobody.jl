## Description #############################################################################
#
#  API implementation for two-body orbit propagator.
#
############################################################################################

# Implement the `Propagators` API for the two-body orbit propagator.
Propagators.epoch(orbp::OrbitPropagatorTwoBody)           = orbp.tbd.orb₀.epoch
Propagators.is_initialized(orbp::OrbitPropagatorTwoBody)  = _is_initialized(orbp.tbd)
Propagators.last_instant(orbp::OrbitPropagatorTwoBody)    = orbp.tbd.Δt
Propagators.mean_elements(orbp::OrbitPropagatorTwoBody)   = orbp.tbd.orbk
Propagators.name(orbp::OrbitPropagatorTwoBody)            = "Two-Body Orbit Propagator"
Propagators.propagator_data(orbp::OrbitPropagatorTwoBody) = orbp.tbd

"""
    Propagators.fit_mean_elements(
        ::Val{:TwoBody},
        vjd::AbstractVector{Tjd},
        vr_i::AbstractVector{Tv},
        vv_i::AbstractVector{Tv};
        kwargs...
    ) where {
        Tjd <: Number,
        Tv <: AbstractVector
    } -> KeplerianElements{MeanAnomaly, Float64, T}, SMatrix{6, 6, T}, NamedTuple

Fit a set of mean Keplerian elements for the two-body orbit propagator using the osculating
elements represented by a set of position vectors `vr_i` [m] and a set of velocity vectors
`vv_i` [m / s] represented in an inertial reference frame at instants in the array `vjd`
[Julian Day].

!!! note

    This algorithm version will allocate a new two-body propagator with the gravitational
    parameter `m0`. If the allocation must be avoided, use the function
    [`Propagators.fit_mean_elements!`](@ref) instead.

# Keywords

- `m0::T`: Standard gravitational parameter of the central body [m³ / s²], whose
    number type `T` is used in the fitting.
    (**Default**: `TBC_M0`)

The other keywords are the same as in [`fit_twobody_mean_elements`](@ref).

# Returns

- `KeplerianElements{MeanAnomaly, Float64, T}`: Fitted Keplerian elements.
- `SMatrix{6, 6, T}`: Final covariance matrix of the least-square algorithm.
- `NamedTuple`: Statistics of the least-square algorithm with the following fields:
    - `converged::Bool`: `true` if the iterations stopped because the residue was lower
        than `atol` or its relative variation was lower than `rtol`, or `false` if they
        stopped by reaching `max_iterations`.
    - `iterations::Int`: Number of iterations performed.
    - `position_rmse::T`: RMSE of the position residue in the last iteration [m].
    - `velocity_rmse::T`: RMSE of the velocity residue in the last iteration [m / s].
    - `total_rmse::T`: Weighted RMSE of the residue in the last iteration.

    The statistics refer to the fitting of the mean elements. If their epoch is updated
    afterward to match `mean_elements_epoch`, the statistics of that update are not
    returned.

# Extended help

## Throws

- `ArgumentError`: If `vjd`, `vr_i`, and `vv_i` do not have the same length, if
    `weight_vector` does not have six elements, or if `max_iterations` is lower than 1.
- `MeanElementsFitDivergenceError`: If the least-square iterations diverge.
"""
function Propagators.fit_mean_elements(
    ::Val{:TwoBody},
    vjd::AbstractVector{Tjd},
    vr_i::AbstractVector{Tv},
    vv_i::AbstractVector{Tv};
    kwargs...,
) where {Tjd <: Number, Tv <: AbstractVector}
    return fit_twobody_mean_elements(vjd, vr_i, vv_i; kwargs...)
end

"""
    Propagators.fit_mean_elements!(
        orbp::OrbitPropagatorTwoBody,
        vjd::AbstractVector{Tjd},
        vr_i::AbstractVector{Tv},
        vv_i::AbstractVector{Tv};
        kwargs...
    ) where {
        Tjd <: Number,
        Tv <: AbstractVector
    } -> KeplerianElements{MeanAnomaly, Tepoch, T}, SMatrix{6, 6, T}, NamedTuple

Fit a set of mean Keplerian elements for the two-body orbit propagator `orbp` using the
osculating elements represented by a set of position vectors `vr_i` [m] and a set of
velocity vectors `vv_i` [m / s] represented in an inertial reference frame at instants in
the array `vjd` [Julian Day].

!!! note

    The orbit propagator `orbp` will be initialized with the Keplerian elements returned
    by the function.

# Keywords

The keywords are the same as in [`fit_twobody_mean_elements`](@ref), except for `m0`, since
the gravitational parameter is that in `orbp`.

# Returns

- `KeplerianElements{MeanAnomaly, Tepoch, T}`: Fitted Keplerian elements.
- `SMatrix{6, 6, T}`: Final covariance matrix of the least-square algorithm.
- `NamedTuple`: Statistics of the least-square algorithm with the following fields:
    - `converged::Bool`: `true` if the iterations stopped because the residue was lower
        than `atol` or its relative variation was lower than `rtol`, or `false` if they
        stopped by reaching `max_iterations`.
    - `iterations::Int`: Number of iterations performed.
    - `position_rmse::T`: RMSE of the position residue in the last iteration [m].
    - `velocity_rmse::T`: RMSE of the velocity residue in the last iteration [m / s].
    - `total_rmse::T`: Weighted RMSE of the residue in the last iteration.

    The statistics refer to the fitting of the mean elements. If their epoch is updated
    afterward to match `mean_elements_epoch`, the statistics of that update are not
    returned.

# Extended help

## Throws

- `ArgumentError`: If `vjd`, `vr_i`, and `vv_i` do not have the same length, if
    `weight_vector` does not have six elements, or if `max_iterations` is lower than 1.
- `MeanElementsFitDivergenceError`: If the least-square iterations diverge.
"""
function Propagators.fit_mean_elements!(
    orbp::OrbitPropagatorTwoBody,
    vjd::AbstractVector{Tjd},
    vr_i::AbstractVector{Tv},
    vv_i::AbstractVector{Tv};
    kwargs...,
) where {Tjd <: Number, Tv <: AbstractVector}
    return fit_twobody_mean_elements!(orbp.tbd, vjd, vr_i, vv_i; kwargs...)
end

"""
    Propagators.init(
        Val(:TwoBody),
        orb₀::KeplerianElements;
        kwargs...
    ) -> OrbitPropagatorTwoBody

Create and initialize the two-body orbit propagator structure using the mean Keplerian
elements `orb₀` [SI units].

!!! note

    The type used in the propagation will be the same as used to define the standard
    gravitational parameter `m0`.

# Keywords

- `m0::T`: Standard gravitational parameter of the central body [m³ / s²].
    (**Default**: `TBC_M0`)
"""
function Propagators.init(::Val{:TwoBody}, orb₀::KeplerianElements; m0::Number = TBC_M0)
    tbd = twobody_init(orb₀; m0 = m0)
    return OrbitPropagatorTwoBody(tbd)
end

"""
    Propagators.init!(orbp::OrbitPropagatorTwoBody, orb₀::KeplerianElements) -> Nothing

Initialize the two-body orbit propagator structure `orbp` using the mean Keplerian elements
`orb₀` [SI units].

!!! warning

    The propagation constant `μ::T` in `orbp.tbd`, set from the keyword `m0`, will not be
    changed. Hence, it must be initialized.

# Arguments

- `orb₀::KeplerianElements`: Initial mean Keplerian elements [SI units].
"""
function Propagators.init!(orbp::OrbitPropagatorTwoBody, orb₀::KeplerianElements)
    twobody_init!(orbp.tbd, orb₀)
    return nothing
end

"""
    Propagators.propagate!(
        orbp::OrbitPropagatorTwoBody{Tepoch, T},
        Δt::Number
    ) where {Tepoch <: Number, T <: Number} -> SVector{3, T}, SVector{3, T}

Propagate the orbit of the two-body orbit propagator `orbp` to `Δt` [s] after the epoch of
the initial mean elements, updating the internal state of `orbp`.

# Returns

- `SVector{3, T}`: Position vector [m] represented in the inertial frame at propagation
    instant.
- `SVector{3, T}`: Velocity vector [m / s] represented in the inertial frame at propagation
    instant.

# Remarks

The output is represented in the inertial reference frame of the input elements.
"""
function Propagators.propagate!(orbp::OrbitPropagatorTwoBody, Δt::Number)
    return twobody!(orbp.tbd, Δt)
end

############################################################################################
#                                        Julia API                                         #
############################################################################################

function Base.copy(
    orbp::OrbitPropagatorTwoBody{Tepoch, T}
) where {Tepoch <: Number, T <: Number}
    return OrbitPropagatorTwoBody{Tepoch, T}(copy(orbp.tbd))
end
