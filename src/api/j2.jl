## Description #############################################################################
#
#  API implementation for J2 orbit propagator.
#
############################################################################################

# Implement the `Propagators` API for the J2 orbit propagator.
Propagators.epoch(orbp::OrbitPropagatorJ2)           = orbp.j2d.orb₀.epoch
Propagators.is_initialized(orbp::OrbitPropagatorJ2)  = _is_initialized(orbp.j2d)
Propagators.last_instant(orbp::OrbitPropagatorJ2)    = orbp.j2d.Δt
Propagators.mean_elements(orbp::OrbitPropagatorJ2)   = orbp.j2d.orbk
Propagators.name(orbp::OrbitPropagatorJ2)            = "J2 Orbit Propagator"
Propagators.propagator_data(orbp::OrbitPropagatorJ2) = orbp.j2d

"""
    Propagators.fit_mean_elements(
        ::Val{:J2},
        vjd::AbstractVector{Tjd},
        vr_i::AbstractVector{Tv},
        vv_i::AbstractVector{Tv};
        kwargs...
    ) -> KeplerianElements{MeanAnomaly, Float64, T}, SMatrix{6, 6, T}, NamedTuple

Fit a set of mean Keplerian elements for the J2 orbit propagator using the osculating
elements represented by a set of position vectors `vr_i` [m] and a set of velocity vectors
`vv_i` [m / s] represented in an inertial reference frame at instants in the array `vjd`
[Julian Day].

!!! note

    This algorithm version will allocate a new J2 propagator with the constants `j2c`. If
    the allocation must be avoided, use the function
    [`Propagators.fit_mean_elements!`](@ref) instead.

# Keywords

- `j2c::J2PropagatorConstants{T}`: J2 orbit propagator constants (see
    [`J2PropagatorConstants`](@ref)), whose number type `T` is used in the fitting.
    (**Default**: `J2C_EGM2008`)

The other keywords are the same as in [`fit_j2_mean_elements`](@ref).

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
    ::Val{:J2},
    vjd::AbstractVector{Tjd},
    vr_i::AbstractVector{Tv},
    vv_i::AbstractVector{Tv};
    kwargs...,
) where {Tjd <: Number, Tv <: AbstractVector}
    return fit_j2_mean_elements(vjd, vr_i, vv_i; kwargs...)
end

"""
    Propagators.fit_mean_elements!(
        orbp::OrbitPropagatorJ2,
        vjd::AbstractVector{Tjd},
        vr_i::AbstractVector{Tv},
        vv_i::AbstractVector{Tv};
        kwargs...
    ) where {
        Tjd <: Number,
        Tv <: AbstractVector
    } -> KeplerianElements{MeanAnomaly, Tepoch, T}, SMatrix{6, 6, T}, NamedTuple

Fit a set of mean Keplerian elements for the J2 orbit propagator `orbp` using the osculating
elements represented by a set of position vectors `vr_i` [m] and a set of velocity vectors
`vv_i` [m / s] represented in an inertial reference frame at instants in the array `vjd`
[Julian Day].

!!! note

    The orbit propagator `orbp` will be initialized with the Keplerian elements returned
    by the function.

# Keywords

The keywords are the same as in [`fit_j2_mean_elements`](@ref), except for `j2c`, since the
constants are those in `orbp`.

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
    orbp::OrbitPropagatorJ2,
    vjd::AbstractVector{Tjd},
    vr_i::AbstractVector{Tv},
    vv_i::AbstractVector{Tv};
    kwargs...,
) where {Tjd <: Number, Tv <: AbstractVector}
    return fit_j2_mean_elements!(orbp.j2d, vjd, vr_i, vv_i; kwargs...)
end

"""
    Propagators.init(Val(:J2), orb₀::KeplerianElements; kwargs...) -> OrbitPropagatorJ2

Create and initialize the J2 orbit propagator structure using the mean Keplerian elements
`orb₀` [SI units].

!!! note

    The type used in the propagation will be the same as used to define the constants in the
    structure `j2c`.

# Keywords

- `j2c::J2PropagatorConstants`: J2 orbit propagator constants (see
    [`J2PropagatorConstants`](@ref)).
    (**Default**: `J2C_EGM2008`)
"""
function Propagators.init(
    ::Val{:J2}, orb₀::KeplerianElements; j2c::J2PropagatorConstants = J2C_EGM2008
)
    j2d = j2_init(orb₀; j2c = j2c)
    return OrbitPropagatorJ2(j2d)
end

"""
    Propagators.init!(orbp::OrbitPropagatorJ2, orb₀::KeplerianElements) -> Nothing

Initialize the J2 orbit propagator structure `orbp` using the mean Keplerian elements
`orb₀` [SI units].

!!! warning

    The propagation constants `j2c::J2PropagatorConstants` in `orbp.j2d` will not be
    changed. Hence, they must be initialized.
"""
function Propagators.init!(orbp::OrbitPropagatorJ2, orb₀::KeplerianElements)
    j2_init!(orbp.j2d, orb₀)
    return nothing
end

"""
    Propagators.propagate!(
        orbp::OrbitPropagatorJ2{Tepoch, T},
        Δt::Number
    ) where {Tepoch <: Number, T <: Number} -> SVector{3, T}, SVector{3, T}

Propagate the orbit of the J2 orbit propagator `orbp` to `Δt` [s] after the epoch of
the initial mean elements, updating the internal state of `orbp`.

# Returns

- `SVector{3, T}`: Position vector [m] represented in the inertial frame at propagation
    instant.
- `SVector{3, T}`: Velocity vector [m / s] represented in the inertial frame at propagation
    instant.

# Remarks

The output is represented in the inertial reference frame of the input elements. The
perturbation theory requires an inertial frame with true equator.
"""
function Propagators.propagate!(orbp::OrbitPropagatorJ2, Δt::Number)
    return j2!(orbp.j2d, Δt)
end

############################################################################################
#                                        Julia API                                         #
############################################################################################

function Base.copy(orbp::OrbitPropagatorJ2{Tepoch, T}) where {Tepoch <: Number, T <: Number}
    return OrbitPropagatorJ2{Tepoch, T}(copy(orbp.j2d))
end
