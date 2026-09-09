## Description #############################################################################
#
#  API implementation for J4 osculating orbit propagator.
#
############################################################################################

# Implement the `Propagators` API for the J4 osculating orbit propagator.
Propagators.epoch(orbp::OrbitPropagatorJ4Osculating)           = orbp.j4oscd.j4d.orb₀.epoch
Propagators.is_initialized(orbp::OrbitPropagatorJ4Osculating)  = _is_initialized(orbp.j4oscd)
Propagators.last_instant(orbp::OrbitPropagatorJ4Osculating)    = _last_instant(orbp.j4oscd)
Propagators.mean_elements(orbp::OrbitPropagatorJ4Osculating)   = orbp.j4oscd.j4d.orbk
Propagators.name(orbp::OrbitPropagatorJ4Osculating)            = "J4 Osculating Orbit Propagator"
Propagators.propagator_data(orbp::OrbitPropagatorJ4Osculating) = orbp.j4oscd

"""
    Propagators.fit_mean_elements(
        ::Val{:J4osc},
        vjd::AbstractVector{Tjd},
        vr_i::AbstractVector{Tv},
        vv_i::AbstractVector{Tv};
        kwargs...
    ) -> KeplerianElements{MeanAnomaly, Float64, T}, SMatrix{6, 6, T}, NamedTuple

Fit a set of mean Keplerian elements for the J4 osculating orbit propagator using the
osculating elements represented by a set of position vectors `vr_i` [m] and a set of
velocity vectors `vv_i` [m / s] represented in an inertial reference frame at instants in
the array `vjd` [Julian Day].

!!! note

    This algorithm version will allocate a new J4 osculating propagator with the constants
    `j4c`. If the allocation must be avoided, use the function
    [`Propagators.fit_mean_elements!`](@ref) instead.

# Keywords

- `j4c::J4PropagatorConstants{T}`: J4 orbit propagator constants (see
    [`J4PropagatorConstants`](@ref)), whose number type `T` is used in the fitting.
    (**Default**: `J4C_EGM2008`)

The other keywords are the same as in [`fit_j4osc_mean_elements`](@ref).

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
    ::Val{:J4osc},
    vjd::AbstractVector{Tjd},
    vr_i::AbstractVector{Tv},
    vv_i::AbstractVector{Tv};
    kwargs...,
) where {Tjd <: Number, Tv <: AbstractVector}
    return fit_j4osc_mean_elements(vjd, vr_i, vv_i; kwargs...)
end

"""
    Propagators.fit_mean_elements!(
        orbp::OrbitPropagatorJ4Osculating,
        vjd::AbstractVector{Tjd},
        vr_i::AbstractVector{Tv},
        vv_i::AbstractVector{Tv};
        kwargs...
    ) where {
        Tjd <: Number,
        Tv <: AbstractVector
    } -> KeplerianElements{MeanAnomaly, Tepoch, T}, SMatrix{6, 6, T}, NamedTuple

Fit a set of mean Keplerian elements for the J4 osculating orbit propagator `orbp` using the
osculating elements represented by a set of position vectors `vr_i` [m] and a set of
velocity vectors `vv_i` [m / s] represented in an inertial reference frame at instants in
the array `vjd` [Julian Day].

!!! note

    The orbit propagator `orbp` will be initialized with the Keplerian elements returned
    by the function.

# Keywords

The keywords are the same as in [`fit_j4osc_mean_elements`](@ref), except for `j4c`, since
the constants are those in `orbp`.

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
    orbp::OrbitPropagatorJ4Osculating,
    vjd::AbstractVector{Tjd},
    vr_i::AbstractVector{Tv},
    vv_i::AbstractVector{Tv};
    kwargs...,
) where {Tjd <: Number, Tv <: AbstractVector}
    return fit_j4osc_mean_elements!(orbp.j4oscd, vjd, vr_i, vv_i; kwargs...)
end

"""
    Propagators.init(
        Val(:J4osc),
        orb₀::KeplerianElements;
        kwargs...
    ) -> OrbitPropagatorJ4Osculating

Create and initialize the J4 osculating orbit propagator structure using the mean Keplerian
elements `orb₀` [SI units].

!!! note

    The type used in the propagation will be the same as used to define the constants in the
    structure `j4c`.

# Keywords

- `j4c::J4PropagatorConstants`: J4 orbit propagator constants (see
    [`J4PropagatorConstants`](@ref)).
    (**Default**: `J4C_EGM2008`)
"""
function Propagators.init(
    ::Val{:J4osc}, orb₀::KeplerianElements; j4c::J4PropagatorConstants = J4C_EGM2008
)
    j4oscd = j4osc_init(orb₀; j4c = j4c)
    return OrbitPropagatorJ4Osculating(j4oscd)
end

"""
    Propagators.init!(orbp::OrbitPropagatorJ4Osculating, orb₀::KeplerianElements) -> Nothing

Initialize the J4 osculating orbit propagator structure `orbp` using the mean Keplerian
elements `orb₀` [SI units].

!!! warning

    The propagation constants `j4c::J4PropagatorConstants` in `orbp.j4oscd.j4d` will not be
    changed. Hence, they must be initialized.
"""
function Propagators.init!(orbp::OrbitPropagatorJ4Osculating, orb₀::KeplerianElements)
    j4osc_init!(orbp.j4oscd, orb₀)
    return nothing
end

"""
    Propagators.propagate!(
        orbp::OrbitPropagatorJ4Osculating{Tepoch, T},
        Δt::Number
    ) where {Tepoch <: Number, T <: Number} -> SVector{3, T}, SVector{3, T}

Propagate the orbit of the J4 osculating orbit propagator `orbp` to `Δt` [s] after the epoch
of the initial mean elements, updating the internal state of `orbp`.

# Returns

- `SVector{3, T}`: Position vector [m] represented in the inertial frame at propagation
    instant.
- `SVector{3, T}`: Velocity vector [m / s] represented in the inertial frame at propagation
    instant.

# Remarks

The output is represented in the inertial reference frame of the input elements. The
perturbation theory requires an inertial frame with true equator.
"""
function Propagators.propagate!(orbp::OrbitPropagatorJ4Osculating, Δt::Number)
    return j4osc!(orbp.j4oscd, Δt)
end

############################################################################################
#                                        Julia API                                         #
############################################################################################

function Base.copy(
    orbp::OrbitPropagatorJ4Osculating{Tepoch, T}
) where {Tepoch <: Number, T <: Number}
    return OrbitPropagatorJ4Osculating{Tepoch, T}(copy(orbp.j4oscd))
end
