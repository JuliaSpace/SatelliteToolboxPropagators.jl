## Description #############################################################################
#
#  API implementation for J2 osculating orbit propagator.
#
############################################################################################

# Implement the `Propagators` API for the J2 osculating orbit propagator.
Propagators.epoch(orbp::OrbitPropagatorJ2Osculating) = orbp.j2oscd.j2d.orb₀.epoch
Propagators.is_initialized(orbp::OrbitPropagatorJ2Osculating) = _is_initialized(orbp.j2oscd)
Propagators.last_instant(orbp::OrbitPropagatorJ2Osculating) = orbp.j2oscd.Δt
Propagators.mean_elements(orbp::OrbitPropagatorJ2Osculating) = orbp.j2oscd.j2d.orbk
Propagators.name(orbp::OrbitPropagatorJ2Osculating) = "J2 Osculating Orbit Propagator"
Propagators.propagator_data(orbp::OrbitPropagatorJ2Osculating) = orbp.j2oscd

"""
    Propagators.fit_mean_elements(
        ::Val{:J2osc},
        vjd::AbstractVector{Tjd},
        vr_i::AbstractVector{Tv},
        vv_i::AbstractVector{Tv};
        kwargs...
    ) -> KeplerianElements{MeanAnomaly, Float64, Float64}, SMatrix{6, 6, Float64}

Fit a set of mean Keplerian elements for the J2 osculating orbit propagator using the
osculating elements represented by a set of position vectors `vr_i` [m] and a set of
velocity vectors `vv_i` [m / s] represented in an inertial reference frame at instants in
the array `vjd` [Julian Day].

!!! note

    This algorithm version will allocate a new J2 osculating propagator with the default
    constants `J2C_EGM2008`. If another set of constants are required, use the function
    [`Propagators.fit_mean_elements!`](@ref) instead.

# Keywords

- `atol::Number`: Tolerance for the residual absolute value. If the residual is lower than
    `atol` at any iteration, the computation loop stops.
    (**Default**: 2e-4)
- `rtol::Number`: Tolerance for the relative difference between the residuals. If the
    relative difference between the residuals in two consecutive iterations is lower than
    `rtol`, the computation loop stops.
    (**Default**: 2e-4)
- `initial_guess::Union{Nothing, KeplerianElements}`: Initial guess for the mean elements
    fitting process. If it is `nothing`, the algorithm will obtain an initial estimate from
    the osculating elements in `vr_i` and `vv_i`.
    (**Default**: nothing)
- `jacobian_method::AbstractJacobianMethod`: Method used to compute the Jacobian matrix. It
    can be `FiniteDiffJacobian()` for finite differences or `ForwardDiffJacobian()` for
    `ForwardDiff.jl` automatic differentiation.
    (**Default**: `FiniteDiffJacobian()`)
- `jacobian_perturbation::Number`: Initial state perturbation to compute the
    finite-difference when calculating the Jacobian matrix.
    (**Default**: 1e-3)
- `jacobian_perturbation_tol::Number`: Tolerance to accept the perturbation when calculating
    the Jacobian matrix. If the computed perturbation is lower than
    `jacobian_perturbation_tol`, we increase it until its absolute value is higher than
    `jacobian_perturbation_tol`.
    (**Default**: 1e-7)
- `max_iterations::Int`: Maximum number of iterations allowed for the least-square fitting.
    (**Default**: 50)
- `mean_elements_epoch::Number`: Epoch for the fitted mean elements.
    (**Default**: vjd[end])
- `verbose::Bool`: If `true`, the algorithm prints debugging information to `stdout`.
    (**Default**: true)
- `weight_vector::AbstractVector`: Vector with the measurements weights for the least-square
    algorithm. We assemble the weight matrix `W` as a diagonal matrix with the elements in
    `weight_vector` at its diagonal.
    (**Default**: `@SVector(ones(Bool, 6))`)

# Returns

- `KeplerianElements{MeanAnomaly, Float64, Float64}`: Fitted Keplerian elements.
- `SMatrix{6, 6, Float64}`: Final covariance matrix of the least-square algorithm.
"""
function Propagators.fit_mean_elements(
    ::Val{:J2osc},
    vjd::AbstractVector{Tjd},
    vr_i::AbstractVector{Tv},
    vv_i::AbstractVector{Tv};
    kwargs...,
) where {Tjd <: Number, Tv <: AbstractVector}
    return fit_j2osc_mean_elements(vjd, vr_i, vv_i; kwargs...)
end

"""
    Propagators.fit_mean_elements!(
        orbp::OrbitPropagatorJ2Osculating,
        vjd::AbstractVector{Tjd},
        vr_i::AbstractVector{Tv},
        vv_i::AbstractVector{Tv};
        kwargs...
    ) where {
        Tjd <: Number,
        Tv <: AbstractVector
    } -> KeplerianElements{MeanAnomaly, Tepoch, T}, SMatrix{6, 6, T}

Fit a set of mean Keplerian elements for the J2 osculating orbit propagator `orbp` using the
osculating elements represented by a set of position vectors `vr_i` [m] and a set of
velocity vectors `vv_i` [m / s] represented in an inertial reference frame at instants in
the array `vjd` [Julian Day].

!!! note

    The orbit propagator `orbp` will be initialized with the Keplerian elements returned
    by the function.

# Keywords

- `atol::Number`: Tolerance for the residual absolute value. If the residual is lower than
    `atol` at any iteration, the computation loop stops.
    (**Default**: 2e-4)
- `rtol::Number`: Tolerance for the relative difference between the residuals. If the
    relative difference between the residuals in two consecutive iterations is lower than
    `rtol`, the computation loop stops.
    (**Default**: 2e-4)
- `initial_guess::Union{Nothing, KeplerianElements}`: Initial guess for the mean elements
    fitting process. If it is `nothing`, the algorithm will obtain an initial estimate from
    the osculating elements in `vr_i` and `vv_i`.
    (**Default**: nothing)
- `jacobian_method::AbstractJacobianMethod`: Method used to compute the Jacobian matrix. It
    can be `FiniteDiffJacobian()` for finite differences or `ForwardDiffJacobian()` for
    `ForwardDiff.jl` automatic differentiation.
    (**Default**: `FiniteDiffJacobian()`)
- `jacobian_perturbation::Number`: Initial state perturbation to compute the
    finite-difference when calculating the Jacobian matrix.
    (**Default**: 1e-3)
- `jacobian_perturbation_tol::Number`: Tolerance to accept the perturbation when calculating
    the Jacobian matrix. If the computed perturbation is lower than
    `jacobian_perturbation_tol`, we increase it until its absolute value is higher than
    `jacobian_perturbation_tol`.
    (**Default**: 1e-7)
- `max_iterations::Int`: Maximum number of iterations allowed for the least-square fitting.
    (**Default**: 50)
- `mean_elements_epoch::Number`: Epoch for the fitted mean elements.
    (**Default**: vjd[end])
- `verbose::Bool`: If `true`, the algorithm prints debugging information to `stdout`.
    (**Default**: true)
- `weight_vector::AbstractVector`: Vector with the measurements weights for the least-square
    algorithm. We assemble the weight matrix `W` as a diagonal matrix with the elements in
    `weight_vector` at its diagonal.
    (**Default**: `@SVector(ones(Bool, 6))`)

# Returns

- `KeplerianElements{MeanAnomaly, Tepoch, T}`: Fitted Keplerian elements.
- `SMatrix{6, 6, T}`: Final covariance matrix of the least-square algorithm.
"""
function Propagators.fit_mean_elements!(
    orbp::OrbitPropagatorJ2Osculating,
    vjd::AbstractVector{Tjd},
    vr_i::AbstractVector{Tv},
    vv_i::AbstractVector{Tv};
    kwargs...,
) where {Tjd <: Number, Tv <: AbstractVector}
    return fit_j2osc_mean_elements!(orbp.j2oscd, vjd, vr_i, vv_i; kwargs...)
end

"""
    Propagators.init(
        Val(:J2osc),
        orb₀::KeplerianElements;
        kwargs...
    ) -> OrbitPropagatorJ2Osculating

Create and initialize the J2 osculating orbit propagator structure using the mean Keplerian
elements `orb₀` [SI units].

!!! note

    The type used in the propagation will be the same as used to define the constants in the
    structure `j2c`.

# Keywords

- `j2c::J2PropagatorConstants`: J2 orbit propagator constants (see
  [`J2PropagatorConstants`](@ref)).
    (**Default**: `J2C_EGM2008`)
"""
function Propagators.init(
    ::Val{:J2osc}, orb₀::KeplerianElements; j2c::J2PropagatorConstants = J2C_EGM2008
)
    j2oscd = j2osc_init(orb₀; j2c = j2c)
    return OrbitPropagatorJ2Osculating(j2oscd)
end

"""
    Propagators.init!(orbp::OrbitPropagatorJ2Osculating, orb₀::KeplerianElements) -> Nothing

Initialize the J2 osculating orbit propagator structure `orbp` using the mean Keplerian
elements `orb₀` [SI units].

!!! warning

    The propagation constants `j2c::J2PropagatorConstants` in `orbp.j2oscd.j2d` will not be
    changed. Hence, they must be initialized.
"""
function Propagators.init!(orbp::OrbitPropagatorJ2Osculating, orb₀::KeplerianElements)
    j2osc_init!(orbp.j2oscd, orb₀)
    return nothing
end

"""
    Propagators.propagate!(
        orbp::OrbitPropagatorJ2Osculating{Tepoch, T},
        t::Number
    ) where {Tepoch <: Number, T <: Number} -> SVector{3, T}, SVector{3, T}

Propagate the orbit of the J2 osculating orbit propagator `orbp` to `t` [s] after the epoch
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
function Propagators.propagate!(orbp::OrbitPropagatorJ2Osculating, t::Number)
    return j2osc!(orbp.j2oscd, t)
end

############################################################################################
#                                        Julia API                                         #
############################################################################################

function Base.copy(
    orbp::OrbitPropagatorJ2Osculating{Tepoch, T}
) where {Tepoch <: Number, T <: Number}
    return OrbitPropagatorJ2Osculating{Tepoch, T}(copy(orbp.j2oscd))
end
