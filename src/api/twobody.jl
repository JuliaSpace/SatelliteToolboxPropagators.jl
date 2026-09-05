## Description #############################################################################
#
#  API implementation for two-body orbit propagator.
#
############################################################################################

# Implement the `Propagators` API for the two-body orbit propagator.
Propagators.epoch(orbp::OrbitPropagatorTwoBody)         = orbp.tbd.orb₀.epoch
Propagators.last_instant(orbp::OrbitPropagatorTwoBody)  = orbp.tbd.Δt
Propagators.mean_elements(orbp::OrbitPropagatorTwoBody) = orbp.tbd.orbk
Propagators.name(orbp::OrbitPropagatorTwoBody)          = "Two-Body Orbit Propagator"

"""
    Propagators.fit_mean_elements(::Val{:TwoBody}, vjd::AbstractVector{Tjd}, vr_i::AbstractVector{Tv}, vv_i::AbstractVector{Tv}; kwargs...) where {Tjd <: Number, Tv <: AbstractVector} -> KeplerianElements{MeanAnomaly, Float64, Float64}, SMatrix{6, 6, Float64}

Fit a set of mean Keplerian elements for the two-body orbit propagator using the osculating
elements represented by a set of position vectors `vr_i` [m] and a set of velocity vectors
`vv_i` [m / s] represented in an inertial reference frame at instants in the array `vjd`
[Julian Day].

!!! note

    This algorithm version will allocate a new two-body propagator with the default
    gravitational parameter `TBC_M0`. If another value is required, use the function
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
    (**Default**: `nothing`)
- `jacobian_method::Union{FiniteDiffJacobian, ForwardDiffJacobian}`: Method used to compute
    the Jacobian matrix. Use `FiniteDiffJacobian()` for finite differences or
    `ForwardDiffJacobian()` for **ForwardDiff.jl** automatic differentiation.
    (**Default**: `FiniteDiffJacobian()`)
- `jacobian_perturbation::Number`: Initial state perturbation to compute the
    finite-difference when calculating the Jacobian matrix. Only used with
    `FiniteDiffJacobian()`.
    (**Default**: 1e-3)
- `jacobian_perturbation_tol::Number`: Tolerance to accept the perturbation when calculating
    the Jacobian matrix. If the computed perturbation is lower than
    `jacobian_perturbation_tol`, we increase it until its absolute value is higher than
    `jacobian_perturbation_tol`. Only used with `FiniteDiffJacobian()`.
    (**Default**: 1e-7)
- `max_iterations::Int`: Maximum number of iterations allowed for the least-square fitting.
    (**Default**: 50)
- `mean_elements_epoch::Number`: Epoch for the fitted mean elements [Julian Day].
    (**Default**: `vjd[end]`)
- `verbose::Bool`: If `true`, the algorithm prints debugging information to `stdout`.
    (**Default**: `true`)
- `weight_vector::AbstractVector`: Vector with the measurements weights for the least-square
    algorithm. We assemble the weight matrix `W` as a diagonal matrix with the elements in
    `weight_vector` at its diagonal.
    (**Default**: `@SVector(ones(Bool, 6))`)

# Returns

- `KeplerianElements{MeanAnomaly, Float64, Float64}`: Fitted Keplerian elements.
- `SMatrix{6, 6, Float64}`: Final covariance matrix of the least-square algorithm.
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
    Propagators.fit_mean_elements!(orbp::OrbitPropagatorTwoBody, vjd::AbstractVector{Tjd}, vr_i::AbstractVector{Tv}, vv_i::AbstractVector{Tv}; kwargs...) where {Tjd <: Number, Tv <: AbstractVector} -> KeplerianElements{MeanAnomaly, Tepoch, T}, SMatrix{6, 6, T}

Fit a set of mean Keplerian elements for the two-body orbit propagator `orbp` using the
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
    (**Default**: `nothing`)
- `jacobian_method::Union{FiniteDiffJacobian, ForwardDiffJacobian}`: Method used to compute
    the Jacobian matrix. Use `FiniteDiffJacobian()` for finite differences or
    `ForwardDiffJacobian()` for **ForwardDiff.jl** automatic differentiation.
    (**Default**: `FiniteDiffJacobian()`)
- `jacobian_perturbation::Number`: Initial state perturbation to compute the
    finite-difference when calculating the Jacobian matrix. Only used with
    `FiniteDiffJacobian()`.
    (**Default**: 1e-3)
- `jacobian_perturbation_tol::Number`: Tolerance to accept the perturbation when calculating
    the Jacobian matrix. If the computed perturbation is lower than
    `jacobian_perturbation_tol`, we increase it until its absolute value is higher than
    `jacobian_perturbation_tol`. Only used with `FiniteDiffJacobian()`.
    (**Default**: 1e-7)
- `max_iterations::Int`: Maximum number of iterations allowed for the least-square fitting.
    (**Default**: 50)
- `mean_elements_epoch::Number`: Epoch for the fitted mean elements [Julian Day].
    (**Default**: `vjd[end]`)
- `verbose::Bool`: If `true`, the algorithm prints debugging information to `stdout`.
    (**Default**: `true`)
- `weight_vector::AbstractVector`: Vector with the measurements weights for the least-square
    algorithm. We assemble the weight matrix `W` as a diagonal matrix with the elements in
    `weight_vector` at its diagonal.
    (**Default**: `@SVector(ones(Bool, 6))`)

# Returns

- `KeplerianElements{MeanAnomaly, Tepoch, T}`: Fitted Keplerian elements.
- `SMatrix{6, 6, T}`: Final covariance matrix of the least-square algorithm.
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
elements `orb₀`.

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
`orb₀`.

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
    Propagators.propagate!(orbp::OrbitPropagatorTwoBody{Tepoch, T}, t::Number) where {Tepoch <: Number, T <: Number} -> SVector{3, T}, SVector{3, T}

Propagate the orbit of the two-body orbit propagator `orbp` to `t` [s] after the epoch of
the initial mean elements, updating the internal state of `orbp`.

# Returns

- `SVector{3, T}`: Position vector [m] represented in the inertial frame at propagation
    instant.
- `SVector{3, T}`: Velocity vector [m / s] represented in the inertial frame at propagation
    instant.

# Remarks

The output is represented in the inertial reference frame of the input elements.
"""
function Propagators.propagate!(orbp::OrbitPropagatorTwoBody, t::Number)
    return twobody!(orbp.tbd, t)
end

############################################################################################
#                                        Julia API                                         #
############################################################################################

function Base.copy(
    orbp::OrbitPropagatorTwoBody{Tepoch, T}
) where {Tepoch <: Number, T <: Number}
    return OrbitPropagatorTwoBody{Tepoch, T}(copy(orbp.tbd))
end
