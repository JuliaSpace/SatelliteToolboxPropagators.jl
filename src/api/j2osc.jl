# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#
# Description
# ==========================================================================================
#
#    API implementation for J2 osculating orbit propagator.
#
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

Propagators.epoch(orbp::OrbitPropagatorJ2Osculating)         = orbp.j2oscd.j2d.orb₀.t
Propagators.last_instant(orbp::OrbitPropagatorJ2Osculating)  = orbp.j2oscd.Δt
Propagators.mean_elements(orbp::OrbitPropagatorJ2Osculating) = orbp.j2oscd.j2d.orbk
Propagators.name(orbp::OrbitPropagatorJ2Osculating)          = "J2 Osculating Orbit Propagator"

"""
    Propagators.fit_mean_elements(::Val{:J2osc}, vjd::AbstractVector{Tjd}, vr_i::AbstractVector{Tv}, vv_i::AbstractVector{Tv}; kwargs...) -> KeplerianElements{Float64, Float64}, SMatrix{6, 6, Float64}

Fit a set of mean Keplerian elements for the J2 osculating orbit propagator using the
osculating elements represented by a set of position vectors `vr_i` [m] and a set of
velocity vectors `vv_i` [m / s] represented in an inertial reference frame at instants in
the array `vjd` [Julian Day].

!!! note
    This algorithm version will allocate a new J2 osculating propagator with the default
    constants `j2c_egm2008`. If another set of constants are required, use the function
    [`Propagators.fit_mean_elements!`](@ref) instead.

# Keywords

- `atol::Number`: Tolerance for the residue absolute value. If the residue is lower than
    `atol` at any iteration, the computation loop stops. (**Default** = 2e-4)
- `rtol::Number`: Tolerance for the relative difference between the residues. If the
    relative difference between the residues in two consecutive iterations is lower than
    `rtol`, the computation loop stops. (**Default** = 2e-4)
- `initial_guess::Union{Nothing, KeplerianElements}`: Initial guess for the mean elements
    fitting process. If it is `nothing`, the algorithm will obtain an initial estimate from
    the osculating elements in `vr_i` and `vv_i`. (**Default** = nothing)
- `jacobian_perturbation::Number`: Initial state perturbation to compute the
    finite-difference when calculating the Jacobian matrix. (**Default** = 1e-3)
- `jacobian_perturbation_tol::Number`: Tolerance to accept the perturbation when calculating
    the Jacobian matrix. If the computed perturbation is lower than
    `jacobian_perturbation_tol`, we increase it until it absolute value is higher than
    `jacobian_perturbation_tol`. (**Default** = 1e-7)
- `max_iterations::Int`: Maximum number of iterations allowed for the least-square fitting.
    (**Default** = 50)
- `mean_elements_epoch::Number`: Epoch for the fitted mean elements.
    (**Default** = vjd[end])
- `verbose::Bool`: If `true`, the algorithm prints debugging information to `stdout`.
    (**Default** = true)
- `weight_vector::AbstractVector`: Vector with the measurements weights for the least-square
    algorithm. We assemble the weight matrix `W` as a diagonal matrix with the elements in
    `weight_vector` at its diagonal. (**Default** = `@SVector(ones(Bool, 6))`)

# Returns

- `KeplerianElements{Float64, Float64}`: Fitted Keplerian elements.
- `SMatrix{6, 6, Float64}`: Final covariance matrix of the least-square algorithm.
"""
function Propagators.fit_mean_elements(
    ::Val{:J2osc},
    vjd::AbstractVector{Tjd},
    vr_i::AbstractVector{Tv},
    vv_i::AbstractVector{Tv};
    kwargs...
) where {Tjd<:Number, Tv<:AbstractVector}
    return fit_j2osc_mean_elements(vjd, vr_i, vv_i; kwargs...)
end

"""
    Propagators.fit_mean_elements!(orbp::OrbitPropagatorJ2Osculating, vjd::AbstractVector{Tjd}, vr_i::AbstractVector{Tv}, vv_i::AbstractVector{Tv}; kwargs...) where {Tjd<:Number, Tv<:AbstractVector}

Fit a set of mean Keplerian elements for the J2 osculating orbit propagator `orbp` using the
osculating elements represented by a set of position vectors `vr_i` [m] and a set of
velocity vectors `vv_i` [m / s] represented in an inertial reference frame at instants in
the array `vjd` [Julian Day].

!!! note
    The orbit propagator `orbp` will be initialized with the Keplerian elements returned
    by the function.

# Keywords

- `atol::Number`: Tolerance for the residue absolute value. If the residue is lower than
    `atol` at any iteration, the computation loop stops. (**Default** = 2e-4)
- `rtol::Number`: Tolerance for the relative difference between the residues. If the
    relative difference between the residues in two consecutive iterations is lower than
    `rtol`, the computation loop stops. (**Default** = 2e-4)
- `initial_guess::Union{Nothing, KeplerianElements}`: Initial guess for the mean elements
    fitting process. If it is `nothing`, the algorithm will obtain an initial estimate from
    the osculating elements in `vr_i` and `vv_i`. (**Default** = nothing)
- `jacobian_perturbation::Number`: Initial state perturbation to compute the
    finite-difference when calculating the Jacobian matrix. (**Default** = 1e-3)
- `jacobian_perturbation_tol::Number`: Tolerance to accept the perturbation when calculating
    the Jacobian matrix. If the computed perturbation is lower than
    `jacobian_perturbation_tol`, we increase it until it absolute value is higher than
    `jacobian_perturbation_tol`. (**Default** = 1e-7)
- `max_iterations::Int`: Maximum number of iterations allowed for the least-square fitting.
    (**Default** = 50)
- `mean_elements_epoch::Number`: Epoch for the fitted mean elements.
    (**Default** = vjd[end])
- `verbose::Bool`: If `true`, the algorithm prints debugging information to `stdout`.
    (**Default** = true)
- `weight_vector::AbstractVector`: Vector with the measurements weights for the least-square
    algorithm. We assemble the weight matrix `W` as a diagonal matrix with the elements in
    `weight_vector` at its diagonal. (**Default** = `@SVector(ones(Bool, 6))`)

# Returns

- `KeplerianElements{Float64, Float64}`: Fitted Keplerian elements.
- `SMatrix{6, 6, Float64}`: Final covariance matrix of the least-square algorithm.
"""
function Propagators.fit_mean_elements!(
    orbp::OrbitPropagatorJ2Osculating,
    vjd::AbstractVector{Tjd},
    vr_i::AbstractVector{Tv},
    vv_i::AbstractVector{Tv};
    kwargs...
) where {Tjd<:Number, Tv<:AbstractVector}
    return fit_j2osc_mean_elements!(orbp.j2oscd, vjd, vr_i, vv_i; kwargs...)
end

"""
    Propagators.init(Val(:J2osc), orb₀::KeplerianElements; kwargs...) -> OrbitPropagatorJ2Osculating

Create and initialize the J2 osculating orbit propagator structure using the mean Keplerian
elements `orb₀` [SI units].

!!! note
    The type used in the propagation will be the same as used to define the constants in the
    structure `j2c`.

# Keywords

- `j2c::J2PropagatorConstants`: J2 orbit propagator constants (see
  [`J2PropagatorConstants`](@ref)). (**Default** = `j2c_egm2008`)
"""
function Propagators.init(
    ::Val{:J2osc},
    orb₀::KeplerianElements;
    j2c::J2PropagatorConstants = j2c_egm2008
)
    j2oscd = j2osc_init(orb₀; j2c = j2c)
    return OrbitPropagatorJ2Osculating(j2oscd)
end

"""
    Propagators.init!(orbp::OrbitPropagatorJ2Osculating, orb₀::KeplerianElements) -> Nothing

Initialize the J2 osculating orbit propagator structure `orbp` using the mean Keplerian
elements `orb₀` [SI units].

!!! warning
    The propagation constants `j2c::J2PropagatorConstants` in `orbp.j2d` will not be
    changed. Hence, they must be initialized.
"""
function Propagators.init!(orbp::OrbitPropagatorJ2Osculating, orb₀::KeplerianElements)
    j2osc_init!(orbp.j2oscd, orb₀)
    return nothing
end

"""
    Propagators.propagate(Val(:J2osc), Δt::Number, orb₀::KeplerianElements; kwargs...) -> SVector{3, T}, SVector{3, T}, OrbitPropagatorJ2Osculating

Initialize the J2 osculating propagator structure using the input elements `orb₀` [SI units]
and propagate the orbit until the time Δt [s].

!!! note
    The type used in the propagation will be the same as used to define the constants in the
    structure `j2c`.

# Keywords

- `j2c::J2PropagatorConstants{T}`: J2 orbit propagator constants (see
  [`J2PropagatorConstants`](@ref)). (**Default** = `j2c_egm2008`)

# Returns

- `SVector{3, T}`: Position vector [m] represented in the inertial frame at propagation
    instant.
- `SVector{3, T}`: Velocity vector [m / s] represented in the inertial frame at propagation
    instant.
- [`OrbitPropagatorJ2Osculating`](@ref): Structure with the initialized parameters.
"""
function Propagators.propagate(
    ::Val{:J2osc},
    Δt::Number,
    orb₀::KeplerianElements;
    j2c::J2PropagatorConstants = j2c_egm2008
)
    r_i, v_i, j2oscd = j2osc(Δt, orb₀; j2c = j2c)
    return r_i, v_i, OrbitPropagatorJ2Osculating(j2oscd)
end

function Propagators.propagate!(orbp::OrbitPropagatorJ2Osculating, t::Number)
    # Auxiliary variables.
    j2oscd = orbp.j2oscd

    # Propagate the orbit.
    return j2osc!(j2oscd, t)
end
