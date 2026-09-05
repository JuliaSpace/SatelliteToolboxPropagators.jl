## Description #############################################################################
#
# J2 orbit propagator algorithm.
#
# This algorithm propagates the orbit considering the perturbed two-body equations as
# presented in [2, p. 372]. It uses the first-order approximation of Kepler's problem,
# considering the effects of secular gravitational perturbations.
#
## References ##############################################################################
#
# [1] Vallado, D. A (2013). Fundamentals of Astrodynamics and Applications. Microcosm
#     Press, Hawthorn, CA, USA.
#
# [2] Kozai, Y (1959). The Motion of a Close Earth Satellite. The Astronomical Journal,
#     v. 64, no. 1274, pp. 367 -- 377.
#
############################################################################################

export J2C_EGM2008, J2C_EGM1996, J2C_JGM02, J2C_JGM03
export J2C_EGM2008_F32, J2C_EGM1996_F32, J2C_JGM02_F32, J2C_JGM03_F32
export j2_init, j2_init!, j2, j2!
export fit_j2_mean_elements, fit_j2_mean_elements!
export update_j2_mean_elements_epoch, update_j2_mean_elements_epoch!

############################################################################################
#                                           TODO                                           #
############################################################################################
#
# 1. Analyze the reference frame representation of the inputs for this algorithm.
#
#   The SGP4 algorithm expects that the input parameters are represented in the TEME (true
#   equator, mean equinox) reference frame. This J2 orbit propagator model requires that the
#   input parameters are consistent with the gravitational perturbation theory in which the
#   `J2` coefficient was computed. Looking at [1, p. 642], it appears that the perturbations
#   are considering a frame in which the Z-axis is aligned with the CIP (Celestial
#   Intermediate Pole, or the Earth rotation axis). Hence, the J2 parameter is defined based
#   on the PEF. Since no rotations or adaptations are programmed, then the input parameters
#   for this propagator should be represented in any reference frame with a true Equator,
#   because of the symmetry.
#
#   This needs to be further analyzed and confirmed.
#
############################################################################################

############################################################################################
#                                        Constants                                         #
############################################################################################

# These constants were obtained from the GFC files. Remember that:
#
#   J_n = -C_n,0 * √(2n + 1)
#

# EGM-08 gravitational constants.
const J2C_EGM2008 = J2PropagatorConstants(
    6378137.0, √(3.986004415e14 / 6378137.0^3), 0.0010826261738522227
)

const J2C_EGM2008_F32 = J2PropagatorConstants{Float32}(
    6378137.0, √(3.986004415e14 / 6378137.0^3), 0.0010826261738522227
)

# EGM-96 gravitational constants.
const J2C_EGM1996 = J2PropagatorConstants(
    6378136.3, √(3.986004415e14 / 6378136.3^3), 0.0010826266835531513
)

const J2C_EGM1996_F32 = J2PropagatorConstants{Float32}(
    6378136.3, √(3.986004415e14 / 6378136.3^3), 0.0010826266835531513
)

# JGM-02 gravitational constants.
const J2C_JGM02 = J2PropagatorConstants(
    6378136.3, √(3.986004415e14 / 6378136.3^3), 0.0010826269256388149
)

const J2C_JGM02_F32 = J2PropagatorConstants{Float32}(
    6378136.3, √(3.986004415e14 / 6378136.3^3), 0.0010826269256388149
)

# JGM-03 gravitational constants.
const J2C_JGM03 = J2PropagatorConstants(
    6378136.3, √(3.986004415e14 / 6378136.3^3), 0.0010826360229829945
)

const J2C_JGM03_F32 = J2PropagatorConstants{Float32}(
    6378136.3, √(3.986004415e14 / 6378136.3^3), 0.0010826360229829945
)

############################################################################################
#                                        Julia API                                         #
############################################################################################

# Define `copy` for the propagator structure.
let
    fields = fieldnames(J2Propagator)
    expressions = [:(new_j2d.$f = j2d.$f) for f in fields]

    @eval begin
        function Base.copy(
            j2d::J2Propagator{Tepoch, T}
        ) where {Tepoch <: Number, T <: Number}
            new_j2d = J2Propagator{Tepoch, T}()
            $(expressions...)
            return new_j2d
        end
    end
end

############################################################################################
#                                        Functions                                         #
############################################################################################

"""
    j2_init(orb₀::KeplerianElements; kwargs...) -> J2Propagator

Create and initialize the J2 orbit propagator structure using the mean Keplerian elements
`orb₀` [SI units].

!!! note

    The type used in the propagation will be the same as used to define the constants in the
    structure `j2c`.

# Keywords

- `j2c::J2PropagatorConstants`: J2 orbit propagator constants (see
    [`J2PropagatorConstants`](@ref)).
    (**Default** = `J2C_EGM2008`)
"""
function j2_init(
    orb₀::KeplerianElements{Tanomaly, Tepoch, Tkepler};
    j2c::J2PropagatorConstants{T} = J2C_EGM2008,
) where {
    Tanomaly <: AbstractAnomaly, Tepoch <: Number, Tkepler <: AbstractFloat, T <: Number
}
    # Allocate the propagator structure.
    j2d = J2Propagator{Tepoch, T}()

    # Assign the constants, which are used in the initialization.
    j2d.j2c = j2c

    # Initialize the propagator and return.
    j2_init!(j2d, orb₀)

    return j2d
end

function j2_init(
    orb₀::KeplerianElements{Tanomaly, Tepoch, Tkepler};
    j2c::J2PropagatorConstants{Tj2c} = J2C_EGM2008,
) where {Tanomaly <: AbstractAnomaly, Tepoch <: Number, Tkepler <: Number, Tj2c <: Number}
    T = promote_type(Tj2c, Tkepler)

    # Allocate the propagator structure.
    j2d = J2Propagator{Tepoch, T}()

    # Assign the constants, which are used in the initialization.
    j2d.j2c = J2PropagatorConstants{T}(j2c.R0, j2c.μm, j2c.J2)

    # Initialize the propagator and return.
    j2_init!(j2d, orb₀)

    return j2d
end

"""
    j2_init!(j2d::J2Propagator, orb₀::KeplerianElements) -> Nothing

Initialize the J2 orbit propagator structure `j2d` using the mean Keplerian elements `orb₀`
[SI units].

!!! warning

    The propagation constants `j2c::J2PropagatorConstants` in `j2d` will not be changed.
    Hence, they must be initialized.
"""
function j2_init!(
    j2d::J2Propagator{Tepoch, T}, orb₀::KeplerianElements
) where {Tepoch <: Number, T <: Number}
    # Unpack the gravitational constants to improve code readability.
    j2c = j2d.j2c
    R₀  = j2c.R0
    μm  = j2c.μm
    J₂  = j2c.J2

    # Unpack orbit elements.
    a₀ = T(orb₀.semi_major_axis)
    e₀ = T(orb₀.eccentricity)
    i₀ = T(orb₀.inclination)

    # The theory implemented here is only valid for elliptical orbits. Without this check,
    # the user would get a `DomainError` from an internal square root, or silently wrong
    # results, instead of a message pointing at the offending element.
    if !(0 <= e₀ < 1)
        throw(
            ArgumentError("The eccentricity must be in the interval [0, 1), but it is $e₀.")
        )
    end

    if a₀ * (1 - e₀) <= 0
        throw(
            ArgumentError(
                "The perigee radius must be positive, but the semi-major axis is $a₀ m and the eccentricity is $e₀.",
            ),
        )
    end

    # Make sure the Keplerian elements use the mean anomaly. The conversion requires a
    # valid eccentricity, so it must happen after the checks above.
    ke₀ = convert(KeplerianElements{MeanAnomaly, Tepoch, T}, orb₀)

    # Initial values and auxiliary variables.
    al₀ = a₀ / R₀            # ............................. Normalized semi-major axis [er]
    e₀² = e₀^2               # .................................... Eccentricity squared [ ]
    n₀  = μm / √(al₀^3)      # ............................. Unperturbed mean motion [rad/s]
    p₀  = al₀ * (1 - e₀²)    # ...................................... Semi-latus rectum [er]
    p₀² = p₀^2               # ............................. Semi-latus rectum squared [er²]
    M₀  = mean_anomaly(ke₀)  # ................................... Initial mean anomaly [rad]

    sin_i₀, cos_i₀ = sincos(T(i₀))
    sin_i₀² = sin_i₀^2
    β² = 1 - e₀²
    β = √β²

    # We use the algorithm provided in [2, p. 372] that consists of updating the Keplerian
    # elements considering only the first order secular terms, i.e., those that depends only
    # on J₂.

    # We need to compute the perturbed mean motion that is used to calculate the first-order
    # time derivative of the orbital elements [2].
    kn₂ = J₂ / p₀² * β
    n̄ = n₀ * (1 + (3//4) * kn₂ * (2 - 3sin_i₀²))

    # First-order time-derivative of the orbital elements.
    k̄₂ = n̄ * J₂ / p₀²
    ∂Ω = -(3//2) * k̄₂ * cos_i₀
    ∂ω = +(3//4) * k̄₂ * (4 - 5sin_i₀²)

    # Initialize the propagator structure with the data.
    j2d.orb₀ = j2d.orbk = ke₀
    j2d.Δt   = 0
    j2d.∂Ω   = ∂Ω
    j2d.∂ω   = ∂ω
    j2d.n̄    = n̄

    return nothing
end

"""
    j2(Δt::Number, orb₀::KeplerianElements; kwargs...) -> SVector{3, T}, SVector{3, T}, J2Propagator

Initialize the J2 propagator structure using the input elements `orb₀` [SI units] and
propagate the orbit until the time Δt [s].

!!! note

    The type used in the propagation will be the same as used to define the constants in the
    structure `j2c`.

# Keywords

- `j2c::J2PropagatorConstants`: J2 orbit propagator constants (see
    [`J2PropagatorConstants`](@ref)).
    (**Default** = `J2C_EGM2008`)

# Returns

- `SVector{3, T}`: Position vector [m] represented in the inertial frame at propagation
    instant.
- `SVector{3, T}`: Velocity vector [m / s] represented in the inertial frame at propagation
    instant.
- [`J2Propagator`](@ref): Structure with the initialized propagator.

# Remarks

The inertial frame in which the output is represented depends on which frame was used to
generate the orbit parameters. Notice that the perturbation theory requires an inertial
frame with true equator.
"""
function j2(Δt::Number, orb₀::KeplerianElements; j2c::J2PropagatorConstants = J2C_EGM2008)
    j2d = j2_init(orb₀; j2c = j2c)
    r_i, v_i = j2!(j2d, Δt)
    return r_i, v_i, j2d
end

"""
    j2!(j2d::J2Propagator{Tepoch, T}, t::Number) where {Tepoch, T} -> SVector{3, T}, SVector{3, T}

Propagate the orbit defined in `j2d` (see [`J2Propagator`](@ref)) to `t` [s] after the
epoch of the input mean elements in `j2d`.

!!! note

    The internal values in `j2d` will be modified.

# Returns

- `SVector{3, T}`: Position vector [m] represented in the inertial frame at propagation
    instant.
- `SVector{3, T}`: Velocity vector [m / s] represented in the inertial frame at propagation
    instant.

# Remarks

The inertial frame in which the output is represented depends on which frame was used to
generate the orbit parameters. Notice that the perturbation theory requires an inertial
frame with true equator.
"""
function j2!(j2d::J2Propagator{Tepoch, T}, t::Number) where {Tepoch <: Number, T <: Number}
    orbk = _j2_mean_elements!(j2d, t)

    # Compute the position and velocity vectors given the orbital elements.
    r_i_k, v_i_k = kepler_to_rv(orbk)

    # Return the position and velocity vector represented in the inertial reference frame.
    return r_i_k, v_i_k
end

# Propagate the mean elements of `j2d` to the instant `t` [s] measured from the epoch of the
# initial elements, update the propagator structure, and return the mean elements. The
# osculating propagator uses this function directly to avoid computing a state vector from
# the mean elements that it would discard.
function _j2_mean_elements!(
    j2d::J2Propagator{Tepoch, T}, t::Number
) where {Tepoch <: Number, T <: Number}
    # Unpack the variables.
    orb₀  = j2d.orb₀
    ∂Ω    = j2d.∂Ω
    ∂ω    = j2d.∂ω
    n̄     = j2d.n̄
    epoch = orb₀.epoch
    a₀    = orb₀.semi_major_axis
    e₀    = orb₀.eccentricity
    i₀    = orb₀.inclination
    Ω₀    = orb₀.raan
    ω₀    = orb₀.argument_of_periapsis
    M₀    = mean_anomaly(orb₀)

    # Time from epoch to propagate the orbit.
    Δt = T(t)

    # Propagate the orbital elements.
    Ω_k = mod(Ω₀ + ∂Ω * Δt, T(2π))
    ω_k = mod(ω₀ + ∂ω * Δt, T(2π))
    M_k = mod(M₀ + n̄ * Δt, T(2π))

    # Assemble the current mean elements.
    orbk = KeplerianElements{MeanAnomaly}(
        epoch + Tepoch(t) / 86400, a₀, e₀, i₀, Ω_k, ω_k, M_k
    )

    # Update the J2 orbit propagator structure.
    j2d.Δt   = Δt
    j2d.orbk = orbk

    return orbk
end

"""
    fit_j2_mean_elements(vjd::AbstractVector{Tjd}, vr_i::AbstractVector{Tv}, vv_i::AbstractVector{Tv}; kwargs...) where {Tjd<:Number, Tv<:AbstractVector} -> KeplerianElements{TrueAnomaly, Float64, Float64}, SMatrix{6, 6, Float64}

Fit a set of mean Keplerian elements for the J2 orbit propagator using the osculating
elements represented by a set of position vectors `vr_i` [m] and a set of velocity vectors
`vv_i` [m / s] represented in an inertial reference frame at instants in the array `vjd`
[Julian Day].

!!! note

    This algorithm version will allocate a new J2 propagator with the default constants
    `J2C_EGM2008`. If another set of constants are required, use the function
    [`fit_j2_mean_elements!`](@ref) instead.

# Keywords

- `atol::Number`: Tolerance for the residual absolute value. If the residual is lower than
    `atol` at any iteration, the computation loop stops.
    (**Default** = 2e-4)
- `rtol::Number`: Tolerance for the relative difference between the residuals. If the
    relative difference between the residuals in two consecutive iterations is lower than
    `rtol`, the computation loop stops.
    (**Default** = 2e-4)
- `initial_guess::Union{Nothing, KeplerianElements}`: Initial guess for the mean elements
    fitting process. If it is `nothing`, the algorithm will obtain an initial estimate from
    the osculating elements in `vr_i` and `vv_i`.
    (**Default** = nothing)
- `jacobian_method::Union{FiniteDiffJacobian, ForwardDiffJacobian}`: Method used to compute
    the Jacobian matrix. Use `FiniteDiffJacobian()` for finite differences or
    `ForwardDiffJacobian()` for `ForwardDiff.jl` automatic differentiation.
    (**Default** = `FiniteDiffJacobian()`)
- `jacobian_perturbation::Number`: Initial state perturbation to compute the
    finite-difference when calculating the Jacobian matrix. Only used with
    `FiniteDiffJacobian()`.
    (**Default** = 1e-3)
- `jacobian_perturbation_tol::Number`: Tolerance to accept the perturbation when calculating
    the Jacobian matrix. If the computed perturbation is lower than
    `jacobian_perturbation_tol`, we increase it until its absolute value is higher than
    `jacobian_perturbation_tol`. Only used with `FiniteDiffJacobian()`.
    (**Default** = 1e-7)
- `max_iterations::Int`: Maximum number of iterations allowed for the least-square fitting.
    (**Default** = 50)
- `mean_elements_epoch::Number`: Epoch for the fitted mean elements.
    (**Default** = vjd[end])
- `verbose::Bool`: If `true`, the algorithm prints debugging information to `stdout`.
    (**Default** = true)
- `weight_vector::AbstractVector`: Vector with the measurements weights for the least-square
    algorithm. We assemble the weight matrix `W` as a diagonal matrix with the elements in
    `weight_vector` at its diagonal.
    (**Default** = `@SVector(ones(Bool, 6))`)

# Returns

- `KeplerianElements{TrueAnomaly, Float64, Float64}`: Fitted Keplerian elements.
- `SMatrix{6, 6, Float64}`: Final covariance matrix of the least-square algorithm.

# Examples

```julia-repl
julia> vr_i = [
           [-6792.402703741442, 2192.6458461287293, 0.18851758695295118]  .* 1000,
           [-1781.214419290065, 1619.7795321872854, 6707.771633846665]    .* 1000,
           [ 5693.643675547716, -1192.342828671633, 4123.976025977494]    .* 1000,
           [ 5291.613719530499, -2354.5417593130833, -4175.561367156414]  .* 1000,
           [-2416.3705905186903, -268.74923235392623, -6715.411357310478] .* 1000,
           [-6795.043410709359, 2184.4414321930635, -0.4327055325971031]  .* 1000,
       ];

julia> vv_i = [
           [0.3445760107690598, 1.0395135806993514, 7.393686131436984]    .* 1000,
           [6.875680282038698, -1.864319399615942, 2.270603214569518]     .* 1000,
           [3.8964090757666496, -2.1887896252945875, -5.9960180359219075] .* 1000,
           [-4.470258022565413, 0.5119576359985208, -5.9608372367141635]  .* 1000,
           [-6.647358060413909, 2.495415251255861, 2.292118747543002]     .* 1000,
           [0.3427096905434428, 1.040125572862349, 7.3936887585116855]    .* 1000,
       ];

julia> vjd = [
           2.46002818657856e6
           2.460028200467449e6
           2.460028214356338e6
           2.4600282282452267e6
           2.4600282421341157e6
           2.4600282560230047e6
       ];

julia> orb, P = fit_j2_mean_elements(vjd, vr_i, vv_i)
ACTION:   Fitting the mean elements for the J2 propagator.
           Iteration        Position RMSE        Velocity RMSE           Total RMSE       RMSE Variation
                                     [km]             [km / s]                  [ ]
PROGRESS:         23               4.3413           0.00540076              4341.31         -2.90721e-06 %

(KeplerianElements{Float64, Float64}: Epoch = 2.46003e6 (2023-03-24T18:08:40.388), [0.16604846233666615 0.06643574144803302 … -3.85541368066423e-5 0.000124032254172014; 0.0664357414470214 0.26633448262787296 … -1.7943612056596758e-5 -1.9567956793856743e-5; … ; -3.855413680493565e-5 -1.7943612058983507e-5 … 4.3971984339116176e-7 -8.092704691911699e-8; 0.0001240322541726098 -1.9567956793417353e-5 … -8.092704692135158e-8 1.2451922454639337e-7])

julia> orb
KeplerianElements{Float64, Float64}:
           Epoch :    2.46003e6 (2023-03-24T18:08:40.388)
 Semi-major axis : 7131.63       km
    Eccentricity :    0.00114299
     Inclination :   98.4366     °
            RAAN :  162.177      °
 Arg. of Perigee :  101.286      °
    True Anomaly :  258.689      °
```
"""
function fit_j2_mean_elements(
    vjd::AbstractVector{Tjd}, vr_i::AbstractVector{Tv}, vv_i::AbstractVector{Tv}; kwargs...
) where {Tjd <: Number, Tv <: AbstractVector}
    # Allocate the J2 propagator structure that will propagate the mean elements.
    j2d = J2Propagator{Float64, Float64}()

    # Assign the constants, which are used in the initialization.
    j2d.j2c = J2C_EGM2008

    return fit_j2_mean_elements!(j2d, vjd, vr_i, vv_i; kwargs...)
end

"""
    fit_j2_mean_elements!(j2d::J2Propagator{Tepoch, T}, vjd::AbstractVector{Tjd}, vr_i::AbstractVector{Tv}, vv_i::AbstractVector{Tv}; kwargs...) where {T<:Number, Tepoch<:Number, Tjd<:Number, Tv<:AbstractVector} -> KeplerianElements{TrueAnomaly, Tepoch, T}, SMatrix{6, 6, T}

Fit a set of mean Keplerian elements for the J2 orbit propagator `j2d` using the osculating
elements represented by a set of position vectors `vr_i` [m] and a set of velocity vectors
`vv_i` [m / s] represented in an inertial reference frame at instants in the array `vjd`
[Julian Day].

!!! note

    The J2 orbit propagator `j2d` will be initialized with the Keplerian elements returned
    by the function.

# Keywords

- `atol::Number`: Tolerance for the residual absolute value. If the residual is lower than
    `atol` at any iteration, the computation loop stops.
    (**Default** = 2e-4)
- `rtol::Number`: Tolerance for the relative difference between the residuals. If the
    relative difference between the residuals in two consecutive iterations is lower than
    `rtol`, the computation loop stops.
    (**Default** = 2e-4)
- `initial_guess::Union{Nothing, KeplerianElements}`: Initial guess for the mean elements
    fitting process. If it is `nothing`, the algorithm will obtain an initial estimate from
    the osculating elements in `vr_i` and `vv_i`.
    (**Default** = nothing)
- `jacobian_method::Union{FiniteDiffJacobian, ForwardDiffJacobian}`: Method used to compute
    the Jacobian matrix. Use `FiniteDiffJacobian()` for finite differences or
    `ForwardDiffJacobian()` for `ForwardDiff.jl` automatic differentiation.
    (**Default** = `FiniteDiffJacobian()`)
- `jacobian_perturbation::Number`: Initial state perturbation to compute the
    finite-difference when calculating the Jacobian matrix. Only used with
    `FiniteDiffJacobian()`.
    (**Default** = 1e-3)
- `jacobian_perturbation_tol::Number`: Tolerance to accept the perturbation when calculating
    the Jacobian matrix. If the computed perturbation is lower than
    `jacobian_perturbation_tol`, we increase it until its absolute value is higher than
    `jacobian_perturbation_tol`. Only used with `FiniteDiffJacobian()`.
    (**Default** = 1e-7)
- `max_iterations::Int`: Maximum number of iterations allowed for the least-square fitting.
    (**Default** = 50)
- `mean_elements_epoch::Number`: Epoch for the fitted mean elements.
    (**Default** = vjd[end])
- `verbose::Bool`: If `true`, the algorithm prints debugging information to `stdout`.
    (**Default** = true)
- `weight_vector::AbstractVector`: Vector with the measurements weights for the least-square
    algorithm. We assemble the weight matrix `W` as a diagonal matrix with the elements in
    `weight_vector` at its diagonal.
    (**Default** = `@SVector(ones(Bool, 6))`)

# Returns

- `KeplerianElements{TrueAnomaly, Tepoch, T}`: Fitted Keplerian elements.
- `SMatrix{6, 6, T}`: Final covariance matrix of the least-square algorithm.

# Examples

```julia-repl
# Allocate a new J2 orbit propagator using a dummy set of Keplerian elements.
julia> j2d = j2_init(KeplerianElements(0, 7000e3, 0, 0, 0, 0, 0));

julia> vr_i = [
           [-6792.402703741442, 2192.6458461287293, 0.18851758695295118]  .* 1000,
           [-1781.214419290065, 1619.7795321872854, 6707.771633846665]    .* 1000,
           [ 5693.643675547716, -1192.342828671633, 4123.976025977494]    .* 1000,
           [ 5291.613719530499, -2354.5417593130833, -4175.561367156414]  .* 1000,
           [-2416.3705905186903, -268.74923235392623, -6715.411357310478] .* 1000,
           [-6795.043410709359, 2184.4414321930635, -0.4327055325971031]  .* 1000,
       ];

julia> vv_i = [
           [0.3445760107690598, 1.0395135806993514, 7.393686131436984]    .* 1000,
           [6.875680282038698, -1.864319399615942, 2.270603214569518]     .* 1000,
           [3.8964090757666496, -2.1887896252945875, -5.9960180359219075] .* 1000,
           [-4.470258022565413, 0.5119576359985208, -5.9608372367141635]  .* 1000,
           [-6.647358060413909, 2.495415251255861, 2.292118747543002]     .* 1000,
           [0.3427096905434428, 1.040125572862349, 7.3936887585116855]    .* 1000,
       ];

julia> vjd = [
           2.46002818657856e6
           2.460028200467449e6
           2.460028214356338e6
           2.4600282282452267e6
           2.4600282421341157e6
           2.4600282560230047e6
       ];

julia> orb, P = fit_j2_mean_elements!(j2d, vjd, vr_i, vv_i)
ACTION:   Fitting the mean elements for the J2 propagator.
           Iteration        Position RMSE        Velocity RMSE           Total RMSE       RMSE Variation
                                     [km]             [km / s]                  [ ]
PROGRESS:         23               4.3413           0.00540076              4341.31         -2.90721e-06 %

(KeplerianElements{Float64, Float64}: Epoch = 2.46003e6 (2023-03-24T18:08:40.388), [0.16604846233666615 0.06643574144803302 … -3.85541368066423e-5 0.000124032254172014; 0.0664357414470214 0.26633448262787296 … -1.7943612056596758e-5 -1.9567956793856743e-5; … ; -3.855413680493565e-5 -1.7943612058983507e-5 … 4.3971984339116176e-7 -8.092704691911699e-8; 0.0001240322541726098 -1.9567956793417353e-5 … -8.092704692135158e-8 1.2451922454639337e-7])

julia> orb
KeplerianElements{Float64, Float64}:
           Epoch :    2.46003e6 (2023-03-24T18:08:40.388)
 Semi-major axis : 7131.63       km
    Eccentricity :    0.00114299
     Inclination :   98.4366     °
            RAAN :  162.177      °
 Arg. of Perigee :  101.286      °
    True Anomaly :  258.689      °
```
"""
function fit_j2_mean_elements!(
    j2d::J2Propagator{Tepoch, T},
    vjd::AbstractVector{Tjd},
    vr_i::AbstractVector{Tv},
    vv_i::AbstractVector{Tv};
    kwargs...,
) where {Tepoch <: Number, T <: Number, Tjd <: Number, Tv <: AbstractVector}
    return _fit_mean_elements!(j2d, vjd, vr_i, vv_i; kwargs...)
end

"""
    update_j2_mean_elements_epoch(orb::KeplerianElements, new_epoch::Union{Number, DateTime}) -> KeplerianElements

Update the epoch of the mean elements `orb` using a J2 orbit propagator to `new_epoch`,
which can be represented by a Julian Day or a `DateTime`.

!!! note

    This algorithm version will allocate a new J2 propagator with the default constants
    `J2C_EGM2008`. If another set of constants are required, use the function
    [`update_j2_mean_elements_epoch!`](@ref) instead.

# Examples

```julia-repl
julia> orb = KeplerianElements(
           DateTime("2023-01-01") |> datetime2julian,
           7190.982e3,
           0.001111,
           98.405 |> deg2rad,
           90     |> deg2rad,
           200    |> deg2rad,
           45     |> deg2rad
       )
KeplerianElements{Float64, Float64}:
           Epoch :    2.45995e6 (2023-01-01T00:00:00)
 Semi-major axis : 7190.98     km
    Eccentricity :    0.001111
     Inclination :   98.405    °
            RAAN :   90.0      °
 Arg. of Perigee :  200.0      °
    True Anomaly :   45.0      °

julia> update_j2_mean_elements_epoch(orb, DateTime("2023-01-02"))
KeplerianElements{Float64, Float64}:
           Epoch :    2.45995e6 (2023-01-02T00:00:00)
 Semi-major axis : 7190.98     km
    Eccentricity :    0.001111
     Inclination :   98.405    °
            RAAN :   90.9565   °
 Arg. of Perigee :  197.078    °
    True Anomaly :  127.291    °
```
"""
function update_j2_mean_elements_epoch(
    orb::KeplerianElements{Tanomaly, Tepoch, T}, new_epoch::Union{Number, DateTime}
) where {Tanomaly <: AbstractAnomaly, Tepoch <: Number, T <: Number}
    # Allocate the J2 propagator structure that will propagate the mean elements.
    j2d = J2Propagator{Tepoch, T}()

    # Assign the constants, which are used in the initialization.
    j2d.j2c = J2C_EGM2008

    return update_j2_mean_elements_epoch!(j2d, orb, new_epoch)
end

"""
    update_j2_mean_elements_epoch!(j2d::J2Propagator, orb::KeplerianElements, new_epoch::Union{Number, DateTime}) -> KeplerianElements

Update the epoch of the mean elements `orb` using the propagator `j2d` to `new_epoch`, which
can be represented by a Julian Day or a `DateTime`.

!!! note

    The J2 orbit propagator `j2d` will be initialized with the Keplerian elements returned
    by the function.

# Examples

```julia-repl
julia> orb = KeplerianElements(
           DateTime("2023-01-01") |> datetime2julian,
           7190.982e3,
           0.001111,
           98.405 |> deg2rad,
           90     |> deg2rad,
           200    |> deg2rad,
           45     |> deg2rad
       )
KeplerianElements{Float64, Float64}:
           Epoch :    2.45995e6 (2023-01-01T00:00:00)
 Semi-major axis : 7190.98     km
    Eccentricity :    0.001111
     Inclination :   98.405    °
            RAAN :   90.0      °
 Arg. of Perigee :  200.0      °
    True Anomaly :   45.0      °

# Allocate a new J2 orbit propagator using the created Keplerian elements. Notice that any
# set of Keplerian elements can be used here.
julia> j2d = j2_init(orb);

julia> update_j2_mean_elements_epoch!(j2d, orb, DateTime("2023-01-02"))
KeplerianElements{Float64, Float64}:
           Epoch :    2.45995e6 (2023-01-02T00:00:00)
 Semi-major axis : 7190.98     km
    Eccentricity :    0.001111
     Inclination :   98.405    °
            RAAN :   90.9565   °
 Arg. of Perigee :  197.078    °
    True Anomaly :  127.291    °
```
"""
function update_j2_mean_elements_epoch!(
    j2d::J2Propagator, orb::KeplerianElements, new_epoch::DateTime
)
    dt = datetime2julian(new_epoch)
    return update_j2_mean_elements_epoch!(j2d, orb, dt)
end

function update_j2_mean_elements_epoch!(
    j2d::J2Propagator, orb::KeplerianElements, new_epoch::Number
)
    # First, we need to initialize the J2 propagator with the mean elements.
    j2_init!(j2d, orb)

    # Now, we just need to propagate the orbit to the desired instant and obtain the mean
    # elements from the J2 propagator structure inside.
    Δt = (new_epoch - j2d.orb₀.epoch) * 86400
    j2!(j2d, Δt)
    orb = j2d.orbk

    # Finally, we initialize the propagator with the new set of mean elements.
    j2_init!(j2d, orb)

    return orb
end

############################################################################################
#                                    Private Functions                                     #
############################################################################################

# Create a propagator that the finite-difference Jacobian can use as scratch space, so it
# does not clobber the propagator kept initialized by the fitting loop.
function _create_j2_fd_propagator(j2d::J2Propagator{Tepoch, T}) where {Tepoch, T}
    fd = J2Propagator{Tepoch, T}()
    fd.j2c = j2d.j2c
    return fd
end

function _create_j2_ad_propagator(j2d::J2Propagator{Tepoch, T}) where {Tepoch, T}
    tag = ForwardDiff.Tag{Nothing, T}
    D   = ForwardDiff.Dual{tag, T, 6}
    j2c = j2d.j2c

    ad = J2Propagator{Tepoch, D}()
    ad.j2c = J2PropagatorConstants{D}(D(j2c.R0), D(j2c.μm), D(j2c.J2))
    return ad
end

function _j2_jacobian(
    ::FiniteDiffJacobian,
    j2d::J2Propagator{Tepoch, T},
    Δt::Number,
    x₁::SVector{6, T},
    y₁::SVector{6, T};
    perturbation::Number = T(1e-3),
    perturbation_tol::Number = T(1e-7),
    pd_ad::Union{Nothing, J2Propagator} = nothing,
    pd_fd::Union{Nothing, J2Propagator} = nothing,
) where {T <: Number, Tepoch <: Number}
    # The perturbed propagations below overwrite the propagator. Use the scratch propagator
    # when the caller provides one, so it can keep `j2d` initialized across the
    # measurements instead of reinitializing it for each one.
    epoch = j2d.orb₀.epoch
    fd::J2Propagator{Tepoch, T} = isnothing(pd_fd) ? j2d : pd_fd

    J = MMatrix{6, 6, T}(undef)
    x₂ = x₁

    # Convert the perturbation parameters to the element type once. Otherwise, `ϵ` would be
    # assigned values of two different types, and the Jacobian columns would be computed in
    # the promoted type instead of in `T`.
    ϵ₀ = T(perturbation)
    ϵ_tol = T(perturbation_tol)

    @inbounds @views for j in 1:6
        α = x₂[j]
        ϵ = α * ϵ₀

        for _ in 1:5
            abs(ϵ) > ϵ_tol && break
            ϵ *= T(1.4)
        end

        if abs(ϵ) < ϵ_tol
            ϵ = signbit(α) ? -ϵ_tol : ϵ_tol
        end

        α += ϵ
        x₂ = setindex(x₂, α, j)

        orb = rv_to_kepler(x₂[SOneTo(3)], x₂[StaticArrays.SUnitRange(4, 6)], epoch)
        j2_init!(fd, orb)
        r_i, v_i = j2!(fd, Δt)
        y₂ = @SVector [r_i[1], r_i[2], r_i[3], v_i[1], v_i[2], v_i[3]]

        J[:, j] .= (y₂ .- y₁) ./ ϵ
        x₂ = setindex(x₂, x₁[j], j)
    end

    return SMatrix{6, 6, T}(J)
end

function _j2_jacobian(
    ::ForwardDiffJacobian,
    j2d::J2Propagator{Tepoch, T},
    Δt::Number,
    x₁::SVector{6, T},
    y₁::SVector{6, T};
    perturbation::Number = T(1e-3),
    perturbation_tol::Number = T(1e-7),
    pd_ad::Union{Nothing, J2Propagator} = nothing,
    pd_fd::Union{Nothing, J2Propagator} = nothing,
) where {T <: Number, Tepoch <: Number}
    epoch = j2d.orb₀.epoch
    N     = 6
    tag   = ForwardDiff.Tag{Nothing, T}
    D     = ForwardDiff.Dual{tag, T, N}

    # Declaring the type of the local makes the propagator concrete regardless of what
    # the caller provides. Otherwise, the unparameterised type in the keyword leaves
    # this call and the ones below dynamically dispatched.
    ad::J2Propagator{Tepoch, D} = isnothing(pd_ad) ? _create_j2_ad_propagator(j2d) : pd_ad

    seeds  = ntuple(i -> ForwardDiff.Partials(ntuple(j -> T(i == j), Val(N))), Val(N))
    x_dual = SVector{N, D}(ntuple(i -> D(x₁[i], seeds[i]), Val(N)))

    orb = rv_to_kepler(x_dual[SOneTo(3)], x_dual[StaticArrays.SUnitRange(4, 6)], epoch)
    j2_init!(ad, orb)
    r, v   = j2!(ad, Δt)
    y_dual = vcat(r, v)

    return SMatrix{6, N, T}(
        ntuple(k -> ForwardDiff.partials(y_dual[mod1(k, 6)], cld(k, 6)), Val(6 * N))
    )
end
