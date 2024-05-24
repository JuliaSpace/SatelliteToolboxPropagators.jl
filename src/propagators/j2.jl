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

export j2c_egm2008, j2c_egm1996, j2c_jgm02, j2c_jgm03
export j2c_egm2008_f32, j2c_egm1996_f32, j2c_jgm02_f32, j2c_jgm03_f32
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
const j2c_egm2008 = J2PropagatorConstants(
    6378137.0,
    √(3.986004415e14 / 6378137.0^3),
    0.0010826261738522227
)

const j2c_egm2008_f32 = J2PropagatorConstants{Float32}(
    6378137.0,
    √(3.986004415e14 / 6378137.0^3),
    0.0010826261738522227
)

# EGM-96 gravitational constants.
const j2c_egm1996 = J2PropagatorConstants(
    6378136.3,
    √(3.986004415e14 / 6378136.3^3),
    0.0010826266835531513
)

const j2c_egm1996_f32 = J2PropagatorConstants{Float32}(
    6378136.3,
    √(3.986004415e14 / 6378136.3^3),
    0.0010826266835531513
)

# JGM-02 gravitational constants.
const j2c_jgm02 = J2PropagatorConstants(
    6378136.3,
    √(3.986004415e14 / 6378136.3^3),
    0.0010826269256388149
)

const j2c_jgm02_f32 = J2PropagatorConstants{Float32}(
    6378136.3,
    √(3.986004415e14 / 6378136.3^3),
    0.0010826269256388149
)

# JGM-03 gravitational constants.
const j2c_jgm03 = J2PropagatorConstants(
    6378136.3,
    √(3.986004415e14 / 6378136.3^3),
    0.0010826360229829945
)

const j2c_jgm03_f32 = J2PropagatorConstants{Float32}(
    6378136.3,
    √(3.986004415e14 / 6378136.3^3),
    0.0010826360229829945
)

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
    (**Default** = `j2c_egm2008`)
"""
function j2_init(
    orb₀::KeplerianElements{Tepoch, Tkepler};
    j2c::J2PropagatorConstants{T} = j2c_egm2008
) where {Tepoch<:Number, Tkepler<:Number, T<:Number}
    # Allocate the propagator structure.
    j2d = J2Propagator{Tepoch, T}()

    # Assign the constants, which are used in the initialization.
    j2d.j2c = j2c

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
    j2d::J2Propagator{Tepoch, T},
    orb₀::KeplerianElements
) where {Tepoch<:Number, T<:Number}
    # Unpack the gravitational constants to improve code readability.
    j2c = j2d.j2c
    R₀  = j2c.R0
    μm  = j2c.μm
    J₂  = j2c.J2

    # Unpack orbit elements.
    a₀ = T(orb₀.a)
    e₀ = T(orb₀.e)
    i₀ = T(orb₀.i)
    f₀ = T(orb₀.f)

    # Initial values and auxiliary variables.
    al₀ = a₀ / R₀                      # ................... Normalized semi-major axis [er]
    e₀² = e₀^2                         # .......................... Eccentricity squared [ ]
    n₀  = μm / √(al₀^3)                # ................... Unperturbed mean motion [rad/s]
    p₀  = al₀ * (1 - e₀²)              # ............................ Semi-latus rectum [er]
    p₀² = p₀^2                         # ................... Semi-latus rectum squared [er²]
    M₀  = true_to_mean_anomaly(e₀, f₀) # ........................ Initial mean anomaly [rad]

    sin_i₀, cos_i₀ = sincos(T(i₀))
    sin_i₀² = sin_i₀^2
    β²      = 1 - e₀²
    β       = √β²

    # We use the algorithm provided in [2, p. 372] that consists on updating the Keplerian
    # elements considering only the first order secular terms, i.e., those that depends only
    # on J₂.

    # We need to compute the perturbed mean motion that is used to calculate the first-order
    # time derivative of the orbital elements [2].
    kn₂ = J₂  / p₀² * β
    n̄   = n₀ * (1 + (3 // 4) * kn₂ * (2 - 3sin_i₀²))

    # First-order time-derivative of the orbital elements.
    k̄₂ = n̄  * J₂  / p₀²
    ∂Ω = -(3 // 2) * k̄₂ * cos_i₀
    ∂ω = +(3 // 4) * k̄₂ * (4 - 5sin_i₀²)

    # Initialize the propagator structure with the data.
    j2d.orb₀ = j2d.orbk = orb₀
    j2d.Δt   = 0
    j2d.M₀   = M₀
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
    (**Default** = `j2c_egm2008`)

# Returns

- `SVector{3, T}`: Position vector [m] represented in the inertial frame at propagation
    instant.
- `SVector{3, T}`: Velocity vector [m / s] represented in the inertial frame at propagation
    instant.
- [`J2Propagator`](@ref): Structure with the initialized parameters.

# Remarks

The inertial frame in which the output is represented depends on which frame it was used to
generate the orbit parameters. Notice that the perturbation theory requires an inertial
frame with true equator.
"""
function j2(Δt::Number, orb₀::KeplerianElements; j2c::J2PropagatorConstants = j2c_egm2008)
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

The inertial frame in which the output is represented depends on which frame it was used to
generate the orbit parameters. Notice that the perturbation theory requires an inertial
frame with true equator.
"""
function j2!(j2d::J2Propagator{Tepoch, T}, t::Number) where {Tepoch<:Number, T<:Number}
    # Unpack the variables.
    orb₀   = j2d.orb₀
    M₀     = j2d.M₀
    ∂Ω     = j2d.∂Ω
    ∂ω     = j2d.∂ω
    n̄      = j2d.n̄
    epoch  = orb₀.t
    a₀     = orb₀.a
    e₀     = orb₀.e
    i₀     = orb₀.i
    Ω₀     = orb₀.Ω
    ω₀     = orb₀.ω

    # Time from epoch to propagate the orbit.
    Δt = T(t)

    # Propagate the orbital elements.
    Ω_k = mod(Ω₀ + ∂Ω * Δt, T(2π))
    ω_k = mod(ω₀ + ∂ω * Δt, T(2π))
    M_k = mod(M₀ + n̄  * Δt, T(2π))

    # Convert the mean anomaly to the true anomaly.
    f_k = mean_to_true_anomaly(e₀, M_k)

    # Assemble the current mean elements.
    orbk = KeplerianElements(epoch + Tepoch(Δt) / 86400, a₀, e₀, i₀, Ω_k, ω_k, f_k)

    # Compute the position and velocity vectors given the orbital elements.
    r_i_k, v_i_k = kepler_to_rv(orbk)

    # Update the J2 orbit propagator structure.
    j2d.Δt   = Δt
    j2d.orbk = orbk

    # Return the position and velocity vector represented in the inertial reference frame.
    return r_i_k, v_i_k
end

"""
    fit_j2_mean_elements(vjd::AbstractVector{Tjd}, vr_i::AbstractVector{Tv}, vv_i::AbstractVector{Tv}; kwargs...) where {Tjd<:Number, Tv<:AbstractVector} -> KeplerianElements{Float64, Float64}, SMatrix{6, 6, Float64}

Fit a set of mean Keplerian elements for the J2 orbit propagator using the osculating
elements represented by a set of position vectors `vr_i` [m] and a set of velocity vectors
`vv_i` [m / s] represented in an inertial reference frame at instants in the array `vjd`
[Julian Day].

!!! note

    This algorithm version will allocate a new J2 propagator with the default constants
    `j2c_egm2008`. If another set of constants are required, use the function
    [`fit_j2_mean_elements!`](@ref) instead.

# Keywords

- `atol::Number`: Tolerance for the residue absolute value. If the residue is lower than
    `atol` at any iteration, the computation loop stops.
    (**Default** = 2e-4)
- `rtol::Number`: Tolerance for the relative difference between the residues. If the
    relative difference between the residues in two consecutive iterations is lower than
    `rtol`, the computation loop stops.
    (**Default** = 2e-4)
- `initial_guess::Union{Nothing, KeplerianElements}`: Initial guess for the mean elements
    fitting process. If it is `nothing`, the algorithm will obtain an initial estimate from
    the osculating elements in `vr_i` and `vv_i`.
    (**Default** = nothing)
- `jacobian_perturbation::Number`: Initial state perturbation to compute the
    finite-difference when calculating the Jacobian matrix.
    (**Default** = 1e-3)
- `jacobian_perturbation_tol::Number`: Tolerance to accept the perturbation when calculating
    the Jacobian matrix. If the computed perturbation is lower than
    `jacobian_perturbation_tol`, we increase it until it absolute value is higher than
    `jacobian_perturbation_tol`.
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

- `KeplerianElements{Float64, Float64}`: Fitted Keplerian elements.
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
    vjd::AbstractVector{Tjd},
    vr_i::AbstractVector{Tv},
    vv_i::AbstractVector{Tv};
    kwargs...
) where {Tjd<:Number, Tv<:AbstractVector}
    # Allocate the J2 propagator structure that will propagate the mean elements.
    j2d = J2Propagator{Float64, Float64}()

    # Assign the constants, which are used in the initialization.
    j2d.j2c = j2c_egm2008

    return fit_j2_mean_elements!(j2d, vjd, vr_i, vv_i; kwargs...)
end

"""
    fit_j2_mean_elements!(j2d::J2Propagator{Tepoch, T}, vjd::AbstractVector{Tjd}, vr_i::AbstractVector{Tv}, vv_i::AbstractVector{Tv}; kwargs...) where {T<:Number, Tepoch<:Number, Tjd<:Number, Tv<:AbstractVector} -> KeplerianElements{Tepoch, T}, SMatrix{6, 6, T}

Fit a set of mean Keplerian elements for the J2 orbit propagator `j2d` using the osculating
elements represented by a set of position vectors `vr_i` [m] and a set of velocity vectors
`vv_i` [m / s] represented in an inertial reference frame at instants in the array `vjd`
[Julian Day].

!!! note

    The J2 orbit propagator `j2d` will be initialized with the Keplerian elements returned
    by the function.

# Keywords

- `atol::Number`: Tolerance for the residue absolute value. If the residue is lower than
    `atol` at any iteration, the computation loop stops.
    (**Default** = 2e-4)
- `rtol::Number`: Tolerance for the relative difference between the residues. If the
    relative difference between the residues in two consecutive iterations is lower than
    `rtol`, the computation loop stops.
    (**Default** = 2e-4)
- `initial_guess::Union{Nothing, KeplerianElements}`: Initial guess for the mean elements
    fitting process. If it is `nothing`, the algorithm will obtain an initial estimate from
    the osculating elements in `vr_i` and `vv_i`.
    (**Default** = nothing)
- `jacobian_perturbation::Number`: Initial state perturbation to compute the
    finite-difference when calculating the Jacobian matrix.
    (**Default** = 1e-3)
- `jacobian_perturbation_tol::Number`: Tolerance to accept the perturbation when calculating
    the Jacobian matrix. If the computed perturbation is lower than
    `jacobian_perturbation_tol`, we increase it until it absolute value is higher than
    `jacobian_perturbation_tol`.
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

- `KeplerianElements{Tepoch, T}`: Fitted Keplerian elements.
- `SMatrix{6, 6, T}`: Final covariance matrix of the least-square algorithm.

# Examples

```julia-repl
# Allocate a new J2 orbit propagator using a dummy Keplerian elements.
julia> j2d = j2_init(KeplerianElements{Float64, Float64}(0, 7000e3, 0, 0, 0, 0, 0));

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
    atol::Number                                     = 2e-4,
    rtol::Number                                     = 2e-4,
    initial_guess::Union{Nothing, KeplerianElements} = nothing,
    jacobian_perturbation::Number                    = 1e-3,
    jacobian_perturbation_tol::Number                = 1e-7,
    max_iterations::Int                              = 50,
    mean_elements_epoch::Number                      = vjd[end],
    verbose::Bool                                    = true,
    weight_vector::AbstractVector                    = @SVector(ones(Bool, 6)),
) where {T<:Number, Tepoch<:Number, Tjd<:Number, Tv<:AbstractVector}
    # Number of available measurements.
    num_measurements = length(vjd)

    # Check the inputs.
    length(vr_i) != num_measurements &&
        throw(ArgumentError("The number of elements in `vjd` and `vr_i` must be the same."))

    length(vv_i) != num_measurements &&
        throw(ArgumentError("The number of elements in `vjd` and `vv_i` must be the same."))

    if length(weight_vector) != 6
        throw(ArgumentError("The weight vector must have 6 elements."))
    end

    # Check if stdout supports colors.
    has_color = get(stdout, :color, false)::Bool
    cd = has_color ? _D : ""
    cb = has_color ? _B : ""
    cy = has_color ? _Y : ""

    # Assemble the weight matrix.
    W = Diagonal(
        @SVector T[
            weight_vector[1],
            weight_vector[2],
            weight_vector[3],
            weight_vector[4],
            weight_vector[5],
            weight_vector[6]
        ]
    )

    # Initial guess of the mean elements.
    #
    # NOTE: x₁ is the previous estimate and x₂ is the current estimate.
    if !isnothing(initial_guess)
        epoch = T(mean_elements_epoch)

        # First, we need to update the mean elements to the desired epoch.
        verbose && println("$(cy)ACTION:$(cd)   Updating the epoch of the initial mean elements guess to match the desired one.")
        orb = update_j2_mean_elements_epoch!(j2d, initial_guess, epoch)

        r_i, v_i = kepler_to_rv(orb)
        x₁ = SVector{6, T}(r_i[1], r_i[2], r_i[3], v_i[1], v_i[2], v_i[3])
    else
        # In this case, we must find the closest osculating vector to the desired epoch.
        id = firstindex(vjd)
        v  = abs(vjd[1] - mean_elements_epoch)

        for k in eachindex(vjd)
            vk = abs(vjd[k] - mean_elements_epoch)
            if vk < v
                id = k
                v  = vk
            end
        end

        epoch = T(vjd[id])
        r_i   = vr_i[id]
        v_i   = vv_i[id]
        x₁    = SVector{6, T}(r_i[1], r_i[2], r_i[3], v_i[1], v_i[2], v_i[3],)
    end

    x₂ = x₁

    # Number of states in the input vector.
    num_states = 6

    # Covariance matrix.
    P = SMatrix{num_states, num_states, T}(I)

    # Variable to store the last residue.
    local σ_i_₁

    # Variable to store how many iterations the residue increased. This is used to account
    # for divergence.
    Δd = 0

    # Header.
    if verbose
        println("$(cy)ACTION:$(cd)   Fitting the mean elements for the J2 propagator.")
        @printf("          %s%10s %20s %20s %20s %20s%s\n", cy, "Iteration", "Position RMSE", "Velocity RMSE", "Total RMSE", "RMSE Variation", cd)
        @printf("          %s%10s %20s %20s %20s %20s%s\n", cb, "", "[km]", "[km / s]", "[ ]", "", cd)
        println()
    end

    # We need a reference to the covariance inverse because we will invert it and return
    # after the iterations.
    local ΣJ′WJ

    # Loop until the maximum allowed iteration.
    @inbounds @views for it in 1:max_iterations
        x₁ = x₂

        # Variables to store the summations to compute the least square fitting algorithm.
        ΣJ′WJ = @SMatrix zeros(T, num_states, num_states)
        ΣJ′Wb = @SVector zeros(T, num_states)

        # Variable to store the RMS errors in this iteration.
        σ_i  = T(0)
        σp_i = T(0)
        σv_i = T(0)

        for k in 1:num_measurements
            # Obtain the measured ephemerides.
            y = vcat(vr_i[k - 1 + begin], vv_i[k - 1 + begin])

            # Initialize the SGP4 with the current estimated mean elements.
            orb = rv_to_kepler(x₁[1:3], x₁[4:6], epoch)
            j2_init!(j2d, orb)

            # Obtain the propagation time for this measurement.
            Δt = (vjd[k - 1 + begin] - epoch) * 86400

            # Propagate the orbit.
            r̂_i, v̂_i = j2!(j2d, Δt)
            ŷ = vcat(r̂_i, v̂_i)

            # Compute the residue.
            b = y - ŷ

            # Compute the Jacobian in-place.
            J = _j2_jacobian(
                j2d,
                Δt,
                x₁,
                ŷ;
                perturbation     = jacobian_perturbation,
                perturbation_tol = jacobian_perturbation_tol
            )

            # Accumulation.
            ΣJ′WJ += J' * W * J
            ΣJ′Wb += J' * W * b
            σ_i   += b' * W * b
            σp_i  += dot(b[1:3], b[1:3])
            σv_i  += dot(b[4:6], b[4:6])
        end

        # Normalize and compute the RMS errors.
        σ_i  = √(σ_i  / num_measurements)
        σp_i = √(σp_i / num_measurements)
        σv_i = √(σv_i / num_measurements)

        # Update the estimate.
        δx = ΣJ′WJ \ ΣJ′Wb

        # Limit the correction to avoid divergence.
        for i in 1:num_states
            threshold = T(0.1)
            if abs(δx[i] / x₁[i]) > threshold
                δx = setindex(δx, threshold * abs(x₁[i]) * sign(δx[i]), i)
            end
        end

        x₂ = x₁ + δx

        # We cannot compute the RMSE variation in the first iteration.
        if it == 1
            verbose &&
                @printf("\x1b[A\x1b[2K\r%sPROGRESS:%s %10d %20g %20g %20g %20s\n", cb, cd, it, σp_i / 1000, σv_i / 1000, σ_i, "---")

        else
            # Compute the RMSE variation.
            Δσ = (σ_i - σ_i_₁) / σ_i_₁

            verbose &&
                @printf("\x1b[A\x1b[2K\r%sPROGRESS:%s %10d %20g %20g %20g %20g %%\n", cb, cd, it, σp_i / 1000, σv_i / 1000, σ_i, 100 * Δσ)

            # Check if the RMSE is increasing.
            if σ_i < σ_i_₁
                Δd = 0
            else
                Δd += 1
            end

            # If the RMSE increased by three iterations and its value is higher than 5e11,
            # we abort because the iterations are diverging.
            ((Δd ≥ 3) && (σ_i > 5e11)) && error("The iterations diverged!")

            # Check if the condition to stop has been reached.
            ((abs(Δσ) < rtol) || (σ_i < atol) || (it ≥ max_iterations)) && break
        end

        σ_i_₁ = σ_i
    end

    verbose && println()

    # Obtain the mean elements.
    orb = @views rv_to_kepler(x₂[1:3], x₂[4:6], epoch)

    # Update the epoch of the fitted mean elements to match the desired one.
    if abs(epoch - mean_elements_epoch) > 0.001 / 86400
        verbose &&
            println("$(cy)ACTION:$(cd)   Updating the epoch of the fitted mean elements to match the desired one.")
        orb = update_j2_mean_elements_epoch!(j2d, orb, mean_elements_epoch)
    end

    # Initialize the propagator with the mean elements.
    j2_init!(j2d, orb)

    # Compute the final covariance.
    P = pinv(ΣJ′WJ)

    # Return the mean elements and the covariance.
    return orb, P
end

"""
    update_j2_mean_elements_epoch(orb::KeplerianElements, new_epoch::Union{Number, DateTime}) -> KepleriranElements

Update the epoch of the mean elements `orb` using a J2 orbit propagator to `new_epoch`,
which can be represented by a Julian Day or a `DateTime`.

!!! note

    This algorithm version will allocate a new J2 propagator with the default constants
    `j2c_egm2008`. If another set of constants are required, use the function
    [`update_j2osc_mean_elements_epoch!`](@ref) instead.

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
    orb::KeplerianElements{Tepoch, T},
    new_epoch::Union{Number, DateTime}
) where {T<:Number, Tepoch<:Number}
    # Allocate the J2 propagator structure that will propagate the mean elements.
    j2d = J2Propagator{Tepoch, T}()

    # Assign the constants, which are used in the initialization.
    j2d.j2c = j2c_egm2008

    return update_j2_mean_elements_epoch!(j2d, orb, new_epoch)
end

"""
    update_j2_mean_elements_epoch!(j2d::J2Propagator, orb::KeplerianElements, new_epoch::Union{Number, DateTime}) -> KepleriranElements

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
    j2d::J2Propagator,
    orb::KeplerianElements,
    new_epoch::DateTime
)
    dt = datetime2julian(new_epoch)
    return update_j2_mean_elements_epoch!(j2d, orb, dt)
end

function update_j2_mean_elements_epoch!(
    j2d::J2Propagator,
    orb::KeplerianElements,
    new_epoch::Number
)
    # First, we need to initialize the J2 propagator with the mean elements.
    j2_init!(j2d, orb)

    # Now, we just need to propagate the orbit to the desired instant and obtain the mean
    # elements from the J2 propagator structure inside.
    Δt = (new_epoch - j2d.orb₀.t) * 86400
    j2!(j2d, Δt)
    orb = j2d.orbk

    # Finally, we initialize the propagator with the new set of mean elements.
    j2_init!(j2d, orb)

    return orb
end

############################################################################################
#                                    Private Functions                                     #
############################################################################################

"""
    _j2_jacobian(j2d::J2OsculatingPropagator{Tepoch, T}, Δt::Number, x₁::SVector{6, T}, y₁::SVector{6, T}; kwargs...)) where {T<:Number, Tepoch<:Number} -> SMatrix{6, 6, T}

Compute the J2 orbit propagator Jacobian by finite-differences using the propagator `j2d`
at instant `Δt` considering the input mean elements `x₁` that must provide the output vector
`y₁`. Hence:

        ∂j2(x, Δt) │
    J = ────────── │
            ∂x     │ x = x₁

# Keywords

- `perturbation::T`: Initial state perturbation to compute the finite-difference:
    `Δx = x * perturbation`.
    (**Default** = 1e-3)
- `perturbation_tol::T`: Tolerance to accept the perturbation. If the computed perturbation
    is lower than `perturbation_tol`, we increase it until it absolute value is higher than
    `perturbation_tol`.
    (**Default** = 1e-7)
"""
function _j2_jacobian(
    j2d::J2Propagator{Tepoch, T},
    Δt::Number,
    x₁::SVector{6, T},
    y₁::SVector{6, T};
    perturbation::Number = T(1e-3),
    perturbation_tol::Number = T(1e-7)
) where {T<:Number, Tepoch<:Number}

    # Allocate the `MMatrix` that will have the Jacobian.
    J = MMatrix{6, 6, T}(undef)

    # Auxiliary variables.
    x₂ = x₁

    @inbounds @views for j in 1:6
        # State that will be perturbed.
        α = x₂[j]

        # Obtain the perturbation, taking care to avoid small values.
        ϵ = α * T(perturbation)

        for _ in 1:5
            abs(ϵ) > perturbation_tol && break
            ϵ *= T(1.4)
        end

        # Avoid division by zero in cases that α is very small. In this situation, we force
        # `|ϵ| = perturbation_tol`.
        if abs(ϵ) < perturbation_tol
            ϵ = signbit(α) ? -perturbation_tol : perturbation_tol
        end

        α += ϵ

        # Modify the perturbed state.
        x₂ = setindex(x₂, α, j)

        # Obtain the Jacobian by finite differentiation.
        orb = rv_to_kepler(x₂[1:3], x₂[4:6], j2d.orb₀.t)
        j2_init!(j2d, orb)
        r_i, v_i = j2!(j2d, Δt)
        y₂ = @SVector [
            r_i[1],
            r_i[2],
            r_i[3],
            v_i[1],
            v_i[2],
            v_i[3],
        ]

        J[:, j] .= (y₂ .- y₁) ./ ϵ

        # Restore the value of the perturbed state for the next iteration.
        x₂ = setindex(x₂, x₁[j], j)
    end

    # Convert the Jacobian to a static matrix to avoid allocations.
    Js = SMatrix{6, 6, T}(J)

    return Js
end
