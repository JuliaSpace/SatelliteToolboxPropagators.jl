# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#
# Description
# ==========================================================================================
#
#   J4 orbit propagator algorithm.
#
#   This algorithm propagates the orbit considering the secular perturbations of central
#   body zonal harmonics as presented in [1, p. 647-654, 692-692] and [2], which is Kozai's
#   method but neglecting long-periodic and short-periodic perturbations.
#
#   The terms J2, J2², and J4 are considered, i.e. J6 is assumed to be 0. This can be used
#   as a propagator of mean elements for mission analysis in which the satellite orbit is
#   maintained.
#
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#
# References
# ==========================================================================================
#
#   [1] Vallado, D. A (2013). Fundamentals of Astrodynamics and Applications. Microcosm
#       Press, Hawthorn, CA, USA.
#
#   [2] Kozai, Y (1959). The Motion of a Close Earth Satellite. The Astronomical Journal,
#       v. 64, no. 1274, pp. 367 -- 377.
#
#   [3] Hoots, F. R., Roehrich, R. L (1980). Models for Propagation of NORAD Elements Set.
#       Spacetrack Report No. 3.
#
#   [4] Blitzer, L. Handbook of Orbital Perturbations. Astronautics 453. University of
#       Arizona.
#
#   [5] https://www.mathworks.com/matlabcentral/fileexchange/43333-sun-synchronous-orbit-design
#
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

export j4c_egm2008, j4c_egm1996, j4c_jgm02, j4c_jgm03
export j4c_egm2008_f32, j4c_egm1996_f32, j4c_jgm02_f32, j4c_jgm03_f32
export j4_init, j4_init!, j4, j4!
export fit_j4_mean_elements, fit_j4_mean_elements!
export update_j4_mean_elements_epoch, update_j4_mean_elements_epoch!

############################################################################################
#                                        Constants
############################################################################################

# These constants were obtained from the GFC files. Remember that:
#
#   J_n = -C_n,0 * √(2n + 1)
#

# EGM-08 gravitational constants.
const j4c_egm2008 = J4PropagatorConstants(
    6378137.0,
    sqrt(3.986004415e14 / 6378137.0^3),
    0.0010826261738522227,
    -1.6198975999169731e-6
)

const j4c_egm2008_f32 = J4PropagatorConstants{Float32}(
    6378137.0,
    sqrt(3.986004415e14 / 6378137.0^3),
    0.0010826261738522227,
    -1.6198975999169731e-6
)

# EGM-96 gravitational constants.
const j4c_egm1996 = J4PropagatorConstants(
    6378136.3,
    sqrt(3.986004415e14 / 6378136.3^3),
    0.0010826266835531513,
    -1.619621591367e-6
)

const j4c_egm1996_f32 = J4PropagatorConstants{Float32}(
    6378136.3,
    sqrt(3.986004415e14 / 6378136.3^3),
    0.0010826266835531513,
    -1.619621591367e-6
)

# JGM-02 gravitational constants.
const j4c_jgm02 = J4PropagatorConstants(
    6378136.3,
    sqrt(3.986004415e14 / 6378136.3^3),
    0.0010826269256388149,
    -1.62042999e-6
)

const j4c_jgm02_f32 = J4PropagatorConstants{Float32}(
    6378136.3,
    sqrt(3.986004415e14 / 6378136.3^3),
    0.0010826269256388149,
    -1.62042999e-6
)

# JGM-03 gravitational constants.
const j4c_jgm03 = J4PropagatorConstants(
    6378136.3,
    sqrt(3.986004415e14 / 6378136.3^3),
    0.0010826360229829945,
    -1.619331205071e-6
)

const j4c_jgm03_f32 = J4PropagatorConstants{Float32}(
    6378136.3,
    sqrt(3.986004415e14 / 6378136.3^3),
    0.0010826360229829945,
    -1.619331205071e-6
)

############################################################################################
#                                        Functions
############################################################################################

"""
    j4_init(orb₀::KeplerianElements; kwargs...) -> J4Propagator

Create and initialize the J4 orbit propagator structure using the mean Keplerian elements
`orb₀`.

!!! note
    The type used in the propagation will be the same as used to define the constants in the
    structure `j4c`.

# Keywords

- `j4c::J4PropagatorConstants`: J4 orbit propagator constants (see
  [`J4PropagatorConstants`](@ref)). (**Default** = `j4c_egm2008`)
"""
function j4_init(
    orb₀::KeplerianElements{Tepoch, Tkepler};
    j4c::J4PropagatorConstants{T} = j4c_egm2008
) where {Tepoch<:Number, Tkepler<:Number, T<:Number}
    # Allocate the propagator structure.
    j4d = J4Propagator{Tepoch, T}()

    # Assign the constants, which are used in the initialization.
    j4d.j4c = j4c

    # Initialize the propagator and return.
    j4_init!(j4d, orb₀)

    return j4d
end

"""
    j4_init!(j4d::J4Propagator, orb₀::KeplerianElements) -> Nothing

Initialize the J4 orbit propagator structure `j4d` using the mean Keplerian elements `orb₀`.

!!! warning
    The propagation constants `j4c::J4PropagatorConstants` in `j4d` will not be changed.
    Hence, they must be initialized.
"""
function j4_init!(
    j4d::J4Propagator{Tepoch, T},
    orb₀::KeplerianElements,
) where {Tepoch<:Number, T<:Number}
    # Unpack the gravitational constants to improve code readability.
    j4c = j4d.j4c
    R₀  = j4c.R0
    μm  = j4c.μm
    J₂  = j4c.J2
    J₄  = j4c.J4

    # Unpack orbit elements.
    a₀ = T(orb₀.a)
    e₀ = T(orb₀.e)
    i₀ = T(orb₀.i)
    f₀ = T(orb₀.f)

    # Initial values and auxiliary variables.
    al₀ = a₀ / R₀                      # ................... Normalized semi-major axis [er]
    e₀² = e₀^2                         # .......................... Eccentricity squared [ ]
    p₀  = al₀ * (1 - e₀²)              # ............................ Semi-latus rectum [er]
    p₀² = p₀^2                         # ................... Semi-latus rectum squared [er²]
    p₀⁴ = p₀^4                         # .......... Semi-latus rectum to the 4th power [er⁴]
    n₀  = μm / √(al₀^3)                # ................. Unperturbed mean motion [rad / s]
    M₀  = true_to_mean_anomaly(e₀, f₀) # ........................ Initial mean anomaly [rad]
    J₂² = J₂^2                         # ............................... J2 constant squared

    sin_i₀, cos_i₀ = sincos(T(i₀))

    sin_i₀² = sin_i₀^2
    sin_i₀⁴ = sin_i₀^4
    cos_i₀⁴ = cos_i₀^4
    β²      = (1 - e₀²)
    β       = √β²

    # We need to compute the perturbed mean motion that is used to calculate the
    # time-derivative of the orbital elements. This expression was obtained from [5] because
    # [2] does not show it completely.
    kn₂  = J₂  / p₀² * β
    kn₂₂ = J₂² / p₀⁴ * β
    kn₄  = J₄  / p₀⁴ * β

    n̄ = n₀ * (
        1 +
        ( 3 // 4  ) * kn₂  * (2 - 3sin_i₀²) +
        ( 3 // 128) * kn₂₂ * (120 + 64β - 40β² + (-240 - 192β + 40β²) * sin_i₀² + (105 + 144β + 25β²) * sin_i₀⁴) -
        (45 // 128) * kn₄  * e₀² * (-8 + 40sin_i₀² - 35sin_i₀⁴)
    )

    # Some auxiliary variables to compute the perturbations.
    k̄₂  = n̄  * J₂  / p₀²
    k̄₂₂ = n̄  * J₂² / p₀⁴
    k₂₂ = n₀ * J₂² / p₀⁴
    k₄  = n₀ * J₄  / p₀⁴

    # TODO: Check J₄ perturbation term sign.
    #
    # We needed to flip the J₄ perturbation term sign from the value in [1] and [2] to
    # obtain values that match those of STK. However, this modification does not seem right
    # if we observe the RAAN secular perturbation term in SGP4 orbit propagator [3, p. 16].
    # Reference [4] also provides the equation with the sign as in [2]. Furthermore, the
    # current version, which does not match STK's, provides lower errors when comparing to a
    # numerical propagator.
    #
    # For more information, see:
    #
    #   https://github.com/JuliaSpace/SatelliteToolbox.jl/issues/91
    #

    δΩ = -( 3 // 2 ) * k̄₂  * cos_i₀ +
          ( 3 // 32) * k̄₂₂ * cos_i₀ * (-36 -  4e₀² + 48β + (40 - 5e₀² - 72β) * sin_i₀²) +
          (15 // 32) * k₄  * cos_i₀ * (8 + 12e₀² - (14 + 21e₀²) * sin_i₀²)

    δω = ( 3 // 4  ) * k̄₂  * (4 - 5sin_i₀²) +
         ( 3 // 128) * k̄₂₂ * (384 + 96e₀² - 384β + (-824 - 116e₀² + 1056β) * sin_i₀² + (430 - 5e₀² - 720β) * sin_i₀⁴) -
         (15 // 16 ) * k₂₂ * e₀² * cos_i₀⁴ -
         (15 // 128) * k₄  * (64 + 72e₀² - (248 + 252e₀²) * sin_i₀² + (196 + 189e₀²) * sin_i₀⁴)

    # Initialize the propagator structure with the data.
    j4d.orb₀ = j4d.orbk = orb₀
    j4d.Δt   = 0
    j4d.M₀   = M₀
    j4d.δΩ   = δΩ
    j4d.δω   = δω
    j4d.n̄    = n̄

    return nothing
end

"""
    j4(Δt::Number, orb₀::KeplerianElements; kwargs...) -> SVector{3, T}, SVector{3, T}, J4Propagator

Initialize the J4 propagator structure using the input elements `orb₀` and propagate the
orbit until the time Δt [s].

!!! note
    The type used in the propagation will be the same as used to define the constants in the
    structure `j4c`.

# Keywords

- `j4c::J4PropagatorConstants`: J4 orbit propagator constants (see
  [`J4PropagatorConstants`](@ref)). (**Default** = `j4c_egm2008`)

# Returns

- `SVector{3, T}`: Position vector [m] represented in the inertial frame at propagation
    instant.
- `SVector{3, T}`: Velocity vector [m / s] represented in the inertial frame at propagation
    instant.
- [`J4Propagator`](@ref): Structure with the initialized parameters.

# Remarks

The inertial frame in which the output is represented depends on which frame it was used to
generate the orbit parameters. Notice that the perturbation theory requires an inertial
frame with true equator.
"""
function j4(Δt::Number, orb₀::KeplerianElements; j4c::J4PropagatorConstants = j4c_egm2008)
    j4d = j4_init(orb₀; j4c = j4c)
    r_i, v_i = j4!(j4d, Δt)
    return r_i, v_i, j4d
end

"""
    j4!(j4d::J4Propagator{Tepoch, T}, t::Number) where {Tepoch<:Number, T<:Number} -> SVector{3, T}, SVector{3, T}

Propagate the orbit defined in `j4d` (see [`J4Propagator`](@ref)) to `t` [s] after the
epoch of the input mean elements in `j4d`.

!!! note
    The internal values in `j4d` will be modified.

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
function j4!(j4d::J4Propagator{Tepoch, T}, t::Number) where {Tepoch<:Number, T<:Number}
    # Unpack the variables.
    orb₀   = j4d.orb₀
    M₀     = j4d.M₀
    δΩ     = j4d.δΩ
    δω     = j4d.δω
    n̄      = j4d.n̄
    epoch  = orb₀.t
    a₀     = orb₀.a
    e₀     = orb₀.e
    i₀     = orb₀.i
    Ω₀     = orb₀.Ω
    ω₀     = orb₀.ω

    # Time elapsed since epoch.
    Δt = T(t)

    # Propagate the orbital elements.
    Ω_k = mod(Ω₀ + δΩ * Δt, T(2π))
    ω_k = mod(ω₀ + δω * Δt, T(2π))
    M_k = mod(M₀ + n̄  * Δt, T(2π))

    # Convert the mean anomaly to the true anomaly.
    f_k = mean_to_true_anomaly(e₀, M_k)

    # Assemble the current mean elements.
    orbk = KeplerianElements(epoch + Tepoch(Δt) / 86400, a₀, e₀, i₀, Ω_k, ω_k, f_k)

    # Compute the position and velocity vectors given the orbital elements.
    r_i_k, v_i_k = kepler_to_rv(orbk)

    # Update the J2 orbit propagator structure.
    j4d.Δt   = Δt
    j4d.orbk = orbk

    # Return the position and velocity vector represented in the inertial reference frame.
    return r_i_k, v_i_k
end

"""
    fit_j4_mean_elements(vjd::AbstractVector{Tjd}, vr_i::AbstractVector{Tv}, vv_i::AbstractVector{Tv}; kwargs...) where {Tjd<:Number, Tv<:AbstractVector} -> KeplerianElements{Float64, Float64}, SMatrix{6, 6, Float64}

Fit a set of mean Keplerian elements for the J4 orbit propagator using the osculating
elements represented by a set of position vectors `vr_i` [m] and a set of velocity vectors
`vv_i` [m / s] represented in an inertial reference frame at instants in the array `vjd`
[Julian Day].

!!! note
    This algorithm version will allocate a new J4 propagator with the default constants
    `j4c_egm2008`. If another set of constants are required, use the function
    [`fit_j4_mean_elements!`](@ref) instead.

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

julia> orb, P = fit_j4_mean_elements(vjd, vr_i, vv_i)
ACTION:   Fitting the mean elements for the J4 propagator.
           Iteration        Position RMSE        Velocity RMSE           Total RMSE       RMSE Variation
                                     [km]             [km / s]                  [ ]
PROGRESS:         23              4.33863           0.00539962              4338.63         -2.90777e-06 %

(KeplerianElements{Float64, Float64}: Epoch = 2.46003e6 (2023-03-24T18:08:40.388), [0.16604866262321472 0.06643593039980161 … -3.8553206984036995e-5 0.00012403204426258087; 0.06643593040175688 0.26633435614867085 … -1.7942563356579497e-5 -1.95681107792859e-5; … ; -3.855320698470177e-5 -1.794256335557519e-5 … 4.3972013198895673e-7 -8.092682623169118e-8; 0.000124032044261828 -1.9568110780915995e-5 … -8.092682623078798e-8 1.2451901466558624e-7])

julia> orb
KeplerianElements{Float64, Float64}:
           Epoch :    2.46003e6 (2023-03-24T18:08:40.388)
 Semi-major axis : 7131.64       km
    Eccentricity :    0.00114298
     Inclination :   98.4366     °
            RAAN :  162.177      °
 Arg. of Perigee :  101.282      °
    True Anomaly :  258.693      °
```
"""
function fit_j4_mean_elements(
    vjd::AbstractVector{Tjd},
    vr_i::AbstractVector{Tv},
    vv_i::AbstractVector{Tv};
    kwargs...
) where {Tjd<:Number, Tv<:AbstractVector}
    # Allocate the J4 propagator structure that will propagate the mean elements.
    j4d = J4Propagator{Float64, Float64}()

    # Assign the constants, which are used in the initialization.
    j4d.j4c = j4c_egm2008

    return fit_j4_mean_elements!(j4d, vjd, vr_i, vv_i; kwargs...)
end

"""
    fit_j4_mean_elements!(j4d::J4Propagator{Tepoch, T}, vjd::AbstractVector{Tjd}, vr_i::AbstractVector{Tv}, vv_i::AbstractVector{Tv}; kwargs...) where {T<:Number, Tepoch<:Number, Tjd<:Number, Tv<:AbstractVector} -> KeplerianElements{Tepoch, T}, SMatrix{6, 6, T}

Fit a set of mean Keplerian elements for the J4 orbit propagator `j4d` using the osculating
elements represented by a set of position vectors `vr_i` [m] and a set of velocity vectors
`vv_i` [m / s] represented in an inertial reference frame at instants in the array `vjd`
[Julian Day].

!!! note
    The J4 orbit propagator `j4d` will be initialized with the Keplerian elements returned
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

- `KeplerianElements{Tepoch, T}`: Fitted Keplerian elements.
- `SMatrix{6, 6, T}`: Final covariance matrix of the least-square algorithm.

# Examples

```julia-repl
# Allocate a new J4 orbit propagator using a dummy Keplerian elements.
julia> j4d = j4_init(KeplerianElements{Float64, Float64}(0, 7000e3, 0, 0, 0, 0, 0));

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

julia> orb, P = fit_j4_mean_elements!(j4d, vjd, vr_i, vv_i)
ACTION:   Fitting the mean elements for the J4 propagator.
           Iteration        Position RMSE        Velocity RMSE           Total RMSE       RMSE Variation
                                     [km]             [km / s]                  [ ]
PROGRESS:         23              4.33863           0.00539962              4338.63         -2.90777e-06 %

(KeplerianElements{Float64, Float64}: Epoch = 2.46003e6 (2023-03-24T18:08:40.388), [0.16604866262321472 0.06643593039980161 … -3.8553206984036995e-5 0.00012403204426258087; 0.06643593040175688 0.26633435614867085 … -1.7942563356579497e-5 -1.95681107792859e-5; … ; -3.855320698470177e-5 -1.794256335557519e-5 … 4.3972013198895673e-7 -8.092682623169118e-8; 0.000124032044261828 -1.9568110780915995e-5 … -8.092682623078798e-8 1.2451901466558624e-7])

julia> orb
KeplerianElements{Float64, Float64}:
           Epoch :    2.46003e6 (2023-03-24T18:08:40.388)
 Semi-major axis : 7131.64       km
    Eccentricity :    0.00114298
     Inclination :   98.4366     °
            RAAN :  162.177      °
 Arg. of Perigee :  101.282      °
    True Anomaly :  258.693      °
```
"""
function fit_j4_mean_elements!(
    j4d::J4Propagator{Tepoch, T},
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
    cd = has_color ? string(_D) : ""
    cb = has_color ? string(_B) : ""
    cy = has_color ? string(_Y) : ""

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
        orb = update_j4_mean_elements_epoch!(j4d, initial_guess, epoch)

        r_i, v_i = kepler_to_rv(orb)
        x₁ = SVector{6, T}(r_i[1], r_i[2], r_i[3], v_i[1], v_i[2], v_i[3])
    else
        # In this case, we must find the closest osculating vector to the desired epoch.
        ~, id = findmin(abs.(vjd .- mean_elements_epoch))
        epoch = T(vjd[id])
        r_i   = vr_i[id]
        v_i   = vv_i[id]
        x₁    = SVector{6, T}(r_i[1], r_i[2], r_i[3], v_i[1], v_i[2], v_i[3],)
    end

    x₂ = x₁

    # Number of states in the input vector.
    num_states = 6

    # Number of observations in each instant.
    num_observations = 6

    # Covariance matrix.
    P = SMatrix{num_states, num_states, T}(I)

    # Variable to store the last residue.
    σ_i_₁ = nothing

    # Variable to store how many iterations the residue increased. This is used to account
    # for divergence.
    Δd = 0

    # Allocate the Jacobian matrix.
    J = zeros(T, num_observations, num_states)

    # Header.
    if verbose
        println("$(cy)ACTION:$(cd)   Fitting the mean elements for the J4 propagator.")
        @printf("          %s%10s %20s %20s %20s %20s%s\n", cy, "Iteration", "Position RMSE", "Velocity RMSE", "Total RMSE", "RMSE Variation", cd)
        @printf("          %s%10s %20s %20s %20s %20s%s\n", cb, "", "[km]", "[km / s]", "[ ]", "", cd)
        println()
    end

    # We need a reference to the covariance inverse because we will invert it and return
    # after the iterations.
    ΣJ′WJ = nothing

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
            y = vcat(vr_i[k], vv_i[k])

            # Initialize the SGP4 with the current estimated mean elements.
            orb = rv_to_kepler(x₁[1:3], x₁[4:6], epoch)
            j4_init!(j4d, orb)

            # Obtain the propagation time for this measurement.
            Δt = (vjd[k] - epoch) * 86400

            # Propagate the orbit.
            r̂_i, v̂_i = j4!(j4d, Δt)
            ŷ = vcat(r̂_i, v̂_i)

            # Compute the residue.
            b = y - ŷ

            # Compute the Jacobian in-place.
            _j4_jacobian!(
                J,
                j4d,
                Δt,
                x₁,
                ŷ;
                perturbation     = jacobian_perturbation,
                perturbation_tol = jacobian_perturbation_tol
            )

            # Convert the Jacobian matrix to a static matrix, leading to substantial
            # performance gains in the following computation.
            Js = SMatrix{num_observations, num_states, T}(J)

            # Accumulation.
            ΣJ′WJ += Js' * W * Js
            ΣJ′Wb += Js' * W * b
            σ_i   += b'  * W * b
            σp_i  += @views b[1:3]' * b[1:3]
            σv_i  += @views b[4:6]' * b[4:6]
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
        orb = update_j4_mean_elements_epoch!(j4d, orb, mean_elements_epoch)
    end

    # Initialize the propagator with the mean elements.
    j4_init!(j4d, orb)

    # Compute the final covariance.
    P = pinv(ΣJ′WJ)

    # Return the mean elements and the covariance.
    return orb, P
end

"""
    update_j4_mean_elements_epoch(orb::KeplerianElements, new_epoch::Union{Number, DateTime}) -> KepleriranElements

Update the epoch of the mean elements `orb` using a J4 orbit propagator to `new_epoch`,
which can be represented by a Julian Day or a `DateTime`.

!!! note
    This algorithm version will allocate a new J4 propagator with the default constants
    `j4c_egm2008`. If another set of constants are required, use the function
    [`update_j4osc_mean_elements_epoch!`](@ref) instead.

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

julia> update_j4_mean_elements_epoch(orb, DateTime("2023-01-02"))
KeplerianElements{Float64, Float64}:
           Epoch :    2.45995e6 (2023-01-02T00:00:00)
 Semi-major axis : 7190.98     km
    Eccentricity :    0.001111
     Inclination :   98.405    °
            RAAN :   90.9555   °
 Arg. of Perigee :  197.079    °
    True Anomaly :  127.293    °
```
"""
function update_j4_mean_elements_epoch(
    orb::KeplerianElements{Tepoch, T},
    new_epoch::Union{Number, DateTime}
) where {T<:Number, Tepoch<:Number}
    # Allocate the J4 propagator structure that will propagate the mean elements.
    j4d = J4Propagator{Tepoch, T}()

    # Assign the constants, which are used in the initialization.
    j4d.j4c = j4c_egm2008

    return update_j4_mean_elements_epoch!(j4d, orb, new_epoch)
end

"""
    update_j4_mean_elements_epoch!(j4d::J4Propagator, orb::KeplerianElements, new_epoch::Union{Number, DateTime}) -> KepleriranElements

Update the epoch of the mean elements `orb` using the propagator `j4d` to `new_epoch`, which
can be represented by a Julian Day or a `DateTime`.

!!! note
    The J4 orbit propagator `j4d` will be initialized with the Keplerian elements returned
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

# Allocate a new J4 orbit propagator using the created Keplerian elements. Notice that any
# set of Keplerian elements can be used here.
julia> j4d = j4_init(orb);

julia> update_j4_mean_elements_epoch!(j4d, orb, DateTime("2023-01-02"))
KeplerianElements{Float64, Float64}:
           Epoch :    2.45995e6 (2023-01-02T00:00:00)
 Semi-major axis : 7190.98     km
    Eccentricity :    0.001111
     Inclination :   98.405    °
            RAAN :   90.9555   °
 Arg. of Perigee :  197.079    °
    True Anomaly :  127.293    °
```
"""
function update_j4_mean_elements_epoch!(
    j4d::J4Propagator,
    orb::KeplerianElements,
    new_epoch::DateTime
)
    dt = datetime2julian(new_epoch)
    return update_j4_mean_elements_epoch!(j4d, orb, dt)
end

function update_j4_mean_elements_epoch!(
    j4d::J4Propagator,
    orb::KeplerianElements,
    new_epoch::Number
)
    # First, we need to initialize the J4 propagator with the mean elements.
    j4_init!(j4d, orb)

    # Now, we just need to propagate the orbit to the desired instant and obtain the mean
    # elements from the J4 propagator structure inside.
    Δt = (new_epoch - j4d.orb₀.t) * 86400
    j4!(j4d, Δt)
    orb = j4d.orbk

    # Finally, we initialize the propagator with the new set of mean elements.
    j4_init!(j4d, orb)

    return orb
end

############################################################################################
#                                    Private Functions
############################################################################################

"""
    _j4_jacobian!(J::AbstractMatrix{T}, j4d::J4OsculatingPropagator{Tepoch, T}, Δt::Number, x₁::SVector{6, T}, y₁::SVector{6, T}; kwargs...)) where {T<:Number, Tepoch<:Number}

Compute the J4 orbit propagator Jacobian by finite-differences using the propagator `j4d`
at instant `Δt` considering the input mean elements `x₁` that must provide the output vector
`y₁`. The result is written to the matrix `J`. Hence:

        ∂j4(x, Δt) │
    J = ────────── │
            ∂x     │ x = x₁

# Keywords

- `perturbation::T`: Initial state perturbation to compute the finite-difference:
    `Δx = x * perturbation`. (**Default** = 1e-3)
- `perturbation_tol::T`: Tolerance to accept the perturbation. If the computed perturbation
    is lower than `perturbation_tol`, we increase it until it absolute value is higher than
    `perturbation_tol`. (**Default** = 1e-7)
"""
function _j4_jacobian!(
    J::AbstractMatrix{T},
    j4d::J4Propagator{Tepoch, T},
    Δt::Number,
    x₁::SVector{6, T},
    y₁::SVector{6, T};
    perturbation::Number = T(1e-3),
    perturbation_tol::Number = T(1e-7)
) where {T<:Number, Tepoch<:Number}

    num_states = 6
    dim_obs    = 6

    # Auxiliary variables.
    x₂ = x₁

    @inbounds @views for j in 1:num_states
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
        orb = rv_to_kepler(x₂[1:3], x₂[4:6], j4d.orb₀.t)
        j4_init!(j4d, orb)
        r_i, v_i = j4!(j4d, Δt)
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

    return nothing
end
