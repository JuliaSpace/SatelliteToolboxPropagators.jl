## Description #############################################################################
#
# J4 orbit propagator algorithm.
#
# This algorithm propagates the orbit considering the secular perturbations of central
# body zonal harmonics as presented in [1, p. 647-654, 692-693] and [2], which is Kozai's
# method but neglecting long-periodic and short-periodic perturbations.
#
# The terms J2, J2², and J4 are considered, i.e. J6 is assumed to be 0. This can be used
# as a propagator of mean elements for mission analysis in which the satellite orbit is
# maintained.
#
## References ##############################################################################
#
# [1] Vallado, D. A (2013). Fundamentals of Astrodynamics and Applications. Microcosm
#     Press, Hawthorn, CA, USA.
#
# [2] Kozai, Y (1959). The Motion of a Close Earth Satellite. The Astronomical Journal,
#     v. 64, no. 1274, pp. 367 -- 377.
#
# [3] Hoots, F. R., Roehrich, R. L (1980). Models for Propagation of NORAD Elements Set.
#     Spacetrack Report No. 3.
#
# [4] Blitzer, L. Handbook of Orbital Perturbations. Astronautics 453. University of
#     Arizona.
#
# [5] https://www.mathworks.com/matlabcentral/fileexchange/43333-sun-synchronous-orbit-design
#
############################################################################################

export j4c_egm2008, j4c_egm1996, j4c_jgm02, j4c_jgm03
export j4c_egm2008_f32, j4c_egm1996_f32, j4c_jgm02_f32, j4c_jgm03_f32
export j4_init, j4_init!, j4, j4!
export fit_j4_mean_elements, fit_j4_mean_elements!
export update_j4_mean_elements_epoch, update_j4_mean_elements_epoch!

############################################################################################
#                                        Constants                                         #
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
    -1.6198975999169731e-6,
)

const j4c_egm2008_f32 = J4PropagatorConstants{Float32}(
    6378137.0,
    sqrt(3.986004415e14 / 6378137.0^3),
    0.0010826261738522227,
    -1.6198975999169731e-6,
)

# EGM-96 gravitational constants.
const j4c_egm1996 = J4PropagatorConstants(
    6378136.3, sqrt(3.986004415e14 / 6378136.3^3), 0.0010826266835531513, -1.619621591367e-6
)

const j4c_egm1996_f32 = J4PropagatorConstants{Float32}(
    6378136.3, sqrt(3.986004415e14 / 6378136.3^3), 0.0010826266835531513, -1.619621591367e-6
)

# JGM-02 gravitational constants.
const j4c_jgm02 = J4PropagatorConstants(
    6378136.3, sqrt(3.986004415e14 / 6378136.3^3), 0.0010826269256388149, -1.62042999e-6
)

const j4c_jgm02_f32 = J4PropagatorConstants{Float32}(
    6378136.3, sqrt(3.986004415e14 / 6378136.3^3), 0.0010826269256388149, -1.62042999e-6
)

# JGM-03 gravitational constants.
const j4c_jgm03 = J4PropagatorConstants(
    6378136.3, sqrt(3.986004415e14 / 6378136.3^3), 0.0010826360229829945, -1.619331205071e-6
)

const j4c_jgm03_f32 = J4PropagatorConstants{Float32}(
    6378136.3, sqrt(3.986004415e14 / 6378136.3^3), 0.0010826360229829945, -1.619331205071e-6
)

############################################################################################
#                                        Julia API                                         #
############################################################################################

# Define `copy` for the propagator structure.
let
    fields = fieldnames(J4Propagator)
    expressions = [:(new_j4d.$f = j4d.$f) for f in fields]

    @eval begin
        function Base.copy(
            j4d::J4Propagator{Tepoch, T}
        ) where {Tepoch <: Number, T <: Number}
            new_j4d = J4Propagator{Tepoch, T}()
            $(expressions...)
            return new_j4d
        end
    end
end

############################################################################################
#                                        Functions                                         #
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
    [`J4PropagatorConstants`](@ref)).
    (**Default** = `j4c_egm2008`)
"""
function j4_init(
    orb₀::KeplerianElements{Tanomaly, Tepoch, Tkepler};
    j4c::J4PropagatorConstants{T} = j4c_egm2008
) where {
    Tanomaly <: AbstractAnomaly,
    Tepoch <: Number,
    Tkepler <: AbstractFloat,
    T <: Number
}
    # Allocate the propagator structure.
    j4d = J4Propagator{Tepoch, T}()

    # Assign the constants, which are used in the initialization.
    j4d.j4c = j4c

    # Initialize the propagator and return.
    j4_init!(j4d, orb₀)

    return j4d
end

function j4_init(
    orb₀::KeplerianElements{Tanomaly, Tepoch, Tkepler};
    j4c::J4PropagatorConstants{Tj4c} = j4c_egm2008
) where {
    Tanomaly <: AbstractAnomaly,
    Tepoch <: Number,
    Tkepler <: Number,
    Tj4c <: Number
}
    T = promote_type(Tj4c, Tkepler)

    # Allocate the propagator structure.
    j4d = J4Propagator{Tepoch, T}()

    # Assign the constants, which are used in the initialization.
    j4d.j4c = J4PropagatorConstants{T}(j4c.R0, j4c.μm, j4c.J2, j4c.J4)

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
    j4d::J4Propagator{Tepoch, T}, orb₀::KeplerianElements
) where {Tepoch <: Number, T <: Number}
    # Unpack the gravitational constants to improve code readability.
    j4c = j4d.j4c
    R₀  = j4c.R0
    μm  = j4c.μm
    J₂  = j4c.J2
    J₄  = j4c.J4

    # Make sure the Keplerian elements use the mean anomaly.
    ke₀ = convert(KeplerianElements{MeanAnomaly, Tepoch, T}, orb₀)

    # Unpack orbit elements.
    a₀ = ke₀.semi_major_axis
    e₀ = ke₀.eccentricity
    i₀ = ke₀.inclination

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

    # Initial values and auxiliary variables.
    al₀ = a₀ / R₀           # .............................. Normalized semi-major axis [er]
    e₀² = e₀^2              # ..................................... Eccentricity squared [ ]
    p₀  = al₀ * (1 - e₀²)   # ....................................... Semi-latus rectum [er]
    p₀² = p₀^2              # .............................. Semi-latus rectum squared [er²]
    p₀⁴ = p₀^4              # ..................... Semi-latus rectum to the 4th power [er⁴]
    n₀  = μm / √(al₀^3)     # ............................ Unperturbed mean motion [rad / s]
    M₀  = mean_anomaly(ke₀) # ................................... Initial mean anomaly [rad]
    J₂² = J₂^2              # .......................................... J2 constant squared

    sin_i₀, cos_i₀ = sincos(T(i₀))

    sin_i₀² = sin_i₀^2
    sin_i₀⁴ = sin_i₀^4
    cos_i₀⁴ = cos_i₀^4
    β²      = (1 - e₀²)
    β       = √β²

    # We need to compute the perturbed mean motion that is used to calculate the
    # time-derivative of the orbital elements. This expression was obtained from [5] because
    # [2] does not show it completely.
    kn₂  = J₂ / p₀² * β
    kn₂₂ = J₂² / p₀⁴ * β
    kn₄  = J₄ / p₀⁴ * β

    n̄ =
        n₀ * (
            1 +
            (3//4) * kn₂ * (2 - 3sin_i₀²) +
            (3//128) *
            kn₂₂ *
            (
                120 + 64β - 40β² +
                (-240 - 192β + 40β²) * sin_i₀² +
                (105 + 144β + 25β²) * sin_i₀⁴
            ) - (45//128) * kn₄ * e₀² * (-8 + 40sin_i₀² - 35sin_i₀⁴)
        )

    # Some auxiliary variables to compute the perturbations.
    k̄₂  = n̄ * J₂ / p₀²
    k̄₂₂ = n̄ * J₂² / p₀⁴
    k₂₂ = n₀ * J₂² / p₀⁴
    k₄  = n₀ * J₄ / p₀⁴

    # TODO: Check J₄ perturbation term sign in RAAN time-derivative.
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

    ∂Ω =
        -(3//2) * k̄₂ * cos_i₀ +
        (3//32) * k̄₂₂ * cos_i₀ * (-36 - 4e₀² + 48β + (40 - 5e₀² - 72β) * sin_i₀²) +
        (15//32) * k₄ * cos_i₀ * (8 + 12e₀² - (14 + 21e₀²) * sin_i₀²)

    ∂ω =
        (3//4) * k̄₂ * (4 - 5sin_i₀²) +
        (3//128) *
        k̄₂₂ *
        (
            384 + 96e₀² - 384β +
            (-824 - 116e₀² + 1056β) * sin_i₀² +
            (430 - 5e₀² - 720β) * sin_i₀⁴
        ) - (15//16) * k₂₂ * e₀² * cos_i₀⁴ -
        (15//128) * k₄ * (64 + 72e₀² - (248 + 252e₀²) * sin_i₀² + (196 + 189e₀²) * sin_i₀⁴)

    # Initialize the propagator structure with the data.
    j4d.orb₀ = j4d.orbk = ke₀
    j4d.Δt   = 0
    j4d.∂Ω   = ∂Ω
    j4d.∂ω   = ∂ω
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
    [`J4PropagatorConstants`](@ref)).
    (**Default** = `j4c_egm2008`)

# Returns

- `SVector{3, T}`: Position vector [m] represented in the inertial frame at propagation
    instant.
- `SVector{3, T}`: Velocity vector [m / s] represented in the inertial frame at propagation
    instant.
- [`J4Propagator`](@ref): Structure with the initialized propagator.

# Remarks

The inertial frame in which the output is represented depends on which frame was used to
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

The inertial frame in which the output is represented depends on which frame was used to
generate the orbit parameters. Notice that the perturbation theory requires an inertial
frame with true equator.
"""
function j4!(j4d::J4Propagator{Tepoch, T}, t::Number) where {Tepoch <: Number, T <: Number}
    orbk = _j4_mean_elements!(j4d, t)

    # Compute the position and velocity vectors given the orbital elements.
    r_i_k, v_i_k = kepler_to_rv(orbk)

    # Return the position and velocity vector represented in the inertial reference frame.
    return r_i_k, v_i_k
end

# Propagate the mean elements of `j4d` to the instant `t` [s] measured from the epoch of the
# initial elements, update the propagator structure, and return the mean elements. The
# osculating propagator uses this function directly to avoid computing a state vector from
# the mean elements that it would discard.
function _j4_mean_elements!(
    j4d::J4Propagator{Tepoch, T}, t::Number
) where {Tepoch <: Number, T <: Number}
    # Unpack the variables.
    orb₀  = j4d.orb₀
    ∂Ω    = j4d.∂Ω
    ∂ω    = j4d.∂ω
    n̄     = j4d.n̄
    epoch = orb₀.epoch
    a₀    = orb₀.semi_major_axis
    e₀    = orb₀.eccentricity
    i₀    = orb₀.inclination
    Ω₀    = orb₀.raan
    ω₀    = orb₀.argument_of_periapsis
    M₀    = mean_anomaly(orb₀)

    # Time elapsed since epoch.
    Δt = T(t)

    # Propagate the orbital elements.
    Ω_k = mod(Ω₀ + ∂Ω * Δt, T(2π))
    ω_k = mod(ω₀ + ∂ω * Δt, T(2π))
    M_k = mod(M₀ + n̄ * Δt, T(2π))

    # Assemble the current mean elements.
    orbk = KeplerianElements{MeanAnomaly}(
        epoch + Tepoch(t) / 86400,
        a₀,
        e₀,
        i₀,
        Ω_k,
        ω_k,
        M_k
    )

    # Update the J4 orbit propagator structure.
    j4d.Δt   = Δt
    j4d.orbk = orbk

    return orbk
end

"""
    fit_j4_mean_elements(vjd::AbstractVector{Tjd}, vr_i::AbstractVector{Tv}, vv_i::AbstractVector{Tv}; kwargs...) where {Tjd<:Number, Tv<:AbstractVector} -> KeplerianElements{TrueAnomaly, Float64, Float64}, SMatrix{6, 6, Float64}

Fit a set of mean Keplerian elements for the J4 orbit propagator using the osculating
elements represented by a set of position vectors `vr_i` [m] and a set of velocity vectors
`vv_i` [m / s] represented in an inertial reference frame at instants in the array `vjd`
[Julian Day].

!!! note

    This algorithm version will allocate a new J4 propagator with the default constants
    `j4c_egm2008`. If another set of constants are required, use the function
    [`fit_j4_mean_elements!`](@ref) instead.

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
    vjd::AbstractVector{Tjd}, vr_i::AbstractVector{Tv}, vv_i::AbstractVector{Tv}; kwargs...
) where {Tjd <: Number, Tv <: AbstractVector}
    # Allocate the J4 propagator structure that will propagate the mean elements.
    j4d = J4Propagator{Float64, Float64}()

    # Assign the constants, which are used in the initialization.
    j4d.j4c = j4c_egm2008

    return fit_j4_mean_elements!(j4d, vjd, vr_i, vv_i; kwargs...)
end

"""
    fit_j4_mean_elements!(j4d::J4Propagator{Tepoch, T}, vjd::AbstractVector{Tjd}, vr_i::AbstractVector{Tv}, vv_i::AbstractVector{Tv}; kwargs...) where {T<:Number, Tepoch<:Number, Tjd<:Number, Tv<:AbstractVector} -> KeplerianElements{TrueAnomaly, Tepoch, T}, SMatrix{6, 6, T}

Fit a set of mean Keplerian elements for the J4 orbit propagator `j4d` using the osculating
elements represented by a set of position vectors `vr_i` [m] and a set of velocity vectors
`vv_i` [m / s] represented in an inertial reference frame at instants in the array `vjd`
[Julian Day].

!!! note

    The J4 orbit propagator `j4d` will be initialized with the Keplerian elements returned
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
# Allocate a new J4 orbit propagator using a dummy set of Keplerian elements.
julia> j4d = j4_init(KeplerianElements(0, 7000e3, 0, 0, 0, 0, 0));

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
    kwargs...,
) where {Tepoch <: Number, T <: Number, Tjd <: Number, Tv <: AbstractVector}
    return _fit_mean_elements!(j4d, vjd, vr_i, vv_i; kwargs...)
end

"""
    update_j4_mean_elements_epoch(orb::KeplerianElements, new_epoch::Union{Number, DateTime}) -> KeplerianElements

Update the epoch of the mean elements `orb` using a J4 orbit propagator to `new_epoch`,
which can be represented by a Julian Day or a `DateTime`.

!!! note

    This algorithm version will allocate a new J4 propagator with the default constants
    `j4c_egm2008`. If another set of constants are required, use the function
    [`update_j4_mean_elements_epoch!`](@ref) instead.

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
    orb::KeplerianElements{Tanomaly, Tepoch, T}, new_epoch::Union{Number, DateTime}
) where {Tanomaly <: AbstractAnomaly, Tepoch <: Number, T <: Number}
    # Allocate the J4 propagator structure that will propagate the mean elements.
    j4d = J4Propagator{Tepoch, T}()

    # Assign the constants, which are used in the initialization.
    j4d.j4c = j4c_egm2008

    return update_j4_mean_elements_epoch!(j4d, orb, new_epoch)
end

"""
    update_j4_mean_elements_epoch!(j4d::J4Propagator, orb::KeplerianElements, new_epoch::Union{Number, DateTime}) -> KeplerianElements

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
    j4d::J4Propagator, orb::KeplerianElements, new_epoch::DateTime
)
    dt = datetime2julian(new_epoch)
    return update_j4_mean_elements_epoch!(j4d, orb, dt)
end

function update_j4_mean_elements_epoch!(
    j4d::J4Propagator, orb::KeplerianElements, new_epoch::Number
)
    # First, we need to initialize the J4 propagator with the mean elements.
    j4_init!(j4d, orb)

    # Now, we just need to propagate the orbit to the desired instant and obtain the mean
    # elements from the J4 propagator structure inside.
    Δt = (new_epoch - j4d.orb₀.epoch) * 86400
    j4!(j4d, Δt)
    orb = j4d.orbk

    # Finally, we initialize the propagator with the new set of mean elements.
    j4_init!(j4d, orb)

    return orb
end

############################################################################################
#                                    Private Functions                                     #
############################################################################################

# Create a propagator that the finite-difference Jacobian can use as scratch space, so it
# does not clobber the propagator kept initialized by the fitting loop.
function _create_j4_fd_propagator(j4d::J4Propagator{Tepoch, T}) where {Tepoch, T}
    fd = J4Propagator{Tepoch, T}()
    fd.j4c = j4d.j4c
    return fd
end

function _create_j4_ad_propagator(j4d::J4Propagator{Tepoch, T}) where {Tepoch, T}
    tag = ForwardDiff.Tag{Nothing, T}
    D   = ForwardDiff.Dual{tag, T, 6}
    j4c = j4d.j4c

    ad = J4Propagator{Tepoch, D}()
    ad.j4c = J4PropagatorConstants{D}(D(j4c.R0), D(j4c.μm), D(j4c.J2), D(j4c.J4))
    return ad
end

function _j4_jacobian(
    ::FiniteDiffJacobian,
    j4d::J4Propagator{Tepoch, T},
    Δt::Number,
    x₁::SVector{6, T},
    y₁::SVector{6, T};
    perturbation::Number = T(1e-3),
    perturbation_tol::Number = T(1e-7),
    pd_ad::Union{Nothing, J4Propagator} = nothing,
    pd_fd::Union{Nothing, J4Propagator} = nothing,
) where {T <: Number, Tepoch <: Number}
    # The perturbed propagations below overwrite the propagator. Use the scratch propagator
    # when the caller provides one, so it can keep `j4d` initialized across the
    # measurements instead of reinitializing it for each one.
    epoch = j4d.orb₀.epoch
    fd::J4Propagator{Tepoch, T} = isnothing(pd_fd) ? j4d : pd_fd

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
        j4_init!(fd, orb)
        r_i, v_i = j4!(fd, Δt)
        y₂ = @SVector [r_i[1], r_i[2], r_i[3], v_i[1], v_i[2], v_i[3]]

        J[:, j] .= (y₂ .- y₁) ./ ϵ
        x₂ = setindex(x₂, x₁[j], j)
    end

    return SMatrix{6, 6, T}(J)
end

function _j4_jacobian(
    ::ForwardDiffJacobian,
    j4d::J4Propagator{Tepoch, T},
    Δt::Number,
    x₁::SVector{6, T},
    y₁::SVector{6, T};
    perturbation::Number = T(1e-3),
    perturbation_tol::Number = T(1e-7),
    pd_ad::Union{Nothing, J4Propagator} = nothing,
    pd_fd::Union{Nothing, J4Propagator} = nothing,
) where {T <: Number, Tepoch <: Number}
    epoch = j4d.orb₀.epoch
    N     = 6
    tag   = ForwardDiff.Tag{Nothing, T}
    D     = ForwardDiff.Dual{tag, T, N}

    # Declaring the type of the local makes the propagator concrete regardless of what
    # the caller provides. Otherwise, the unparameterised type in the keyword leaves
    # this call and the ones below dynamically dispatched.
    ad::J4Propagator{Tepoch, D} = isnothing(pd_ad) ? _create_j4_ad_propagator(j4d) : pd_ad

    seeds  = ntuple(i -> ForwardDiff.Partials(ntuple(j -> T(i == j), Val(N))), Val(N))
    x_dual = SVector{N, D}(ntuple(i -> D(x₁[i], seeds[i]), Val(N)))

    orb = rv_to_kepler(x_dual[SOneTo(3)], x_dual[StaticArrays.SUnitRange(4, 6)], epoch)
    j4_init!(ad, orb)
    r, v   = j4!(ad, Δt)
    y_dual = vcat(r, v)

    return SMatrix{6, N, T}(
        ntuple(k -> ForwardDiff.partials(y_dual[mod1(k, 6)], cld(k, 6)), Val(6 * N))
    )
end
