# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#
# Description
# ==========================================================================================
#
#   J4 orbit propagator algorithm.
#
#   This algorithm propagates the orbit considering the secular perturbations of central
#   body zonal harmonics as presented in [1, p. 647-654, 692-692], which is Kozai's method
#   but neglecting long-periodic and short-periodic perturbations.
#
#   The terms J2, J2², and J4 are considered, i.e. the J6 is assumed to be 0.  The effect of
#   the drag is also taken into account. This can be used as a propagator of mean elements
#   for mission analysis in which the satellite orbit is maintained.
#
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#
# References
# ==========================================================================================
#
#   [1] Vallado, D. A (2013). Fundamentals of Astrodynamics and Applications. Microcosm
#       Press, Hawthorn, CA, USA.
#
#   [2] Hoots, F. R., Roehrich, R. L (1980). Models for Propagation of NORAD Elements Set.
#       Spacetrack Report No. 3.
#
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

export j4c_egm08, j4c_egm96, j4c_jgm02, j4c_jgm03
export j4c_egm08_f32, j4c_egm96_f32, j4c_jgm02_f32, j4c_jgm03_f32
export j4_init, j4!

############################################################################################
#                                        Constants
############################################################################################

# These constants were obtained from the GFC files. Remember that:
#
#   J_n = -C_n,0 * √(2n + 1)
#

# EGM-08 gravitational constants.
const j4c_egm08 = J4PropagatorConstants(
    6378137.0,
    sqrt(3.986004415e14 / 6378137.0^3),
    0.0010826261738522227,
    -1.6198975999169731e-6
)

const j4c_egm08_f32 = J4PropagatorConstants{Float32}(
    6378137.0,
    sqrt(3.986004415e14 / 6378137.0^3),
    0.0010826261738522227,
    -1.6198975999169731e-6
)

# EGM-96 gravitational constants.
const j4c_egm96 = J4PropagatorConstants(
    6378136.3,
    sqrt(3.986004415e14 / 6378136.3^3),
    0.0010826266835531513,
    -1.619621591367e-6
)

const j4c_egm96_f32 = J4PropagatorConstants{Float32}(
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
    j4_init(orb₀::KeplerianElements, dn_o2::Number = 0, ddn_o6::Number = 0; kwargs...) where T<:Number -> J4Propagator

Initialize the J4 orbit propagator algorithm using the mean Keplerian elements `orb₀`.

!!! note
    The type used in the propagation will be the same as used to define the constants in the
    structure `j4c`.

# Arguments

- `orb₀::KeplerianElements`: Initial mean Keplerian elements [SI units].
- `dn_o2::Number`: First time derivative of the mean motion divided by two [rad/s^2].
    (**Default** = 0)
- `ddn_o6::Number`: Second time derivative of the mean motion divided by six [rad/s^3].
    (**Default** = 0)

# Keywords

- `j4c::J4PropagatorConstants`: J4 orbit propagator constants (see
  [`J4PropagatorConstants`](@ref)). (**Default** = `j4c_egm08`)
"""
function j4_init(
    orb₀::KeplerianElements,
    dn_o2::Number = 0,
    ddn_o6::Number = 0;
    j4c::J4PropagatorConstants{T} = j4c_egm08
) where T<:Number
    # Unpack the gravitational constants to improve code readability.
    R0 = j4c.R0
    μm = j4c.μm
    J2 = j4c.J2
    J4 = j4c.J4

    # Unpack orbit elements.
    epoch = orb₀.t
    a₀    = T(orb₀.a)
    e₀    = T(orb₀.e)
    i₀    = T(orb₀.i)
    Ω₀    = T(orb₀.Ω)
    ω₀    = T(orb₀.ω)
    f₀    = T(orb₀.f)

    # Initial values and auxiliary variables.
    al₀ = a₀ / R0                      # ................... Normalized semi-major axis [er]
    e₀² = e₀^2                         # .......................... Eccentricity squared [ ]
    p₀  = al₀ * (1 - e₀^2)             # ............................ Semi-latus rectum [er]
    p₀² = p₀^2                         # ................... Semi-latus rectum squared [er²]
    p₀⁴ = p₀^4                         # .......... Semi-latus rectum to the 4th power [er⁴]
    n₀  = μm / al₀^(T(3 / 2))          # ................. Unperturbed mean motion [rad / s]
    M₀  = true_to_mean_anomaly(e₀, f₀) # ........................ Initial mean anomaly [rad]
    dn  = 2T(dn_o2)                    # ..... Time-derivative of the mean motion [rad / s²]
    J2² = J2^2                         # ............................... J2 constant squared

    sin_i₀, cos_i₀ = sincos(T(i₀))

    sin_i₀² = sin_i₀^2
    sin_i₀⁴ = sin_i₀^4
    aux     = (1 - e₀²)
    saux    = sqrt(aux)

    # We need to compute the "mean" mean motion that is used to calculate the first-order
    # time derivative of the orbital elements.
    #
    # NOTE: Description of J4 propagator in [1, p. 648-653].
    #
    # Using the equations in [1], we could not match the results from STK. After analyzing
    # the perturbation equations, it turns out that the time-derivative depends on the mean
    # motion instead of the unperturbed mean motion. We can see this by looking at the
    # algorithm in Kozai's method in [1, p. 693].
    #
    # Notice that using the full expression here, with the J2² and J4 terms, yields a
    # solution with much higher error compared with STK result.
    āl = al₀ * (1 - T(3 / 4) * J2 / p₀² * saux * (2 - 3sin_i₀²))
    p̄  = āl   * aux
    n̄  = n₀  * (1 + T(3 / 4) * J2 / p̄^2 * saux * (2 - 3sin_i₀²))

    # First-order time-derivative of the orbital elements.
    δa = -T(2 / 3) * al₀ * dn / n₀
    δe = -T(2 / 3) * (1 - e₀) * dn / n₀

    # TODO: Check J4 perturbation term sign.
    #
    # We needed to flip the J4 perturbation term sign to obtain values that match those of
    # STK. However, this modification does not seem right if we observe the RAAN secular
    # perturbation term in SGP4 orbit propagator [2, p. 16]. For more information, see:
    #
    #   https://github.com/JuliaSpace/SatelliteToolbox.jl/issues/91
    #
    δΩ = -T( 3 / 2  ) * n̄ * J2  / p₀² * cos_i₀ +
          T( 3 / 32 ) * n̄ * J2² / p₀⁴ * cos_i₀ * (-36 -  4e₀² + 48saux + (40 - 5e₀² - 72saux) * sin_i₀²) +
          T(15 / 32 ) * n̄ * J4  / p₀⁴ * cos_i₀ * (  8 + 12e₀² - (14 + 21e₀²) * sin_i₀²)

    δω = +T( 3 / 4  ) * n̄ * J2  / p₀² * (4 - 5sin_i₀²) +
          T( 9 / 384) * n̄ * J2² / p₀⁴ * (192 + 56e₀² - 192saux + (-172 + 288saux) * sin_i₀² + e₀² * sin_i₀⁴) -
          T(15 / 128) * n̄ * J4  / p₀⁴ * (64 + 72e₀² - (248 + 252e₀²) * sin_i₀² + (196 + 189e₀²) * sin_i₀⁴)

    # Create a new instant of `KeplerianElements` with the converted types.
    orb = KeplerianElements(epoch, a₀, e₀, i₀, Ω₀, ω₀, f₀)

    # Create the output structure with the data.
    return J4Propagator(
        orb,
        orb,
        T(dn_o2),
        T(ddn_o6),
        j4c,
        T(0),
        al₀,
        M₀,
        δa,
        δe,
        δΩ,
        δω,
        n̄
    )
end

"""
    j4!(j4d::J4Propagator{Tepoch, T}, t::Number) where {Tepoch, T} -> SVector{3, T}, SVector{3, T}

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
function j4!(j4d::J4Propagator{Tepoch, T}, t::Number) where {Tepoch, T}
    # Unpack the variables.
    orb₀   = j4d.orb₀
    j4c    = j4d.j4c
    al₀    = j4d.al₀
    M₀     = j4d.M₀
    δa     = j4d.δa
    δe     = j4d.δe
    δΩ     = j4d.δΩ
    δω     = j4d.δω
    n̄      = j4d.n̄
    dn_o2  = j4d.dn_o2
    ddn_o6 = j4d.ddn_o6
    R0     = j4c.R0
    epoch  = orb₀.t
    e₀     = orb₀.e
    i₀     = orb₀.i
    Ω₀     = orb₀.Ω
    ω₀     = orb₀.ω

    # Time elapsed since epoch.
    Δt = T(t)

    # Propagate the orbital elements.
    al_k = al₀ + δa * Δt
    e_k  = e₀  + δe * Δt
    i_k  = i₀
    Ω_k  = mod(Ω₀ + δΩ * Δt, T(2π))
    ω_k  = mod(ω₀ + δω * Δt, T(2π))

    # The mean anomaly update equation can be seen in [1, p. 693]. However, we add the terms
    # related with the time-derivative of the mean motion as in [1, p. 692].
    M_k = mod(@evalpoly(Δt, M₀, n̄, dn_o2, ddn_o6), T(2π))

    # Convert the mean anomaly to the true anomaly.
    f_k  = mean_to_true_anomaly(e_k, M_k)

    # Make sure that eccentricity is not lower than 0.
    e_k = max(e_k, T(0))

    # Assemble the current mean elements.
    orbk = KeplerianElements(epoch + Tepoch(Δt) / 86400, al_k * R0, e_k, i_k, Ω_k, ω_k, f_k)

    # Compute the position and velocity vectors given the orbital elements.
    r_i_k, v_i_k = kepler_to_rv(orbk)

    # Update the J2 orbit propagator structure.
    j4d.Δt   = Δt
    j4d.orbk = orbk

    # Return the position and velocity vector represented in the inertial reference frame.
    return r_i_k, v_i_k
end
