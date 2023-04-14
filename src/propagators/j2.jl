# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#
# Description
# ==========================================================================================
#
#   J2 orbit propagator algorithm.
#
#   This algorithm propagates the orbit considering the perturbed two-body equations as
#   presented in [1, p. 690-692]. It uses the first-order approximation of Kepler's problem,
#   considering the effects of secular gravitational and drag perturbations.
#
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#
# References
# ==========================================================================================
#
#   [1] Vallado, D. A (2013). Fundamentals of Astrodynamics and Applications. Microcosm
#       Press, Hawthorn, CA, USA.
#
#   [2] Wertz, J. R (1978). Spacecraft attitude determination and control. Kluwer Academic
#       Publishers, Dordrecht, The Netherlands.
#
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

export j2c_egm08, j2c_egm96, j2c_jgm02, j2c_jgm03
export j2c_egm08_f32, j2c_egm96_f32, j2c_jgm02_f32, j2c_jgm03_f32
export j2_init, j2_init!, j2!

############################################################################################
#                                           TODO
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
#                                        Constants
############################################################################################

# These constants were obtained from the GFC files. Remember that:
#
#   J_n = -C_n,0 * √(2n + 1)
#

# EGM-08 gravitational constants.
const j2c_egm08 = J2PropagatorConstants(
    6378137.0,
    √(3.986004415e14 / 6378137.0^3),
    0.0010826261738522227
)

const j2c_egm08_f32 = J2PropagatorConstants{Float32}(
    6378137.0,
    √(3.986004415e14 / 6378137.0^3),
    0.0010826261738522227
)

# EGM-96 gravitational constants.
const j2c_egm96 = J2PropagatorConstants(
    6378136.3,
    √(3.986004415e14 / 6378136.3^3),
    0.0010826266835531513
)

const j2c_egm96_f32 = J2PropagatorConstants{Float32}(
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
#                                        Functions
############################################################################################

"""
    j2_init(orb₀::KeplerianElements, dn_o2::Number = 0, ddn_o6::Number = 0; kwargs...) where T<:Number -> J2Propagator

Create and initialize the J2 orbit propagator structure using the mean Keplerian elements
`orb₀`.

!!! note
    The type used in the propagation will be the same as used to define the constants in the
    structure `j2c`.

# Arguments

- `orb₀::KeplerianElements`: Initial mean Keplerian elements [SI units].
- `dn_o2::Number`: First time derivative of the mean motion divided by two [rad/s^2].
    (**Default** = 0)
- `ddn_o6::Number`: Second time derivative of the mean motion divided by six [rad/s^3].
    (**Default** = 0)

# Keywords

- `j2c::J2PropagatorConstants`: J2 orbit propagator constants (see
  [`J2PropagatorConstants`](@ref)). (**Default** = `j2c_egm08`)
"""
function j2_init(
    orb₀::KeplerianElements{Tepoch, Tkepler},
    dn_o2::Number = 0,
    ddn_o6::Number = 0;
    j2c::J2PropagatorConstants{T} = j2c_egm08
) where {Tepoch<:Number, Tkepler<:Number, T<:Number}
    # Allocate the propagator structure.
    j2d = J2Propagator{Tepoch, T}()

    # Assign the constants, which are used in the initialization.
    j2d.j2c = j2c

    # Initialize the propagator and return.
    j2_init!(j2d, orb₀, dn_o2, ddn_o6)

    return j2d
end

"""
    j2_init!(j2d::J2Propagator, orb₀::KeplerianElements, dn_o2::Number = 0, ddn_o6::Number = 0) -> Nothing

Initialize the J2 orbit propagator structure `j2d` using the mean Keplerian elements `orb₀`.

!!! warning
    The propagation constants `j2c::J2PropagatorConstants` in `j2d` will not be changed.
    Hence, they must be initialized.

# Arguments

- `orb₀::KeplerianElements`: Initial mean Keplerian elements [SI units].
- `dn_o2::Number`: First time derivative of the mean motion divided by two [rad/s^2].
    (**Default** = 0)
- `ddn_o6::Number`: Second time derivative of the mean motion divided by six [rad/s^3].
    (**Default** = 0)
"""
function j2_init!(
    j2d::J2Propagator{Tepoch, T},
    orb₀::KeplerianElements,
    dn_o2::Number = 0,
    ddn_o6::Number = 0
) where {Tepoch<:Number, T<:Number}
    # Unpack the gravitational constants to improve code readability.
    j2c = j2d.j2c
    R0  = j2c.R0
    μm  = j2c.μm
    J2  = j2c.J2

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
    n₀  = μm / al₀^(T(3 / 2))          # ................... Unperturbed mean motion [rad/s]
    p₀  = al₀ * (1 - e₀²)              # ............................ Semi-latus rectum [er]
    p₀² = p₀^2                         # ................... Semi-latus rectum squared [er²]
    M₀  = true_to_mean_anomaly(e₀, f₀) # ........................ Initial mean anomaly [rad]
    dn  = 2T(dn_o2)                    # ..... . Time-derivative of the mean motion [rad/s²]

    sin_i₀, cos_i₀ = sincos(T(i₀))
    sin_i₀² = sin_i₀^2

    # We need to compute the "mean" mean motion that is used to calculate the first-order
    # time derivative of the orbital elements.
    #
    # NOTE: Description of J2 propagator in [1, p. 691].
    #
    # Using the equations in [1, p. 691], we could not match the results from STK as
    # mentioned in the issue:
    #
    #   https://github.com/JuliaSpace/SatelliteToolbox.jl/issues/91
    #
    # After analyzing the perturbation equations, it turns out that the time-derivative
    # depends on the mean motion instead of the unperturbed mean motion as in the algorithm
    # 65. We can see this by looking at the algorithm in Kozai's method in [1, p. 693].
    n̄ = n₀ * (1 + T(3 / 4) * J2 / p₀² * √(1 - e₀²) * (2 - 3sin_i₀²))

    # First-order time-derivative of the orbital elements.
    δa = -T(2 / 3) * al₀ * dn / n₀
    δe = -T(2 / 3) * (1 - T(e₀)) * dn / n₀
    δΩ = -T(3 / 2) * n̄ * J2 / p₀² * cos_i₀
    δω = +T(3 / 4) * n̄ * J2 / p₀² * (4 - 5sin_i₀²)

    # Initialize the propagator structure with the data.
    j2d.orb₀   = j2d.orbk = orb₀
    j2d.dn_o2  = dn_o2
    j2d.ddn_o6 = ddn_o6
    j2d.Δt     = 0
    j2d.al₀    = al₀
    j2d.M₀     = M₀
    j2d.δa     = δa
    j2d.δe     = δe
    j2d.δΩ     = δΩ
    j2d.δω     = δω
    j2d.n̄      = n̄

    return nothing
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
    j2c    = j2d.j2c
    al₀    = j2d.al₀
    M₀     = j2d.M₀
    δa     = j2d.δa
    δe     = j2d.δe
    δΩ     = j2d.δΩ
    δω     = j2d.δω
    n̄      = j2d.n̄
    dn_o2  = j2d.dn_o2
    ddn_o6 = j2d.ddn_o6
    R0     = j2c.R0
    epoch  = orb₀.t
    e₀     = orb₀.e
    i₀     = orb₀.i
    Ω₀     = orb₀.Ω
    ω₀     = orb₀.ω

    # Time from epoch to propagate the orbit.
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
    f_k = mean_to_true_anomaly(e_k, M_k)

    # Make sure that eccentricity is not lower than 0.
    e_k = max(e_k, T(0))

    # Assemble the current mean elements.
    orbk = KeplerianElements(epoch + Tepoch(Δt) / 86400, al_k * R0, e_k, i_k, Ω_k, ω_k, f_k)

    # Compute the position and velocity vectors given the orbital elements.
    r_i_k, v_i_k = kepler_to_rv(orbk)

    # Update the J2 orbit propagator structure.
    j2d.Δt   = Δt
    j2d.orbk = orbk

    # Return the position and velocity vector represented in the inertial reference frame.
    return r_i_k, v_i_k
end
