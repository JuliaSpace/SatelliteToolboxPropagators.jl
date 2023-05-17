# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#
# Description
# ==========================================================================================
#
#   J2 osculating orbit propagator algorithm.
#
#   This algorithm propagates the orbit considering the secular and short-period
#   perturbations introduced by the J2 gravitational term. The algorithm is based on Kwok
#   version as indicated in [1, p. 708-710].
#
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#
# References
# ==========================================================================================
#
#   [1] Vallado, D. A (2013). Fundamentals of Astrodynamics and Applications. Microcosm
#       Press, Hawthorn, CA, USA.
#
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

export j2osc_init, j2osc_init!, j2osc, j2osc!

############################################################################################
#                                        Functions
############################################################################################

"""
    j2osc_init(orb₀::KeplerianElements, dn_o2::Number = 0, ddn_o6::Number = 0; kwargs...) where T<:Number -> J2OsculatingPropagator

Create and initialize the J2 osculating orbit propagator structure using the mean Keplerian
elements `orb₀`.

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
function j2osc_init(
    orb₀::KeplerianElements{Tepoch, Tkepler},
    dn_o2::Number = 0,
    ddn_o6::Number = 0;
    j2c::J2PropagatorConstants{T} = j2c_egm08
) where {Tepoch<:Number, Tkepler<:Number, T<:Number}
    # Allocate the J2 propagator structure that will propagate the mean elements.
    j2d = J2Propagator{Tepoch, T}()

    # Assign the constants, which are used in the initialization.
    j2d.j2c = j2c

    # Allocate the J2 osculating propagator structure.
    j2oscd = J2OsculatingPropagator{Tepoch, T}()
    j2oscd.j2d = j2d

    # Initialize the propagator and return.
    j2osc_init!(j2oscd, orb₀, dn_o2, ddn_o6)

    return j2oscd
end

"""
    j2osc_init!(j2oscd::J2OsculatingPropagator, orb₀::KeplerianElements, dn_o2::Number = 0, ddn_o6::Number = 0) -> Nothing

Initialize the J2 osculating orbit propagator structure `j2oscd` using the mean Keplerian
elements `orb₀`.

!!! warning
    The propagation constants `j2c::J2PropagatorConstants` in `j2oscd.j2d` will not be
    changed. Hence, they must be initialized.

# Arguments

- `orb₀::KeplerianElements`: Initial mean Keplerian elements [SI units].
- `dn_o2::Number`: First time derivative of the mean motion divided by two [rad/s^2].
    (**Default** = 0)
- `ddn_o6::Number`: Second time derivative of the mean motion divided by six [rad/s^3].
    (**Default** = 0)
"""
function j2osc_init!(
    j2oscd::J2OsculatingPropagator{Tepoch, T},
    orb₀::KeplerianElements,
    dn_o2::Number = 0,
    ddn_o6::Number = 0
) where {Tepoch<:Number, T<:Number}
    # Initialize the J2 propagator that will propagate the mean elements.
    j2_init!(j2oscd.j2d, orb₀, dn_o2, ddn_o6)

    # Call the propagation one time to update the osculating elements.
    j2osc!(j2oscd, 0)

    return nothing
end

"""
    j2osc(Δt::Number, orb₀::KeplerianElements, dn_o2::Number = 0, ddn_o6::Number = 0; kwargs...)

Initialize the J2 osculating propagator structure using the input elements `orb₀` and
propagate the orbit until the time Δt [s].

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

- `j2c::J2PropagatorConstants{T}`: J2 orbit propagator constants (see
  [`J2PropagatorConstants`](@ref)). (**Default** = `j2c_egm08`)

# Returns

- `SVector{3, T}`: Position vector [m] represented in the inertial frame at propagation
    instant.
- `SVector{3, T}`: Velocity vector [m / s] represented in the inertial frame at propagation
    instant.
- [`J2OsculatingPropagator`](@ref): Structure with the initialized parameters.

# Remarks

The inertial frame in which the output is represented depends on which frame it was used to
generate the orbit parameters. Notice that the perturbation theory requires an inertial
frame with true equator.
"""
function j2osc(
    Δt::Number,
    orb₀::KeplerianElements,
    dn_o2::Number = 0,
    ddn_o6::Number = 0;
    j2c::J2PropagatorConstants{T} = j2c_egm08
) where T<:Number
    j2oscd = j2osc_init(orb₀, dn_o2, ddn_o6; j2c = j2c)
    r_i, v_i = j2osc!(j2oscd, Δt)
    return r_i, v_i, j2oscd
end

"""
    j2osc!(j2oscd::J2OsculatingPropagator{Tepoch, T}, t::Number) where {Tepoch, T} -> SVector{3, T}, SVector{3, T}

Propagate the orbit defined in `j2oscd` (see [`J2OsculatingPropagator`](@ref)) to `t` [s]
after the epoch of the input mean elements in `j2d`.

!!! note
    The internal values in `j2oscd` will be modified.

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
function j2osc!(j2oscd::J2OsculatingPropagator{Tepoch, T}, t::Number) where {Tepoch<:Number, T<:Number}
    # First, we need to propagate the mean elements since they are necessary to compute the
    # short-periodic perturbations.
    j2d = j2oscd.j2d
    j2!(j2d, t)

    # Unpack the propagator constants.
    j2c = j2d.j2c
    R0  = j2c.R0
    μm  = j2c.μm
    J2  = j2c.J2

    # Time from epoch to propagate the orbit.
    Δt = T(t)

    # Obtain the mean elements at this time instant.
    mean_orbk = j2d.orbk

    a_k  = mean_orbk.a
    e_k  = mean_orbk.e
    e_k² = e_k * e_k
    i_k  = mean_orbk.i
    Ω_k  = mean_orbk.Ω
    ω_k  = mean_orbk.ω
    f_k  = mean_orbk.f
    M_k  = true_to_mean_anomaly(e_k, f_k)
    p_k  = a_k * (1 - e_k²)
    p_k² = p_k * p_k
    u_k  = ω_k + f_k

    # Auxiliary variables to reduce the computational burden.
    KJ2 = J2 * R0 * R0

    sin_i_k, cos_i_k         = sincos(i_k)
    sin_f_k, cos_f_k         = sincos(f_k)
    sin_2u_k, cos_2u_k       = sincos(2u_k)
    sin_2ω_f_k, cos_2ω_f_k   = sincos(2ω_k + f_k)
    sin_2ω_3f_k, cos_2ω_3f_k = sincos(2ω_k + 3f_k)

    sin_i_k²  = sin_i_k * sin_i_k
    cos_i_k²  = cos_i_k * cos_i_k
    e_cos_f_k = e_k * cos_f_k
    e_sin_f_k = e_k * sin_f_k

    aux1 = 3cos_2u_k + 3e_k * cos_2ω_f_k + e_k * cos_2ω_3f_k
    aux2 = √(1 - e_k²)
    aux3 = 3cos_i_k² - 1
    aux4 = -aux3 / (1 + aux2)

    # Compute the short-periodic perturbations considering only the J2 gravitational term.
    δisp_k = +KJ2 * sin_i_k * cos_i_k / (4p_k²) * aux1

    δpsp_k = +KJ2 * sin_i_k² / (2p_k) * aux1

    δΩsp_k = -KJ2 * cos_i_k / (4p_k²) * (
        6 * (f_k - M_k + e_sin_f_k) - 3sin_2u_k - 3e_k * sin_2ω_f_k - e_k * sin_2ω_3f_k
    )

    δrsp_k = -KJ2 / (4p_k) * (
        aux3 * (2aux2 / (1 + e_cos_f_k) + e_cos_f_k / (1 + aux2) + 1) - sin_i_k² * cos_2u_k
    )

    δṙsp_k = +KJ2 * √μm / (4 * √(p_k^5)) * (
        aux3 * e_sin_f_k * (aux2 + ((1 + e_cos_f_k)^2) / (1 + aux2)) -
        2sin_i_k² * (1 - e_cos_f_k)^2 * sin_2u_k
    )

    δusp_k = +KJ2 / (8p_k²) * (
        (6 - 30cos_i_k²) * (f_k - M_k) +
        4e_sin_f_k * (1 - 6cos_i_k² + aux4) +
        aux4 * e_k² * sin(2f_k) +
        (5cos_i_k² - 2) * (2e_k) * sin_2ω_f_k +
        (7cos_i_k² - 1) * sin_2u_k +
        2cos_i_k² * e_k * sin_2ω_3f_k
    )

    r_k = p_k / (1 + e_cos_f_k)
    ṙ_k = √(μm / p_k) * e_sin_f_k

    r_osc_k = r_k + δrsp_k
    ṙ_osc_k = ṙ_k + δṙsp_k
    p_osc_k = p_k + δpsp_k

    A_k = p_osc_k / r_osc_k - 1
    B_k = √(p_osc_k / μm) * ṙ_osc_k

    e_osc_k² = A_k^2 + B_k^2
    e_osc_k  = √e_osc_k²
    a_osc_k  = p_osc_k / (1 - e_osc_k²)
    i_osc_k  = i_k + δisp_k
    Ω_osc_k  = Ω_k + δΩsp_k
    u_osc_k  = u_k + δusp_k
    f_osc_k  = atan(B_k, A_k)
    ω_osc_k  = u_osc_k - f_osc_k

    # Assemble the current osculating elements.
    orbk = KeplerianElements(
        j2d.orb₀.t + Tepoch(Δt) / 86400,
        a_osc_k,
        e_osc_k,
        i_osc_k,
        Ω_osc_k,
        ω_osc_k,
        f_osc_k
    )

    # Compute the position and velocity considering the osculating elements.
    r_i_k, v_i_k = kepler_to_rv(orbk)

    # Update the J2 orbit propagator structure.
    j2oscd.Δt   = T(t)
    j2oscd.orbk = orbk

    return r_i_k, v_i_k
end
