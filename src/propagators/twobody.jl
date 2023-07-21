# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#
# Description
# ==========================================================================================
#
#   Two-Body orbit propagator.
#
#   This algorithm considers a perfect Keplerian orbit. In other words, no perturbation is
#   considered during the propagation and the Earth is modeled as a perfect sphere.
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

export tbc_m0, tbc_m0_f32
export twobody_init, twobody_init!, twobody, twobody!

############################################################################################
#                                        Constants
############################################################################################

# Earth's standard gravitational parameter [m³/s²]
const tbc_m0     = 3.986004415e14
const tbc_m0_f32 = 3.986004415f14

############################################################################################
#                                        Functions
############################################################################################

"""
    twobody_init(orb₀::KeplerianElements; kwargs...) -> TwoBodyPropagator

Create and initialize the two-body propagator structure using the mean Keplerian elements
`orb₀`.

!!! note
    The type used in the propagation will be the same as used to define the gravitational
    constant `μ`.

# Arguments

- `orb₀::KeplerianElements`: Initial mean Keplerian elements [SI units].

# Keywords

- `μ::T`: Standard gravitational parameter of the central body [m³/s²].
    (**Default** = `tbc_m0`)
"""
function twobody_init(
    orb₀::KeplerianElements{Tepoch, Tkepler};
    μ::T = tbc_m0
) where {Tepoch<:Number, Tkepler<:Number, T<:Number}
    # Allocate the propagator structure.
    tbd = TwoBodyPropagator{Tepoch, T}()

    # Assign the constant, which are used in initialization.
    tbd.μ = μ

    # Initialize the propagator and return.
    twobody_init!(tbd, orb₀)

    return tbd
end

"""
    twobody_init!(tbd::TwoBodyPropagator, orb₀::KeplerianElements; kwargs...) -> Nothing

Initialize the two-body propagator structure `tbd` using the mean Keplerian elements `orb₀`.

!!! warning
    The propagation constant `μ::Number` in `tbd` will not be changed. Hence, it must be
    initialized.

# Arguments

- `orb₀::KeplerianElements`: Initial mean Keplerian elements [SI units].
"""
function twobody_init!(
    tbd::TwoBodyPropagator{Tepoch, T},
    orb₀::KeplerianElements
) where {Tepoch<:Number, T<:Number}
    # Compute the mean motion using the semi-major axis.
    n₀ = √(tbd.μ / T(orb₀.a)^3)

    # Compute the initial mean anomaly.
    M₀ = true_to_mean_anomaly(T(orb₀.e), T(orb₀.f))

    # Create and return the two-body orbit propagator structure.
    tbd.orb₀ = orb₀
    tbd.orbk = orb₀
    tbd.Δt   = 0
    tbd.M₀   = M₀
    tbd.n₀   = n₀

    return nothing
end

"""
    twobody(Δt::Number, orb₀::KeplerianElements; kwargs...)

Initialize the two-body propagator structure using the input elements `orb₀` and propagate
the orbit until the time Δt [s].

!!! note
    The type used in the propagation will be the same as used to define the gravitational
    constant `μ`.

# Arguments

- `orb₀::KeplerianElements`: Initial mean Keplerian elements [SI units].

# Keywords

- `μ::T`: Standard gravitational parameter of the central body [m³/s²].
    (**Default** = `tbc_m0`)

# Returns

- `SVector{3, T}`: Position vector [m] represented in the inertial frame at propagation
    instant.
- `SVector{3, T}`: Velocity vector [m / s] represented in the inertial frame at propagation
    instant.
- [`TwoBodyPropagator`](@ref): Structure with the initialized parameters.

# Remarks

The inertial frame in which the output is represented depends on which frame it was used to
generate the orbit parameters. Notice that the perturbation theory requires an inertial
frame with true equator.
"""
function twobody(Δt::Number, orb₀::KeplerianElements; μ::T = tbc_m0) where T<:Number
    tbd = twobody_init(orb₀; μ = μ)
    r_i, v_i = twobody!(tbd, Δt)
    return r_i, v_i, tbd
end

"""
    twobody!(tbd::TwoBodyPropagator{Tepoch, T}, t::Number) where {Tepoch, T}

Propagate the orbit defined in `tbd` (see [`TwoBodyPropagator`](@ref)) to `t` [s] after the
epoch of the input mean elements in `tbd`.

!!! note
    The internal values in `tbd` will be modified.

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
function twobody!(tbd::TwoBodyPropagator{Tepoch, T}, t::Number) where {Tepoch<:Number, T<:Number}
    # Unpack.
    orb₀ = tbd.orb₀

    # Time elapsed since epoch.
    epoch  = orb₀.t
    Δt     = T(t)

    # Propagate the orbital elements.
    a_k = orb₀.a
    e_k = orb₀.e
    i_k = orb₀.i
    Ω_k = orb₀.Ω
    ω_k = orb₀.ω
    M_k = tbd.M₀ + tbd.n₀ * Δt

    # Convert the mean anomaly to true anomaly.
    f_k = mean_to_true_anomaly(e_k, M_k)

    # Assemble the current mean elements.
    orbk = KeplerianElements(epoch + Tepoch(Δt) / 86400, a_k, e_k, i_k, Ω_k, ω_k, f_k)

    # Compute the position and velocity vectors given the orbital elements.
    r_i_k, v_i_k = kepler_to_rv(orbk)

    # Update the propagator structure.
    tbd.Δt   = Δt
    tbd.orbk = orbk

    # Return the position and velocity vector represented in the inertial
    # reference frame.
    return r_i_k, v_i_k
end
