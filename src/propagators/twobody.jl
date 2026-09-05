## Description #############################################################################
#
#   Two-Body orbit propagator.
#
#   This algorithm considers a perfect Keplerian orbit. In other words, no perturbation is
#   considered during the propagation and the Earth is modeled as a perfect sphere.
#
## References ##############################################################################
#
# [1] Vallado, D. A (2013). Fundamentals of Astrodynamics and Applications. Microcosm
#     Press, Hawthorn, CA, USA.
#
############################################################################################

export tbc_m0, tbc_m0_f32
export twobody_init, twobody_init!, twobody, twobody!

############################################################################################
#                                        Constants                                         #
############################################################################################

# Earth's standard gravitational parameter [m³/s²]
const tbc_m0     = 3.986004415e14
const tbc_m0_f32 = 3.986004415f14

############################################################################################
#                                        Julia API                                         #
############################################################################################

# Define `copy` for the propagator structure.
let
    fields = fieldnames(TwoBodyPropagator)
    expressions = [:(new_tbd.$f = tbd.$f) for f in fields]

    @eval begin
        function Base.copy(
            tbd::TwoBodyPropagator{Tepoch, T}
        ) where {Tepoch <: Number, T <: Number}
            new_tbd = TwoBodyPropagator{Tepoch, T}()
            $(expressions...)
            return new_tbd
        end
    end
end

############################################################################################
#                                        Functions                                         #
############################################################################################

"""
    twobody_init(orb₀::KeplerianElements; kwargs...) -> TwoBodyPropagator

Create and initialize the two-body propagator structure using the mean Keplerian elements
`orb₀`.

!!! note

    The type used in the propagation will be the same as used to define the standard
    gravitational parameter `m0`.

# Keywords

- `m0::T`: Standard gravitational parameter of the central body [m³ / s²].
    (**Default** = `tbc_m0`)
"""
function twobody_init(
    orb₀::KeplerianElements{Tanomaly, Tepoch, Tkepler}; m0::T = tbc_m0
) where {
    Tanomaly <: AbstractAnomaly, Tepoch <: Number, Tkepler <: AbstractFloat, T <: Number
}
    # Allocate the propagator structure.
    tbd = TwoBodyPropagator{Tepoch, T}()

    # Assign the constant, which is used in the initialization.
    tbd.μ = m0

    # Initialize the propagator and return.
    twobody_init!(tbd, orb₀)

    return tbd
end

function twobody_init(
    orb₀::KeplerianElements{Tanomaly, Tepoch, Tkepler}; m0::Tm0 = tbc_m0
) where {Tanomaly <: AbstractAnomaly, Tepoch <: Number, Tkepler <: Number, Tm0 <: Number}
    T = promote_type(Tm0, Tkepler)

    # Allocate the propagator structure.
    tbd = TwoBodyPropagator{Tepoch, T}()

    # Assign the constant, which is used in the initialization.
    tbd.μ = T(m0)

    # Initialize the propagator and return.
    twobody_init!(tbd, orb₀)

    return tbd
end

"""
    twobody_init!(tbd::TwoBodyPropagator, orb₀::KeplerianElements) -> Nothing

Initialize the two-body propagator structure `tbd` using the mean Keplerian elements `orb₀`.

!!! warning

    The propagation constant `μ::Number` in `tbd` will not be changed. Hence, it must be
    initialized.
"""
function twobody_init!(
    tbd::TwoBodyPropagator{Tepoch, T}, orb₀::KeplerianElements
) where {Tepoch <: Number, T <: Number}
    # Make sure the Keplerian elements use the mean anomaly.
    ke₀ = convert(KeplerianElements{MeanAnomaly, Tepoch, T}, orb₀)

    # Unpack elements.
    a₀ = ke₀.semi_major_axis
    e₀ = ke₀.eccentricity
    M₀ = mean_anomaly(ke₀)

    # Without this check, the user would get a `DomainError` from an internal square root,
    # or silently wrong results, instead of a message pointing at the offending element.
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

    # Compute the mean motion using the semi-major axis.
    n₀ = √(tbd.μ / a₀^3)

    # Create and return the two-body orbit propagator structure.
    tbd.orb₀ = ke₀
    tbd.orbk = ke₀
    tbd.Δt   = 0
    tbd.n₀   = n₀

    return nothing
end

"""
    twobody(Δt::Number, orb₀::KeplerianElements; kwargs...) -> SVector{3, T}, SVector{3, T}, TwoBodyPropagator

Initialize the two-body propagator structure using the input elements `orb₀` and propagate
the orbit until the time Δt [s].

!!! note

    The type used in the propagation will be the same as used to define the standard
    gravitational parameter `m0`.

# Keywords

- `m0::T`: Standard gravitational parameter of the central body [m³ / s²].
    (**Default** = `tbc_m0`)

# Returns

- `SVector{3, T}`: Position vector [m] represented in the inertial frame at propagation
    instant.
- `SVector{3, T}`: Velocity vector [m / s] represented in the inertial frame at propagation
    instant.
- [`TwoBodyPropagator`](@ref): Structure with the initialized propagator.

# Remarks

The inertial frame in which the output is represented depends on which frame was used to
generate the orbit parameters.
"""
function twobody(Δt::Number, orb₀::KeplerianElements; m0::T = tbc_m0) where {T <: Number}
    tbd = twobody_init(orb₀; m0 = m0)
    r_i, v_i = twobody!(tbd, Δt)
    return r_i, v_i, tbd
end

"""
    twobody!(tbd::TwoBodyPropagator{Tepoch, T}, t::Number) where {Tepoch, T} -> SVector{3, T}, SVector{3, T}

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

The inertial frame in which the output is represented depends on which frame was used to
generate the orbit parameters.
"""
function twobody!(
    tbd::TwoBodyPropagator{Tepoch, T}, t::Number
) where {Tepoch <: Number, T <: Number}
    # Unpack.
    orb₀ = tbd.orb₀
    a₀   = orb₀.semi_major_axis
    e₀   = orb₀.eccentricity
    i₀   = orb₀.inclination
    Ω₀   = orb₀.raan
    ω₀   = orb₀.argument_of_periapsis
    M₀   = mean_anomaly(orb₀)

    # Time elapsed since epoch.
    epoch = orb₀.epoch
    Δt    = T(t)

    # Propagate the orbital elements.
    M_k = mod(M₀ + tbd.n₀ * Δt, T(2π))

    # Assemble the current mean elements.
    orbk = KeplerianElements{MeanAnomaly}(
        epoch + Tepoch(t) / 86400, a₀, e₀, i₀, Ω₀, ω₀, M_k
    )

    # Compute the position and velocity vectors given the orbital elements.
    r_i_k, v_i_k = kepler_to_rv(orbk)

    # Update the propagator structure.
    tbd.Δt   = Δt
    tbd.orbk = orbk

    # Return the position and velocity vector represented in the inertial
    # reference frame.
    return r_i_k, v_i_k
end
