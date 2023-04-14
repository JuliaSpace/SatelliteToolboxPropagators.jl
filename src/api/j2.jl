# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#
# Description
# ==========================================================================================
#
#    API implementation for J2 orbit propagator.
#
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

Propagators.epoch(orbp::OrbitPropagatorJ2)         = orbp.j2d.orb₀.t
Propagators.last_instant(orbp::OrbitPropagatorJ2)  = orbp.j2d.Δt
Propagators.mean_elements(orbp::OrbitPropagatorJ2) = orbp.j2d.orbk
Propagators.name(orbp::OrbitPropagatorJ2)          = "J2 Orbit Propagator"

"""
    Propagators.init(Val(:J2), orb₀::KeplerianElements, dn_o2::Number = 0, ddn_o6::Number = 0; kwargs...) -> OrbitPropagatorJ2

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
function Propagators.init(
    ::Val{:J2},
    orb₀::KeplerianElements,
    dn_o2::Number = 0,
    ddn_o6::Number = 0;
    j2c::J2PropagatorConstants = j2c_egm08
)
    j2d = j2_init(orb₀, dn_o2, ddn_o6; j2c = j2c)
    return OrbitPropagatorJ2(j2d)
end

"""
    Propagators.init!(orbp::OrbitPropagatorJ2, orb₀::KeplerianElements, dn_o2::Number = 0, ddn_o6::Number = 0) -> Nothing

Initialize the J2 orbit propagator structure `orbp` using the mean Keplerian elements
`orb₀`.

!!! warning
    The propagation constants `j2c::J2PropagatorConstants` in `orbp.j2d` will not be
    changed. Hence, they must be initialized.

# Arguments

- `orb₀::KeplerianElements`: Initial mean Keplerian elements [SI units].
- `dn_o2::Number`: First time derivative of the mean motion divided by two [rad/s^2].
    (**Default** = 0)
- `ddn_o6::Number`: Second time derivative of the mean motion divided by six [rad/s^3].
    (**Default** = 0)
"""
function Propagators.init!(
    orbp::OrbitPropagatorJ2,
    orb₀::KeplerianElements,
    dn_o2::Number = 0,
    ddn_o6::Number = 0
)
    j2_init!(orbp.j2d, orb₀, dn_o2, ddn_o6)
    return nothing
end

function Propagators.propagate!(orbp::OrbitPropagatorJ2, t::Number)
    # Auxiliary variables.
    j2d = orbp.j2d

    # Propagate the orbit.
    return j2!(j2d, t)
end
