# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#
# Description
# ==========================================================================================
#
#    API implementation for J4 orbit propagator.
#
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

Propagators.epoch(orbp::OrbitPropagatorJ4)         = orbp.j4d.orb₀.t
Propagators.last_instant(orbp::OrbitPropagatorJ4)  = orbp.j4d.Δt
Propagators.mean_elements(orbp::OrbitPropagatorJ4) = orbp.j4d.orbk
Propagators.name(orbp::OrbitPropagatorJ4)          = "J4 Orbit Propagator"

"""
    Propagators.init(Val(:J4), orb₀::Orbit, dn_o4::Number = 0, ddn_o6::Number = 0; kwargs...) -> OrbitPropagatorJ4

Create and initialize the J4 orbit propagator structure using the mean Keplerian elements
`orb₀`.

!!! note
    The type used in the propagation will be the same as used to define the constants in the
    structure `j4c`.

# Arguments

- `orb₀::KeplerianElements`: Initial mean Keplerian elements [SI units].
- `dn_o4::Number`: First time derivative of the mean motion divided by two [rad/s^4].
    (**Default** = 0)
- `ddn_o6::Number`: Second time derivative of the mean motion divided by six [rad/s^3].
    (**Default** = 0)

# Keywords

- `j4c::J4PropagatorConstants`: J4 orbit propagator constants (see
  [`J4PropagatorConstants`](@ref)). (**Default** = `j4c_egm08`)
"""
function Propagators.init(
    ::Val{:J4},
    orb₀::KeplerianElements,
    dn_o4::Number = 0,
    ddn_o6::Number = 0;
    j4c::J4PropagatorConstants = j4c_egm08
)
    j4d = j4_init(orb₀, dn_o4, ddn_o6; j4c = j4c)
    return OrbitPropagatorJ4(j4d)
end

"""
    Propagators.init!(orbp::OrbitPropagatorJ4, orb₀::KeplerianElements, dn_o2::Number = 0, ddn_o6::Number = 0) -> Nothing

Initialize the J4 orbit propagator structure `orbp` using the mean Keplerian elements
`orb₀`.

!!! warning
    The propagation constants `j4c::J4PropagatorConstants` in `j4d` will not be changed.
    Hence, they must be initialized.

# Arguments

- `orb₀::KeplerianElements`: Initial mean Keplerian elements [SI units].
- `dn_o2::Number`: First time derivative of the mean motion divided by two [rad/s^2].
    (**Default** = 0)
- `ddn_o6::Number`: Second time derivative of the mean motion divided by six [rad/s^3].
    (**Default** = 0)
"""
function Propagators.init!(
    orbp::OrbitPropagatorJ4,
    orb₀::KeplerianElements,
    dn_o2::Number = 0,
    ddn_o6::Number = 0
)
    j4_init!(orbp.j4d, orb₀, dn_o2, ddn_o6)
    return nothing
end

function Propagators.propagate!(orbp::OrbitPropagatorJ4, t::Number)
    # Auxiliary variables.
    j4d = orbp.j4d

    # Propagate the orbit.
    return j4!(j4d, t)
end
