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
    init_orbit_propagator(Val(:J4), orb₀::Orbit, dn_o4::Number = 0, ddn_o6::Number = 0; kwargs...) -> OrbitPropagatorJ4

Initialize the J4 orbit propagator algorithm using the mean Keplerian elements `orb₀`.

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
    # Create the new Two Body propagator structure.
    j4d = j4_init(orb₀, dn_o4, ddn_o6; j4c = j4c)

    # Create and return the orbit propagator structure.
    return OrbitPropagatorJ4(j4d)
end

function Propagators.propagate!(orbp::OrbitPropagatorJ4, t::Number)
    # Auxiliary variables.
    j4d = orbp.j4d

    # Propagate the orbit.
    return j4!(j4d, t)
end
