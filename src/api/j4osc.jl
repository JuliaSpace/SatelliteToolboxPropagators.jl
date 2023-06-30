# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#
# Description
# ==========================================================================================
#
#    API implementation for J4 osculating orbit propagator.
#
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

Propagators.epoch(orbp::OrbitPropagatorJ4Osculating)         = orbp.j4oscd.j4d.orb₀.t
Propagators.last_instant(orbp::OrbitPropagatorJ4Osculating)  = orbp.j4oscd.Δt
Propagators.mean_elements(orbp::OrbitPropagatorJ4Osculating) = orbp.j4oscd.j4d.orbk
Propagators.name(orbp::OrbitPropagatorJ4Osculating)          = "J4 Osculating Orbit Propagator"

"""
    Propagators.init(Val(:J4osc), orb₀::KeplerianElements, dn_o2::Number = 0, ddn_o6::Number = 0; kwargs...) -> OrbitPropagatorJ4Osculating

Create and initialize the J4 osculating orbit propagator structure using the mean Keplerian
elements `orb₀`.

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
  [`J4PropagatorConstants`](@ref)). (**Default** = `j4c_egm2008`)
"""
function Propagators.init(
    ::Val{:J4osc},
    orb₀::KeplerianElements,
    dn_o2::Number = 0,
    ddn_o6::Number = 0;
    j4c::J4PropagatorConstants = j4c_egm2008
)
    j4oscd = j4osc_init(orb₀, dn_o2, ddn_o6; j4c = j4c)
    return OrbitPropagatorJ4Osculating(j4oscd)
end

"""
    Propagators.init!(orbp::OrbitPropagatorJ4Osculating, orb₀::KeplerianElements, dn_o2::Number = 0, ddn_o6::Number = 0) -> Nothing

Initialize the J4 osculating orbit propagator structure `orbp` using the mean Keplerian
elements `orb₀`.

!!! warning
    The propagation constants `j4c::J4PropagatorConstants` in `orbp.j4d` will not be
    changed. Hence, they must be initialized.

# Arguments

- `orb₀::KeplerianElements`: Initial mean Keplerian elements [SI units].
- `dn_o2::Number`: First time derivative of the mean motion divided by two [rad/s^2].
    (**Default** = 0)
- `ddn_o6::Number`: Second time derivative of the mean motion divided by six [rad/s^3].
    (**Default** = 0)
"""
function Propagators.init!(
    orbp::OrbitPropagatorJ4Osculating,
    orb₀::KeplerianElements,
    dn_o2::Number = 0,
    ddn_o6::Number = 0
)
    j4osc_init!(orbp.j4oscd, orb₀, dn_o2, ddn_o6)
    return nothing
end

"""
    Propagators.propagate(Val(:J4osc), Δt::Number, orb₀::KeplerianElements, dn_o2::Number = 0, ddn_o6::Number = 0; kwargs...)

Initialize the J4 osculating propagator structure using the input elements `orb₀` and
propagate the orbit until the time Δt [s].

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

- `j4c::J4PropagatorConstants{T}`: J4 orbit propagator constants (see
  [`J4PropagatorConstants`](@ref)). (**Default** = `j4c_egm2008`)

# Returns

- `SVector{3, T}`: Position vector [m] represented in the inertial frame at propagation
    instant.
- `SVector{3, T}`: Velocity vector [m / s] represented in the inertial frame at propagation
    instant.
- [`OrbitPropagatorJ4Osculating`](@ref): Structure with the initialized parameters.
"""
function Propagators.propagate(
    ::Val{:J4osc},
    Δt::Number,
    orb₀::KeplerianElements,
    dn_o2::Number = 0,
    ddn_o6::Number = 0;
    j4c::J4PropagatorConstants = j4c_egm2008
)
    r_i, v_i, j4oscd = j4osc(Δt, orb₀, dn_o2, ddn_o6; j4c = j4c)
    return r_i, v_i, OrbitPropagatorJ4Osculating(j4oscd)
end

function Propagators.propagate!(orbp::OrbitPropagatorJ4Osculating, t::Number)
    # Auxiliary variables.
    j4oscd = orbp.j4oscd

    # Propagate the orbit.
    return j4osc!(j4oscd, t)
end
