# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#
# Description
# ==========================================================================================
#
#    API implementation for J2 osculating orbit propagator.
#
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

Propagators.epoch(orbp::OrbitPropagatorJ2Osculating)         = orbp.j2oscd.j2d.orb₀.t
Propagators.last_instant(orbp::OrbitPropagatorJ2Osculating)  = orbp.j2oscd.Δt
Propagators.mean_elements(orbp::OrbitPropagatorJ2Osculating) = orbp.j2oscd.j2d.orbk
Propagators.name(orbp::OrbitPropagatorJ2Osculating)          = "J2 Osculating Orbit Propagator"

"""
    Propagators.init(Val(:J2osc), orb₀::KeplerianElements, dn_o2::Number = 0, ddn_o6::Number = 0; kwargs...) -> OrbitPropagatorJ2Osculating

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
  [`J2PropagatorConstants`](@ref)). (**Default** = `j2c_egm2008`)
"""
function Propagators.init(
    ::Val{:J2osc},
    orb₀::KeplerianElements,
    dn_o2::Number = 0,
    ddn_o6::Number = 0;
    j2c::J2PropagatorConstants = j2c_egm2008
)
    j2oscd = j2osc_init(orb₀, dn_o2, ddn_o6; j2c = j2c)
    return OrbitPropagatorJ2Osculating(j2oscd)
end

"""
    Propagators.init!(orbp::OrbitPropagatorJ2Osculating, orb₀::KeplerianElements, dn_o2::Number = 0, ddn_o6::Number = 0) -> Nothing

Initialize the J2 osculating orbit propagator structure `orbp` using the mean Keplerian
elements `orb₀`.

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
    orbp::OrbitPropagatorJ2Osculating,
    orb₀::KeplerianElements,
    dn_o2::Number = 0,
    ddn_o6::Number = 0
)
    j2osc_init!(orbp.j2oscd, orb₀, dn_o2, ddn_o6)
    return nothing
end

"""
    Propagators.propagate(Val(:J2osc), Δt::Number, orb₀::KeplerianElements, dn_o2::Number = 0, ddn_o6::Number = 0; kwargs...)

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
  [`J2PropagatorConstants`](@ref)). (**Default** = `j2c_egm2008`)

# Returns

- `SVector{3, T}`: Position vector [m] represented in the inertial frame at propagation
    instant.
- `SVector{3, T}`: Velocity vector [m / s] represented in the inertial frame at propagation
    instant.
- [`OrbitPropagatorJ2Osculating`](@ref): Structure with the initialized parameters.
"""
function Propagators.propagate(
    ::Val{:J2osc},
    Δt::Number,
    orb₀::KeplerianElements,
    dn_o2::Number = 0,
    ddn_o6::Number = 0;
    j2c::J2PropagatorConstants = j2c_egm2008
)
    r_i, v_i, j2oscd = j2osc(Δt, orb₀, dn_o2, ddn_o6; j2c = j2c)
    return r_i, v_i, OrbitPropagatorJ2Osculating(j2oscd)
end

function Propagators.propagate!(orbp::OrbitPropagatorJ2Osculating, t::Number)
    # Auxiliary variables.
    j2oscd = orbp.j2oscd

    # Propagate the orbit.
    return j2osc!(j2oscd, t)
end
