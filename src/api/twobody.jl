# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#
# Description
# ==========================================================================================
#
#    API implementation for two body orbit propagator.
#
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

Propagators.epoch(orbp::OrbitPropagatorTwoBody)         = orbp.tbd.orb₀.t
Propagators.last_instant(orbp::OrbitPropagatorTwoBody)  = orbp.tbd.Δt
Propagators.mean_elements(orbp::OrbitPropagatorTwoBody) = orbp.tbd.orbk
Propagators.name(orbp::OrbitPropagatorTwoBody)          = "Two Body Orbit Propagator"

"""
    Propagators.init(Val(:TwoBody), orb₀::KeplerianElements; kwargs...) -> OrbitPropagatorTwoBody

Create and initialize the two body orbit propagator structure using the mean Keplerian
elements `orb₀`.

!!! note
    The type used in the propagation will be the same as used to define the gravitational
    constant `μ`.

# Arguments

- `orb₀::KeplerianElements`: Initial mean Keplerian elements [SI units].

# Keywords

- `μ::T`: Standard gravitational parameter of the central body [m³/s²].
    (**Default** = `tbc_m0`)
"""
function Propagators.init(::Val{:TwoBody}, orb₀::KeplerianElements; μ::Number = tbc_m0)
    tbd = twobody_init(orb₀; μ = μ)
    return OrbitPropagatorTwoBody(tbd)
end

"""
    Propagators.init!(orbp::OrbitPropagatorTwoBody, orb₀::KeplerianElements) -> Nothing

Initialize the two body orbit propagator structure `orbp` using the mean Keplerian elements
`orb₀`.

!!! warning
    The propagation constant `μ::Number` in `tbd` will not be changed. Hence, it must be
    initialized.

# Arguments

- `orb₀::KeplerianElements`: Initial mean Keplerian elements [SI units].
"""
function Propagators.init!(orbp::OrbitPropagatorTwoBody, orb₀::KeplerianElements)
    twobody_init!(orbp.tbd, orb₀)
    return nothing
end

"""
    Propagators.propagate(Val(:TwoBody), Δt::Number, orb₀::KeplerianElements; kwargs...)

Initialize the two body propagator structure using the input elements `orb₀` and propagate
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
"""
function Propagators.propagate(
    ::Val{:TwoBody},
    Δt::Number,
    orb₀::KeplerianElements;
    μ::Number = tbc_m0
)
    r_i, v_i, tbd = twobody(Δt, orb₀; μ = μ)
    return r_i, v_i, OrbitPropagatorTwoBody(tbd)
end

function Propagators.propagate!(orbp::OrbitPropagatorTwoBody, t::Number)
    # Auxiliary variables.
    tbd = orbp.tbd

    # Propagate the orbit.
    return twobody!(tbd, t)
end
