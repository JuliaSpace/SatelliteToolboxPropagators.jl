# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#
# Description
# ==========================================================================================
#
#    API for the orbit propagators.
#
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

export init_orbit_propagator, get_epoch, get_mean_elements
export propagate!, propagate_to_epoch!, step!

"""
    init_orbit_propagator(T, args...; kwargs...)

Initialize the orbit propagator of type `T`. The arguments `args` and keywords `kwargs`
depends of the propagator type.
"""
init_orbit_propagator

"""
    get_epoch(orbp)

Return the epoch of the propagator `orbp` [JD].
"""
get_epoch

"""
    get_mean_elements(orbp)

Return the mean elements of the latest propagation performed by `orbp`. This is an optinal
function in the API. It will return `nothing` if the propagator does not support it.
"""
get_mean_elements(orbp::OrbitPropagator) = nothing

"""
    propagate!(orbp::OrbitPropagator{Tepoch, T}, t::Number) where {Tepoch, T} -> SVector{3, T}, SVector{3, T}

Propagate the orbit using `orbp` by `t` [s] from the initial orbit epoch.

# Returns

- `SVector{3, T}`: Position vector [m] represented in the inertial frame at propagation
    instant.
- `SVector{3, T}`: Velocity vector [m / s] represented in the inertial frame at propagation
    instant.
"""
propagate!

"""
    propagate_to_epoch!(orbp::OrbitPropagator{Tepoch, T}, t::Number) where {Tepoch, T} -> SVector{3, T}, SVector{3, T}

Propagate the orbit using `orbp` until the epoch `JD` [Julian Day].

# Returns

- `SVector{3, T}`: Position vector [m] represented in the inertial frame at propagation
    instant.
- `SVector{3, T}`: Velocity vector [m / s] represented in the inertial frame at propagation
    instant.
"""
function propagate_to_epoch!(orbp::OrbitPropagator, jd::Number)
    return propagate!(orbp, 86400 * (jd - get_epoch(orbp)))
end

"""
    step!(orbp::OrbitPropagator{Tepoch, T}, t::Number) where {Tepoch, T} -> SVector{3, T}, SVector{3, T}

Propagate the orbit using `orbp` by `Δt` [s] from the current orbit epoch.

# Returns

- `SVector{3, T}`: Position vector [m] represented in the inertial frame at propagation
    instant.
- `SVector{3, T}`: Velocity vector [m / s] represented in the inertial frame at propagation
    instant.
"""
step!

############################################################################################
#                                    Iterator Interface
############################################################################################

# There functions allow broadcast when using the orbit propagators.
Base.iterate(orbp::OrbitPropagator) = (orbp, nothing)
Base.iterate(orbp::OrbitPropagator, ::Nothing) = nothing
Base.length(orbp::OrbitPropagator) = 1
Base.eltype(orbp::T) where T<:OrbitPropagator = T
