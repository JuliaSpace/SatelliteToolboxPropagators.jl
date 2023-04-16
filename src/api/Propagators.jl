# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#
# Description
# ==========================================================================================
#
#    API for the orbit propagators.
#
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

module Propagators

using Dates
using Crayons

import Base: eltype, length, iterate, show

export OrbitPropagator

const _D = Crayon(reset = true)
const _B = crayon"bold"

"""
    abstract type OrbitPropagator{Tepoch<:Number, T<:Number}

Abstract type for the orbit propagators.
"""
abstract type OrbitPropagator{Tepoch<:Number, T<:Number} end

"""
    init(::Val{:propagator}, args...; kwargs...) -> OrbitPropagator

Create and initialize the orbit `propagator`. The arguments `args` and keywords `kwargs`
depends of the propagator type.
"""
function init end

"""
    init!(orbp::OrbitPropagator, args...; kwargs...) -> Nothing

Initialize the orbit propagator `orbp`. The arguments `args` and keywords `kwargs` depends
of the propagator type.
"""
function init! end

"""
    epoch(orbp::OrbitPropagator{Tepoch, T}) where {Tepoch<:Number, T<:Number} -> Tepoch

Return the initial elements' epoch of the propagator `orbp` [JD].
"""
function epoch end

"""
    last_instant(orbp::OrbitPropagator{Tepoch, T}) where {Tepoch<:Number, T<:Number} -> T

Return the last propagation instant [s] measured from the epoch. This action must be
performed by the function:
"""
function last_instant end

"""
    mean_elements(orbp::OrbitPropagator) -> Union{Nothing, KeplerianElements}

Return the mean elements using the structure [`KeplerianElements`](@ref) of the latest
propagation performed by `orbp`. This is an optinal function in the API. It will return
`nothing` if the propagator does not support it.
"""
mean_elements(orbp::OrbitPropagator) = nothing

"""
    name(orbp::OrbitPropagator) -> String

Return the name of the orbit propagator `orbp`. If this function is not defined, the
structure name is used: `typeof(orbp) |> string`.
"""
name(orbp::OrbitPropagator) = typeof(orbp) |> string

"""
    propagate(::Val{:propagator}, Δt::Number, args...; kwargs...)

Initialize the orbit `propagator` and propagate the orbit by `t` [s] from the initial orbit
epoch. The initialization arguments `args...` and `kwargs...` are the same as in the
initialization function [`Propagators.init`](@ref).

# Returns

- `SVector{3, T}`: Position vector [m] represented in the inertial frame at propagation
    instant.
- `SVector{3, T}`: Velocity vector [m / s] represented in the inertial frame at propagation
    instant.
- [`OrbitPropagator`](@ref): Structure with the initialized parameters.
"""
function propagate end

"""
    propagate!(orbp::OrbitPropagator{Tepoch, T}, t::Number) where {Tepoch, T} -> SVector{3, T}, SVector{3, T}

Propagate the orbit using `orbp` by `t` [s] from the initial orbit epoch.

# Returns

- `SVector{3, T}`: Position vector [m] represented in the inertial frame at propagation
    instant.
- `SVector{3, T}`: Velocity vector [m / s] represented in the inertial frame at propagation
    instant.
"""
function propagate! end

"""
    propagate_to_epoch(::Val{:propagator}, jd::Number, args...; kwargs...)

Initialize the orbit `propagator` and propagate the orbit until the epoch `jd` [s] from the
initial orbit epoch. The initialization arguments `args...` and `kwargs...` are the same as
in the initialization function [`Propagators.init`](@ref).

# Returns

- `SVector{3, T}`: Position vector [m] represented in the inertial frame at propagation
    instant.
- `SVector{3, T}`: Velocity vector [m / s] represented in the inertial frame at propagation
    instant.
- [`OrbitPropagator`](@ref): Structure with the initialized parameters.
"""
function propagate_to_epoch(T::Val, jd::Number, args...; kwargs...)
    orbp = init(T, args...; kwargs...)
    r_i, v_i = propagate_to_epoch!(orbp, jd)
    return r_i, v_i, orbp
end

"""
    propagate_to_epoch!(orbp::OrbitPropagator{Tepoch, T}, jd::Number) where {Tepoch, T} -> SVector{3, T}, SVector{3, T}

Propagate the orbit using `orbp` until the epoch `jd` [Julian Day].

# Returns

- `SVector{3, T}`: Position vector [m] represented in the inertial frame at propagation
    instant.
- `SVector{3, T}`: Velocity vector [m / s] represented in the inertial frame at propagation
    instant.
"""
function propagate_to_epoch!(orbp::OrbitPropagator, jd::Number)
    return propagate!(orbp, 86400 * (jd - epoch(orbp)))
end

"""
    step!(orbp::OrbitPropagator{Tepoch, T}, Δt::Number) where {Tepoch, T} -> SVector{3, T}, SVector{3, T}

Propagate the orbit using `orbp` by `Δt` [s] from the current orbit epoch.

# Returns

- `SVector{3, T}`: Position vector [m] represented in the inertial frame at propagation
    instant.
- `SVector{3, T}`: Velocity vector [m / s] represented in the inertial frame at propagation
    instant.
"""
function step!(orbp::OrbitPropagator, Δt::Number)
    return propagate!(orbp, last_instant(orbp) + Δt)
end

############################################################################################
#                                    Iterator Interface
############################################################################################

# There functions allow broadcast when using the orbit propagators.
iterate(orbp::OrbitPropagator) = (orbp, nothing)
iterate(orbp::OrbitPropagator, ::Nothing) = nothing
length(orbp::OrbitPropagator) = 1
eltype(orbp::T) where T<:OrbitPropagator = T

############################################################################################
#                                           Show
############################################################################################

function show(io::IO, orbp::T) where T<:OrbitPropagator
    prop_epoch = epoch(orbp) |> julian2datetime
    Δt         = last_instant(orbp)
    print(io, name(orbp), " (Epoch = ", prop_epoch, ", Δt = ", Δt, " s)")
    return nothing
end

function show(io::IO, mime::MIME"text/plain", orbp::T) where T<:OrbitPropagator
    # Check for color support in the `io`.
    color = get(io, :color, false)
    b = color ? string(_B) : ""
    d = color ? string(_D) : ""

    prop_name       = name(orbp)
    prop_epoch      = epoch(orbp)
    prop_epoch_dt   = prop_epoch |> julian2datetime
    last_instant_dt = prop_epoch + last_instant(orbp) / 86400 |> julian2datetime

    println(io, string(T), ":")
    println(io, "$(b)   Propagator name :$(d) ", prop_name)
    println(io, "$(b)  Propagator epoch :$(d) ", prop_epoch_dt)
    print(  io, "$(b)  Last propagation :$(d) ", last_instant_dt)

    return nothing
end

end # module Propagator
