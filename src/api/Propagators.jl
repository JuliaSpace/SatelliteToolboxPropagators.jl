## Description #############################################################################
#
#  API for the orbit propagators.
#
############################################################################################

module Propagators

using Dates
using Crayons
using StaticArrays

import Base: copy, eltype, length, iterate, show
import SatelliteToolboxBase: @maybe_threads, get_partition

export OrbitPropagator

############################################################################################
#                                         Contants                                         #
############################################################################################

# Escape sequences related to the crayons.
const _D = string(Crayon(reset = true))
const _B = string(crayon"bold")

############################################################################################
#                                          Types                                           #
############################################################################################

"""
    abstract type OrbitPropagator{Tepoch<:Number, T<:Number}

Abstract type for the orbit propagators.
"""
abstract type OrbitPropagator{Tepoch<:Number, T<:Number} end

############################################################################################
#                                     Public Functions                                     #
############################################################################################

"""
    fit_mean_elements(::Val{:propagator}, vjd::AbstractVector{Tjd}, vr_i::AbstractVector{Tv}, vv_i::AbstractVector{Tv}; kwargs...) where {Tjd<:Number, Tv<:AbstractVector} -> <Mean elements>

Fit a set of mean elements for the `propagator` using the osculating elements represented by
a set of position vectors `vr_i` [m] and a set of velocity vectors `vv_i` [m / s]
represented in an inertial reference frame at instants in the array `vjd` [Julian Day]. The
keywords `kwargs` depends on the propagator type.

This function returns the set of mean elements used to initialize the `propagator`.
"""
function fit_mean_elements end

"""
    fit_mean_elements!(orbp::OrbitPropagator, vjd::AbstractVector{Tjd}, vr_i::AbstractVector{Tv}, vv_i::AbstractVector{Tv}; kwargs...) where {Tjd<:Number, Tv<:AbstractVector} -> <Mean elements>

Fit a set of mean elements for the propagator `orbp` using the osculating elements
represented by a set of position vectors `vr_i` [m] and a set of velocity vectors `vv_i` [m
/ s] represented in an inertial reference frame at instants in the array `vjd` [Julian Day].
The keywords `kwargs` depends on the propagator type.

This function returns the set of mean elements used to initialize the `propagator` and also
initializes `orbp` with the fitted mean elements.
"""
function fit_mean_elements! end

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

Return the last propagation instant [s] measured from the epoch.
"""
function last_instant end

"""
    mean_elements(orbp::OrbitPropagator) -> Union{Nothing, KeplerianElements}

Return the mean elements using the structure `KeplerianElements` of the latest propagation
performed by `orbp`. This is an optinal function in the API. It will return `nothing` if the
propagator does not support it.
"""
mean_elements(orbp::OrbitPropagator) = nothing

"""
    name(orbp::OrbitPropagator) -> String

Return the name of the orbit propagator `orbp`. If this function is not defined, the
structure name is used: `typeof(orbp) |> string`.
"""
name(orbp::OrbitPropagator) = typeof(orbp) |> string

"""
    propagate(::Val{:propagator}, Δt::Number, args...; kwargs...) -> SVector{3, T}, SVector{3, T}, OrbitPropagator{Tepoch, T}

Initialize the orbit `propagator` and propagate the orbit by `t` [s] from the initial orbit
epoch. The initialization arguments `args...` and `kwargs...` are the same as in the
initialization function [`Propagators.init`](@ref).

!!! note

    `T` is the propagator number type. For more information, see [`Propagators.init`](@ref).

# Returns

- `SVector{3, T}`: Position vector [m] represented in the inertial frame at propagation
    instant.
- `SVector{3, T}`: Velocity vector [m / s] represented in the inertial frame at propagation
    instant.
- [`OrbitPropagator{Tepoch, T}`](@ref): Structure with the initialized propagator.
"""
function propagate(prop, Δt::Number, args...; kwargs...)
    orbp = Propagators.init(prop, args...; kwargs...)
    r_i, v_i = Propagators.propagate!(orbp, Δt)
    return r_i, v_i, orbp
end

"""
    propagate(::Val{:propagator}, vt::AbstractVector, args...; kwargs...) -> Vector{SVector{3, T}}, Vector{SVector{3, T}}, OrbitPropagator{Tepoch, T}

Initialize the orbit `propagator` and propagate the orbit for every instant defined in `vt`
[s] from the initial orbit epoch. The initialization arguments `args...` and `kwargs...`
(except for `ntasks`) are the same as in the initialization function
[`Propagators.init`](@ref).

!!! note

    `T` is the propagator number type. For more information, see [`Propagators.init`](@ref).

# Keywords

- `ntasks::Integer`: Number of parallel tasks to propagate the orbit. If it is set to a
    number equal or lower than 1, the function will propagate the orbit sequentially.
    (**Default** = `Threads.nthreads()`)

# Returns

- `Vector{SVector{3, T}}`: Array with the position vectors [m] in the inertial frame at each
    propagation instant defined in `vt`.
- `Vector{SVector{3, T}}`: Array with the velocity vectors [m / s] in the inertial frame at
    each propagation instant defined in `vt`.
- [`OrbitPropagator{Tepoch, T}`](@ref): Structure with the initialized propagator.
"""
function propagate(prop, vt::AbstractVector, args...; kwargs...)
    orbp = Propagators.init(prop, args...; kwargs...)
    vr_i, vv_i = Propagators.propagate!(orbp, vt)
    return vr_i, vv_i, orbp
end

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
    propagate!(orbp::OrbitPropagator{Tepoch, T}, vt::AbstractVector; kwargs...) where {Tepoch <: Number, T <: Number} -> Vector{SVector{3, T}}, Vector{SVector{3, T}}

Propagate the orbit using `orbp` for every instant defined in `vt` [s].

# Keywords

- `ntasks::Integer`: Number of parallel tasks to propagate the orbit. If it is set to a
    number equal or lower than 1, the function will propagate the orbit sequentially.
    (**Default** = `Threads.nthreads()`)

# Returns

- `Vector{SVector{3, T}}`: Array with the position vectors [m] in the inertial frame at each
    propagation instant defined in `vt`.
- `Vector{SVector{3, T}}`: Array with the velocity vectors [m / s] in the inertial frame at
    each propagation instant defined in `vt`.
"""
function propagate!(
    orbp::OrbitPropagator{Tepoch, T},
    vt::AbstractVector;
    ntasks::Integer = Threads.nthreads()
) where {Tepoch<:Number, T<:Number}
    # We need to perform the first propagation to obtain the return type of the propagator.
    r₀, v₀ = Propagators.propagate!(orbp, first(vt))

    # Number of propagation points.
    len_vt = length(vt)

    # Allocate the output vectors.
    vr = Vector{typeof(r₀)}(undef, len_vt)
    vv = Vector{typeof(v₀)}(undef, len_vt)

    vr[begin] = r₀
    vv[begin] = v₀

    len_vt == 1 && return vr, vv

    inds = eachindex(vt)

    # We need to store the first index offset of `vt` to allow filling the output vectors
    # correctly.
    Δi = firstindex(vt) - 1

    # Make sure the number of tasks is not higher than the number of propagation points.
    ntasks = min(ntasks, len_vt)

    # If we have only two instants in the time vector, we will not spawn any threads,
    # because the first and the last instants are propagated separately.
    if len_vt > 2
        @maybe_threads ntasks for c in 1:ntasks
            # We already propagated for the first instant, and we must ensure we propagate
            # the last instant at the end of the function.
            i₀, i₁ = @views get_partition(c, inds[(1 + begin):(end - 1)], ntasks)

            # The propagation usually modifies the structure. Hence we need to copy it for each
            # task.
            corbp = c == 1 ? orbp : copy(orbp)

            @inbounds for i in i₀:i₁
                vr[i - Δi], vv[i - Δi] = Propagators.propagate!(corbp, vt[i])
            end
        end
    end

    # We must ensure that the last propagation instant is the obtained at the end to keep
    # the internal data of the propagation consistent.
    vr[end], vv[end] = Propagators.propagate!(orbp, last(vt))

    return vr, vv
end

"""
    propagate_to_epoch(::Val{:propagator}, jd::Number, args...; kwargs...)
    propagate_to_epoch(::Val{:propagator}, dt::DateTime, args...; kwargs...)

Initialize the orbit `propagator` and propagate the orbit until the epoch defined by either
the Julian Day `jd` [UTC] or by a `DateTime` object `dt` [UTC] from the initial orbit epoch.
The initialization arguments `args...` and `kwargs...` are the same as in the initialization
function [`Propagators.init`](@ref).

# Returns

- `SVector{3, T}`: Position vector [m] represented in the inertial frame at propagation
    instant.
- `SVector{3, T}`: Velocity vector [m / s] represented in the inertial frame at propagation
    instant.
- [`OrbitPropagator`](@ref): Structure with the initialized parameters.
"""
function propagate_to_epoch(T::Val, dt::DateTime, args...; kwargs...)
    jd = datetime2julian(dt)
    return propagate_to_epoch(T, jd, args...; kwargs...)
end

function propagate_to_epoch(T::Val, jd::Number, args...; kwargs...)
    orbp = init(T, args...; kwargs...)
    r_i, v_i = propagate_to_epoch!(orbp, jd)
    return r_i, v_i, orbp
end

"""
    propagate_to_epoch!(orbp::OrbitPropagator{Tepoch, T}, jd::Number) where {Tepoch, T} -> SVector{3, T}, SVector{3, T}
    propagate_to_epoch!(orbp::OrbitPropagator{Tepoch, T}, dt::DateTime) where {Tepoch, T} -> SVector{3, T}, SVector{3, T}

Propagate the orbit using `orbp` until the epoch defined either by the Julian Day `jd`
[UTC] or by the `DateTime` object `dt` [UTC].

# Returns

- `SVector{3, T}`: Position vector [m] represented in the inertial frame at propagation
    instant.
- `SVector{3, T}`: Velocity vector [m / s] represented in the inertial frame at propagation
    instant.
"""
function propagate_to_epoch!(orbp::OrbitPropagator, dt::DateTime)
    jd = datetime2julian(dt)
    return propagate_to_epoch!(orbp, jd)
end

function propagate_to_epoch!(orbp::OrbitPropagator, jd::Number)
    return propagate!(orbp, 86400 * (jd - epoch(orbp)))
end

"""
    propagate_to_epoch!(orbp::OrbitPropagator{Tepoch, T}, vjd::AbstractVector; kwargs...) where {Tepoch, T} -> SVector{3, T}, SVector{3, T}
    propagate_to_epoch!(orbp::OrbitPropagator{Tepoch, T}, vdt::AbstractVector{DateTime}; kwargs...) where {Tepoch, T} -> SVector{3, T}, SVector{3, T}

Propagate the orbit using `orbp` for every epoch defined in the vector of Julian Days `vjd`
[UTC] or in the vector of `DateTime` objects `vdt` [UTC].

# Keywords

- `ntasks::Integer`: Number of parallel tasks to propagate the orbit. If it is set to a
    number equal or lower than 1, the function will propagate the orbit sequentially.
    (**Default** = `Threads.nthreads()`)

# Returns

- `Vector{SVector{3, T}}`: Array with the position vectors [m] in the inertial frame at each
    propagation instant defined in `vt`.
- `Vector{SVector{3, T}}`: Array with the velocity vectors [m / s] in the inertial frame at
    each propagation instant defined in `vt`.
"""
function propagate_to_epoch!(orbp::OrbitPropagator, vdt::AbstractVector{T}) where T<:DateTime
    vjd = datetime2julian.(vdt)
    return propagate_to_epoch!(orbp, vjd)
end

function propagate_to_epoch!(orbp::OrbitPropagator, vjd::AbstractVector)
    return propagate!(orbp, 86400 .* (vjd .- epoch(orbp)))
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
#                                           Copy                                           #
############################################################################################

# Define a fallback copy method if the propagator has not implemented one.
Base.copy(orbp::OrbitPropagator) = deepcopy(orbp)

############################################################################################
#                                    Iterator Interface                                    #
############################################################################################

# There functions allow broadcast when using the orbit propagators.
iterate(orbp::OrbitPropagator) = (orbp, nothing)
iterate(orbp::OrbitPropagator, ::Nothing) = nothing
length(orbp::OrbitPropagator) = 1
eltype(orbp::T) where T<:OrbitPropagator = T

############################################################################################
#                                           Show                                           #
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
    b = color ? _B : ""
    d = color ? _D : ""

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
