module Propagators

using Dates

import Base: copy, eltype, length, iterate, show
import SatelliteToolboxBase
import SatelliteToolboxBase: @maybe_threads, get_partition, OrbitStateVector

export OrbitPropagator

############################################################################################
#                                          Types                                           #
############################################################################################

"""
    abstract type OrbitPropagator{Tepoch <: Number, T <: Number}

Abstract type for the orbit propagators, where `Tepoch` is the type used to represent the
epoch of the input elements and `T` is the type used for the internal variables.
"""
abstract type OrbitPropagator{Tepoch <: Number, T <: Number} end

# Types accepted as a propagation instant measured from the epoch: a number [s] or a period.
const PropagationInstant = Union{Number, Dates.Period, Dates.CompoundPeriod}

# Types accepted as an epoch: a Julian Day [UTC] or a `DateTime` [UTC].
const PropagationEpoch = Union{Number, DateTime}

# Types accepted as the output sink of the propagation functions.
const PropagationSink = Union{Type{Tuple}, Type{OrbitStateVector}}

############################################################################################
#                                     Public Functions                                     #
############################################################################################

"""
    fit_mean_elements(
        [sink::Type, ]::Val{:propagator},
        vjd::AbstractVector{Tjd},
        vr_i::AbstractVector{Tv},
        vv_i::AbstractVector{Tv};
        kwargs...
    ) where {Tjd <: Number, Tv <: AbstractVector} -> <Mean elements>, <Covariance>

    fit_mean_elements(
        [sink::Type, ]::Val{:propagator},
        vsv::AbstractVector{OrbitStateVector{Tepoch, T}};
        kwargs...
    ) where {Tepoch <: Number, T <: Number} -> <Mean elements>, <Covariance>

Fit a set of mean elements for the `propagator` using the osculating state vector
represented in an inertial reference frame. The state vector can be represented using a set
of position vectors `vr_i` [m] and a set of velocity vectors `vv_i` [m / s] obtained at the
instants in the array `vjd` [Julian Day], or an array of `OrbitStateVector` `vsv` [SI],
containing the same information. The keywords `kwargs` depend on the propagator type.

The optional parameter `sink` selects the representation of the mean elements for the
propagators that support more than one. SGP4 accepts `TLE` and `OrbitMeanElementsMessage`,
defaulting to the latter.

# Returns

- `<Mean elements>`: Set of mean elements used to initialize the `propagator`. The concrete
    type depends on the propagator, and it is `KeplerianElements{MeanAnomaly}` for every
    propagator except SGP4, which returns an `OrbitMeanElementsMessage` or a `TLE` selected
    by `sink`.
- `<Covariance>`: Final covariance matrix of the least-square algorithm. It is a
    `SMatrix{6, 6}` for every propagator except SGP4, which returns a `SMatrix{7, 7}`
    because it also fits the B* parameter.
"""
function fit_mean_elements end

function fit_mean_elements(
    prop::Val, vsv::AbstractVector{OrbitStateVector{Tepoch, T}}; kwargs...
) where {Tepoch <: Number, T <: Number}
    vjd  = map(x -> x.epoch, vsv)
    vr_i = map(x -> x.r, vsv)
    vv_i = map(x -> x.v, vsv)

    return fit_mean_elements(prop, vjd, vr_i, vv_i; kwargs...)
end

function fit_mean_elements(
    sink::Type, prop::Val, vsv::AbstractVector{OrbitStateVector{Tepoch, T}}; kwargs...
) where {Tepoch <: Number, T <: Number}
    vjd  = map(x -> x.epoch, vsv)
    vr_i = map(x -> x.r, vsv)
    vv_i = map(x -> x.v, vsv)

    return fit_mean_elements(sink, prop, vjd, vr_i, vv_i; kwargs...)
end

"""
    fit_mean_elements!(
        orbp::OrbitPropagator,
        vjd::AbstractVector{Tjd},
        vr_i::AbstractVector{Tv},
        vv_i::AbstractVector{Tv}[, sink::Type];
        kwargs...
    ) where {Tjd <: Number, Tv <: AbstractVector} -> <Mean elements>, <Covariance>

    fit_mean_elements!(
        orbp::OrbitPropagator,
        vsv::AbstractVector{OrbitStateVector{Tepoch, T}}[, sink::Type];
        kwargs...
    ) where {Tepoch <: Number, T <: Number} -> <Mean elements>, <Covariance>

Fit a set of mean elements for the propagator `orbp` using the osculating state vector
represented in an inertial reference frame. The state vector can be represented using a set
of position vectors `vr_i` [m] and a set of velocity vectors `vv_i` [m / s] obtained at the
instants in the array `vjd` [Julian Day], or an array of `OrbitStateVector` `vsv` [SI],
containing the same information. The keywords `kwargs` depend on the propagator type.

The optional parameter `sink` selects the representation of the mean elements for the
propagators that support more than one. SGP4 accepts `TLE` and `OrbitMeanElementsMessage`,
defaulting to the latter.

This function also initializes `orbp` with the fitted mean elements.

# Returns

- `<Mean elements>`: Set of mean elements used to initialize `orbp`. The concrete type
    depends on the propagator, and it is `KeplerianElements{MeanAnomaly}` for every
    propagator except SGP4, which returns an `OrbitMeanElementsMessage` or a `TLE` selected
    by `sink`.
- `<Covariance>`: Final covariance matrix of the least-square algorithm. It is a
    `SMatrix{6, 6}` for every propagator except SGP4, which returns a `SMatrix{7, 7}`
    because it also fits the B* parameter.
"""
function fit_mean_elements! end

function fit_mean_elements!(
    orbp::OrbitPropagator, vsv::AbstractVector{OrbitStateVector{Tepoch, T}}; kwargs...
) where {Tepoch <: Number, T <: Number}
    vjd  = map(x -> x.epoch, vsv)
    vr_i = map(x -> x.r, vsv)
    vv_i = map(x -> x.v, vsv)

    return fit_mean_elements!(orbp, vjd, vr_i, vv_i; kwargs...)
end

function fit_mean_elements!(
    orbp::OrbitPropagator,
    vsv::AbstractVector{OrbitStateVector{Tepoch, T}},
    sink::Type;
    kwargs...,
) where {Tepoch <: Number, T <: Number}
    vjd  = map(x -> x.epoch, vsv)
    vr_i = map(x -> x.r, vsv)
    vv_i = map(x -> x.v, vsv)

    return fit_mean_elements!(orbp, vjd, vr_i, vv_i, sink; kwargs...)
end

"""
    init(::Val{:propagator}, args...; kwargs...) -> OrbitPropagator

Create and initialize the orbit `propagator`. The arguments `args` and keywords `kwargs`
depend on the propagator type.
"""
function init end

"""
    init!(orbp::OrbitPropagator, args...; kwargs...) -> Nothing

Initialize the orbit propagator `orbp`. The arguments `args` and keywords `kwargs` depend
on the propagator type.
"""
function init! end

"""
    epoch(orbp::OrbitPropagator{Tepoch, T}) where {Tepoch <: Number, T <: Number} -> Tepoch

Return the initial elements' epoch of the propagator `orbp` [Julian Day, UTC].
"""
function epoch end

"""
    last_instant(
        orbp::OrbitPropagator{Tepoch, T}
    ) where {Tepoch <: Number, T <: Number} -> T

Return the last propagation instant [s] measured from the epoch.
"""
function last_instant end

"""
    is_initialized(orbp::OrbitPropagator) -> Bool

Return whether the orbit propagator `orbp` has been initialized and can be propagated. This
is an optional function in the API, used to print a propagator created without initial
elements without accessing its undefined fields. It returns `true` if the propagator does
not overload it.
"""
is_initialized(orbp::OrbitPropagator) = true

"""
    mean_elements(orbp::OrbitPropagator) -> Union{Nothing, KeplerianElements{MeanAnomaly}}

Return the mean elements using the structure `KeplerianElements{MeanAnomaly}` of the latest
propagation performed by `orbp`. This is an optional function in the API. It will return
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
    propagator_data(orbp::OrbitPropagator) -> Any

Return the structure that stores the data of the propagation theory wrapped by the orbit
propagator `orbp`. This is an optional function in the API, used to print the rich
representation of that structure after the header of `orbp`. It returns `nothing` if the
propagator does not overload it, in which case only the epoch and the last propagation
instant are printed.
"""
propagator_data(orbp::OrbitPropagator) = nothing

"""
    propagate(
        [sink::Type, ]::Val{:propagator},
        t::Union{Number, Dates.Period, Dates.CompoundPeriod},
        args...;
        kwargs...
    ) -> SVector{3, T}, SVector{3, T}, OrbitPropagator{Tepoch, T}

    propagate(
        ::Type{OrbitStateVector},
        ::Val{:propagator},
        t::Union{Number, Dates.Period, Dates.CompoundPeriod},
        args...;
        kwargs...
    ) -> OrbitStateVector{Tepoch, T}, OrbitPropagator{Tepoch, T}

Initialize the orbit `propagator` and propagate the orbit by `t` from the initial orbit
epoch, where `t` is either a number [s] or a period defined using **Dates.jl**. The
initialization arguments `args...` and `kwargs...` are the same as in the initialization
function [`Propagators.init`](@ref). The output type depends on the parameter `sink`. If it
is omitted, it defaults to `Tuple` and the output is a tuple with the position and velocity
vectors.

!!! note

    `T` is the propagator number type. For more information, see [`Propagators.init`](@ref).

# Returns

If `sink` is `Tuple`:

- `SVector{3, T}`: Position vector [m] represented in the inertial frame at propagation
    instant.
- `SVector{3, T}`: Velocity vector [m / s] represented in the inertial frame at propagation
    instant.
- [`OrbitPropagator{Tepoch, T}`](@ref): Structure with the initialized propagator.

If `sink` is `OrbitStateVector`:

- `OrbitStateVector{Tepoch, T}`: Structure with the orbit state vector [SI] at the
    propagation instant.
- [`OrbitPropagator{Tepoch, T}`](@ref): Structure with the initialized propagator.
"""
function propagate(prop::Val, t::PropagationInstant, args...; kwargs...)
    return propagate(Tuple, prop, t, args...; kwargs...)
end

function propagate(
    sink::PropagationSink, prop::Val, t::PropagationInstant, args...; kwargs...
)
    orbp = init(prop, args...; kwargs...)
    return _append_propagator(propagate!(orbp, t, sink), orbp)
end

"""
    propagate(
        [sink::Type, ]::Val{:propagator},
        vt::AbstractVector,
        args...;
        kwargs...
    ) -> Vector{SVector{3, T}}, Vector{SVector{3, T}}, OrbitPropagator{Tepoch, T}

    propagate(
        ::Type{OrbitStateVector},
        ::Val{:propagator},
        vt::AbstractVector,
        args...;
        kwargs...
    ) -> Vector{OrbitStateVector{Tepoch, T}}, OrbitPropagator{Tepoch, T}

Initialize the orbit `propagator` and propagate the orbit for every instant defined in `vt`
from the initial orbit epoch, where the elements of `vt` are either numbers [s] or periods
defined using **Dates.jl**. The initialization arguments `args...` and `kwargs...` are the
same as in the initialization function [`Propagators.init`](@ref). The output type depends
on the parameter `sink`. If it is omitted, it defaults to `Tuple`, and the output is a tuple
with the arrays containing the position and velocity vectors.

!!! note

    `T` is the propagator number type. For more information, see [`Propagators.init`](@ref).

# Keywords

- `ntasks::Integer`: Number of parallel tasks to propagate the orbit. If it is set to a
    number equal or lower than 1, the function will propagate the orbit sequentially. The
    number of tasks is also limited to the number of propagation instants that are actually
    partitioned among them.
    (**Default**: `Threads.nthreads()`)

# Returns

If `sink` is `Tuple`:

- `Vector{SVector{3, T}}`: Array with the position vectors [m] in the inertial frame at each
    propagation instant defined in `vt`.
- `Vector{SVector{3, T}}`: Array with the velocity vectors [m / s] in the inertial frame at
    each propagation instant defined in `vt`.
- [`OrbitPropagator{Tepoch, T}`](@ref): Structure with the initialized propagator.

If `sink` is `OrbitStateVector`:

- `Vector{OrbitStateVector{Tepoch, T}}`: Array with the orbit state vectors [SI] at each
    propagation instant defined in `vt`.
- [`OrbitPropagator{Tepoch, T}`](@ref): Structure with the initialized propagator.
"""
function propagate(prop::Val, vt::AbstractVector, args...; kwargs...)
    return propagate(Tuple, prop, vt, args...; kwargs...)
end

function propagate(
    sink::PropagationSink,
    prop::Val,
    vt::AbstractVector,
    args...;
    ntasks::Integer = Threads.nthreads(),
    kwargs...,
)
    orbp = init(prop, args...; kwargs...)
    return _append_propagator(propagate!(orbp, vt, sink; ntasks = ntasks), orbp)
end

"""
    propagate!(
        orbp::OrbitPropagator{Tepoch, T},
        t::Union{Number, Dates.Period, Dates.CompoundPeriod}[, sink::Type]
    ) where {Tepoch <: Number, T <: Number} -> SVector{3, T}, SVector{3, T}

    propagate!(
        orbp::OrbitPropagator{Tepoch, T},
        t::Union{Number, Dates.Period, Dates.CompoundPeriod},
        ::Type{OrbitStateVector}
    ) where {Tepoch <: Number, T <: Number} -> OrbitStateVector{Tepoch, T}

Propagate the orbit using `orbp` by `t` from the initial orbit epoch, where `t` is either a
number [s] or a period defined using **Dates.jl**. The output type depends on the parameter
`sink`. If it is omitted, it defaults to `Tuple` and the output is a tuple with the position
and velocity vectors.

# Returns

If `sink` is `Tuple`:

- `SVector{3, T}`: Position vector [m] represented in the inertial frame at propagation
    instant.
- `SVector{3, T}`: Velocity vector [m / s] represented in the inertial frame at propagation
    instant.

If `sink` is `OrbitStateVector`:

- `OrbitStateVector{Tepoch, T}`: Structure with the orbit state vector [SI] at the
    propagation instant.
"""
function propagate! end

function propagate!(orbp::OrbitPropagator, p::Union{Dates.Period, Dates.CompoundPeriod})
    return propagate!(orbp, _to_seconds(p))
end

function propagate!(orbp::OrbitPropagator, t::PropagationInstant, ::Type{Tuple})
    return propagate!(orbp, t)
end

function propagate!(orbp::OrbitPropagator, t::PropagationInstant, ::Type{OrbitStateVector})
    Δt = _to_seconds(t)
    r_i, v_i = propagate!(orbp, Δt)
    return OrbitStateVector(epoch(orbp) + Δt / 86400, r_i, v_i)
end

"""
    propagate!(
        orbp::OrbitPropagator{Tepoch, T},
        vt::AbstractVector[, sink::Type];
        kwargs...
    ) where {Tepoch <: Number, T <: Number} -> Vector{SVector{3, T}}, Vector{SVector{3, T}}

    propagate!(
        orbp::OrbitPropagator{Tepoch, T},
        vt::AbstractVector,
        ::Type{OrbitStateVector};
        kwargs...
    ) where {Tepoch <: Number, T <: Number} -> Vector{OrbitStateVector{Tepoch, T}}

Propagate the orbit using `orbp` for every instant defined in `vt` from the initial orbit
epoch, where the elements of `vt` are either numbers [s] or periods defined using
**Dates.jl**. The output type depends on the parameter `sink`. If it is omitted, it defaults
to `Tuple`, and the output is a tuple with the arrays containing the position and velocity
vectors.

The instants are partitioned among parallel tasks, each one propagating a copy of `orbp`.
After the call, `orbp` is left at the last instant in `vt`.

# Keywords

- `ntasks::Integer`: Number of parallel tasks to propagate the orbit. If it is set to a
    number equal or lower than 1, the function will propagate the orbit sequentially. The
    number of tasks is also limited to the number of propagation instants that are actually
    partitioned among them.
    (**Default**: `Threads.nthreads()`)

# Returns

If `sink` is `Tuple`:

- `Vector{SVector{3, T}}`: Array with the position vectors [m] in the inertial frame at each
    propagation instant defined in `vt`.
- `Vector{SVector{3, T}}`: Array with the velocity vectors [m / s] in the inertial frame at
    each propagation instant defined in `vt`.

If `sink` is `OrbitStateVector`:

- `Vector{OrbitStateVector{Tepoch, T}}`: Array with the orbit state vectors [SI] at each
    propagation instant defined in `vt`.
"""
function propagate!(
    orbp::OrbitPropagator, vt::AbstractVector; ntasks::Integer = Threads.nthreads()
)
    return _propagate_vector!(_to_seconds, orbp, vt, ntasks)
end

function propagate!(orbp::OrbitPropagator, vt::AbstractVector, ::Type{Tuple}; kwargs...)
    return propagate!(orbp, vt; kwargs...)
end

function propagate!(
    orbp::OrbitPropagator, vt::AbstractVector, ::Type{OrbitStateVector}; kwargs...
)
    jd₀ = epoch(orbp)
    vr_i, vv_i = propagate!(orbp, vt; kwargs...)
    return _to_state_vectors(t -> jd₀ + _to_seconds(t) / 86400, vt, vr_i, vv_i)
end

"""
    propagate_to_epoch(
        [sink::Type, ]::Val{:propagator},
        epoch::Union{Number, DateTime},
        args...;
        kwargs...
    ) -> SVector{3, T}, SVector{3, T}, OrbitPropagator{Tepoch, T}

    propagate_to_epoch(
        ::Type{OrbitStateVector},
        ::Val{:propagator},
        epoch::Union{Number, DateTime},
        args...;
        kwargs...
    ) -> OrbitStateVector{Tepoch, T}, OrbitPropagator{Tepoch, T}

Initialize the orbit `propagator` and propagate the orbit until `epoch`, which is either a
Julian Day [UTC] or a `DateTime` [UTC]. The initialization arguments `args...` and
`kwargs...` are the same as in the initialization function [`Propagators.init`](@ref). The
output type depends on the parameter `sink`. If it is omitted, it defaults to `Tuple` and
the output is a tuple with the position and velocity vectors.

!!! note

    `T` is the propagator number type. For more information, see [`Propagators.init`](@ref).

# Returns

If `sink` is `Tuple`:

- `SVector{3, T}`: Position vector [m] represented in the inertial frame at propagation
    instant.
- `SVector{3, T}`: Velocity vector [m / s] represented in the inertial frame at propagation
    instant.
- [`OrbitPropagator{Tepoch, T}`](@ref): Structure with the initialized propagator.

If `sink` is `OrbitStateVector`:

- `OrbitStateVector{Tepoch, T}`: Structure with the orbit state vector [SI] at the
    propagation instant.
- [`OrbitPropagator{Tepoch, T}`](@ref): Structure with the initialized propagator.
"""
function propagate_to_epoch(prop::Val, epoch::PropagationEpoch, args...; kwargs...)
    return propagate_to_epoch(Tuple, prop, epoch, args...; kwargs...)
end

function propagate_to_epoch(
    sink::PropagationSink, prop::Val, epoch::PropagationEpoch, args...; kwargs...
)
    orbp = init(prop, args...; kwargs...)
    return _append_propagator(propagate_to_epoch!(orbp, epoch, sink), orbp)
end

"""
    propagate_to_epoch(
        [sink::Type, ]::Val{:propagator},
        vepoch::AbstractVector,
        args...;
        kwargs...
    ) -> Vector{SVector{3, T}}, Vector{SVector{3, T}}, OrbitPropagator{Tepoch, T}

    propagate_to_epoch(
        ::Type{OrbitStateVector},
        ::Val{:propagator},
        vepoch::AbstractVector,
        args...;
        kwargs...
    ) -> Vector{OrbitStateVector{Tepoch, T}}, OrbitPropagator{Tepoch, T}

Initialize the orbit `propagator` and propagate the orbit for every epoch defined in
`vepoch`, whose elements are either Julian Days [UTC] or `DateTime` objects [UTC]. The
initialization arguments `args...` and `kwargs...` are the same as in the initialization
function [`Propagators.init`](@ref). The output type depends on the parameter `sink`. If it
is omitted, it defaults to `Tuple`, and the output is a tuple with the arrays containing the
position and velocity vectors.

!!! note

    `T` is the propagator number type. For more information, see [`Propagators.init`](@ref).

# Keywords

- `ntasks::Integer`: Number of parallel tasks to propagate the orbit. If it is set to a
    number equal or lower than 1, the function will propagate the orbit sequentially. The
    number of tasks is also limited to the number of propagation instants that are actually
    partitioned among them.
    (**Default**: `Threads.nthreads()`)

# Returns

If `sink` is `Tuple`:

- `Vector{SVector{3, T}}`: Array with the position vectors [m] in the inertial frame at each
    epoch defined in `vepoch`.
- `Vector{SVector{3, T}}`: Array with the velocity vectors [m / s] in the inertial frame at
    each epoch defined in `vepoch`.
- [`OrbitPropagator{Tepoch, T}`](@ref): Structure with the initialized propagator.

If `sink` is `OrbitStateVector`:

- `Vector{OrbitStateVector{Tepoch, T}}`: Array with the orbit state vectors [SI] at each
    epoch defined in `vepoch`.
- [`OrbitPropagator{Tepoch, T}`](@ref): Structure with the initialized propagator.
"""
function propagate_to_epoch(prop::Val, vepoch::AbstractVector, args...; kwargs...)
    return propagate_to_epoch(Tuple, prop, vepoch, args...; kwargs...)
end

function propagate_to_epoch(
    sink::PropagationSink,
    prop::Val,
    vepoch::AbstractVector,
    args...;
    ntasks::Integer = Threads.nthreads(),
    kwargs...,
)
    orbp = init(prop, args...; kwargs...)
    return _append_propagator(
        propagate_to_epoch!(orbp, vepoch, sink; ntasks = ntasks), orbp
    )
end

"""
    propagate_to_epoch!(
        orbp::OrbitPropagator{Tepoch, T},
        epoch::Union{Number, DateTime}[, sink::Type]
    ) where {Tepoch <: Number, T <: Number} -> SVector{3, T}, SVector{3, T}

    propagate_to_epoch!(
        orbp::OrbitPropagator{Tepoch, T},
        epoch::Union{Number, DateTime},
        ::Type{OrbitStateVector}
    ) where {Tepoch <: Number, T <: Number} -> OrbitStateVector{Tepoch, T}

Propagate the orbit using `orbp` until `epoch`, which is either a Julian Day [UTC] or a
`DateTime` [UTC]. The output type depends on the parameter `sink`. If it is omitted, it
defaults to `Tuple` and the output is a tuple with the position and velocity vectors.

# Returns

If `sink` is `Tuple`:

- `SVector{3, T}`: Position vector [m] represented in the inertial frame at propagation
    instant.
- `SVector{3, T}`: Velocity vector [m / s] represented in the inertial frame at propagation
    instant.

If `sink` is `OrbitStateVector`:

- `OrbitStateVector{Tepoch, T}`: Structure with the orbit state vector [SI] at the
    propagation instant.
"""
function propagate_to_epoch!(
    orbp::OrbitPropagator, epoch::PropagationEpoch, sink::PropagationSink = Tuple
)
    return propagate!(orbp, _epoch_to_seconds(orbp, epoch), sink)
end

"""
    propagate_to_epoch!(
        orbp::OrbitPropagator{Tepoch, T},
        vepoch::AbstractVector[, sink::Type];
        kwargs...
    ) where {Tepoch <: Number, T <: Number} -> Vector{SVector{3, T}}, Vector{SVector{3, T}}

    propagate_to_epoch!(
        orbp::OrbitPropagator{Tepoch, T},
        vepoch::AbstractVector,
        ::Type{OrbitStateVector};
        kwargs...
    ) where {Tepoch <: Number, T <: Number} -> Vector{OrbitStateVector{Tepoch, T}}

Propagate the orbit using `orbp` for every epoch defined in `vepoch`, whose elements are
either Julian Days [UTC] or `DateTime` objects [UTC]. The output type depends on the
parameter `sink`. If it is omitted, it defaults to `Tuple`, and the output is a tuple with
the arrays containing the position and velocity vectors.

The epochs are partitioned among parallel tasks, each one propagating a copy of `orbp`.
After the call, `orbp` is left at the last epoch in `vepoch`.

# Keywords

- `ntasks::Integer`: Number of parallel tasks to propagate the orbit. If it is set to a
    number equal or lower than 1, the function will propagate the orbit sequentially. The
    number of tasks is also limited to the number of propagation instants that are actually
    partitioned among them.
    (**Default**: `Threads.nthreads()`)

# Returns

If `sink` is `Tuple`:

- `Vector{SVector{3, T}}`: Array with the position vectors [m] in the inertial frame at each
    epoch defined in `vepoch`.
- `Vector{SVector{3, T}}`: Array with the velocity vectors [m / s] in the inertial frame at
    each epoch defined in `vepoch`.

If `sink` is `OrbitStateVector`:

- `Vector{OrbitStateVector{Tepoch, T}}`: Array with the orbit state vectors [SI] at each
    epoch defined in `vepoch`.
"""
function propagate_to_epoch!(
    orbp::OrbitPropagator, vepoch::AbstractVector; ntasks::Integer = Threads.nthreads()
)
    jd₀ = epoch(orbp)
    return _propagate_vector!(t -> 86400 * (_to_julian_day(t) - jd₀), orbp, vepoch, ntasks)
end

function propagate_to_epoch!(
    orbp::OrbitPropagator, vepoch::AbstractVector, ::Type{Tuple}; kwargs...
)
    return propagate_to_epoch!(orbp, vepoch; kwargs...)
end

function propagate_to_epoch!(
    orbp::OrbitPropagator, vepoch::AbstractVector, ::Type{OrbitStateVector}; kwargs...
)
    vr_i, vv_i = propagate_to_epoch!(orbp, vepoch; kwargs...)
    return _to_state_vectors(_to_julian_day, vepoch, vr_i, vv_i)
end

"""
    step!(
        orbp::OrbitPropagator{Tepoch, T},
        Δt::Union{Number, Dates.Period, Dates.CompoundPeriod}[, sink::Type]
    ) where {Tepoch <: Number, T <: Number} -> SVector{3, T}, SVector{3, T}

    step!(
        orbp::OrbitPropagator{Tepoch, T},
        Δt::Union{Number, Dates.Period, Dates.CompoundPeriod},
        ::Type{OrbitStateVector}
    ) where {Tepoch <: Number, T <: Number} -> OrbitStateVector{Tepoch, T}

Propagate the orbit using `orbp` by `Δt` from the last propagation instant, where `Δt` is
either a number [s] or a period defined using **Dates.jl**. The output type depends on the
parameter `sink`. If it is omitted, it defaults to `Tuple` and the output is a tuple with
the position and velocity vectors.

# Returns

If `sink` is `Tuple`:

- `SVector{3, T}`: Position vector [m] represented in the inertial frame at propagation
    instant.
- `SVector{3, T}`: Velocity vector [m / s] represented in the inertial frame at propagation
    instant.

If `sink` is `OrbitStateVector`:

- `OrbitStateVector{Tepoch, T}`: Structure with the orbit state vector [SI] at the
    propagation instant.
"""
function step!(orbp::OrbitPropagator, Δt::PropagationInstant, sink::PropagationSink = Tuple)
    return propagate!(orbp, last_instant(orbp) + _to_seconds(Δt), sink)
end

############################################################################################
#                                           Copy                                           #
############################################################################################

Base.copy(orbp::OrbitPropagator) = deepcopy(orbp)

############################################################################################
#                                    Iterator Interface                                    #
############################################################################################

iterate(orbp::OrbitPropagator) = (orbp, nothing)
iterate(orbp::OrbitPropagator, ::Nothing) = nothing
length(orbp::OrbitPropagator) = 1
eltype(orbp::T) where {T <: OrbitPropagator} = T

############################################################################################
#                                           Show                                           #
############################################################################################

function show(io::IO, orbp::OrbitPropagator)
    header = _header(orbp)

    if !is_initialized(orbp)
        print(io, header, " (not initialized)")
        return nothing
    end

    SatelliteToolboxBase.print_compact(io, header, epoch(orbp))

    return nothing
end

function show(io::IO, ::MIME"text/plain", orbp::OrbitPropagator)
    header = _header(orbp)
    data   = propagator_data(orbp)

    # If the propagator provides its data structure, we print its body under the header of
    # the wrapper.
    if !isnothing(data)
        SatelliteToolboxBase.print_tree(io, header, data)
        return nothing
    end

    # Otherwise, we can only print the information obtained through the API.
    if is_initialized(orbp)
        epoch_str = SatelliteToolboxBase.epoch_string(epoch(orbp))
        Δt_str    = SatelliteToolboxBase.format_value(last_instant(orbp))

        fields = SatelliteToolboxBase.PrintedField[
            ("Epoch",            epoch_str, ""),
            ("Last Propagation", Δt_str,    "s"),
        ]
    else
        fields = SatelliteToolboxBase.PrintedField[("Status", "not initialized", "")]
    end

    sections = SatelliteToolboxBase.PrintedSection[]
    SatelliteToolboxBase.print_tree(io, header, fields, sections)

    return nothing
end

############################################################################################
#                                    Private Functions                                     #
############################################################################################

"""
    _header(orbp::OrbitPropagator) -> String

Return the header of the printed representations of `orbp`: its type with the parameters,
followed by the propagator name in parentheses if the propagator defines one.
"""
function _header(orbp::OrbitPropagator)
    type_name = _type_name(orbp)
    prop_name = name(orbp)
    prop_name == type_name && return type_name
    return string(type_name, " (", prop_name, ")")
end

"""
    _type_name(x) -> String

Return the name of the type of `x` with its parameters, as used in the headers of the
printed representations.
"""
function _type_name(x)
    T = typeof(x)
    return string(nameof(T), "{", join(T.parameters, ", "), "}")
end

"""
    _append_propagator(result, orbp::OrbitPropagator) -> Tuple

Append the propagator `orbp` to the `result` of a propagation so that the functions that
initialize and propagate the orbit at once return the initialized propagator as their last
element. A `result` that is a tuple is splatted, whereas any other object becomes the first
element of the returned tuple.
"""
_append_propagator(result::Tuple, orbp::OrbitPropagator) = (result..., orbp)
_append_propagator(result, orbp::OrbitPropagator) = (result, orbp)

"""
    _epoch_to_seconds(orbp::OrbitPropagator, epoch::Union{Number, DateTime}) -> Number

Convert `epoch`, which is either a Julian Day [UTC] or a `DateTime` [UTC], to the elapsed
time [s] from the initial elements' epoch of `orbp`.
"""
function _epoch_to_seconds(orbp::OrbitPropagator, epoch_::PropagationEpoch)
    return 86400 * (_to_julian_day(epoch_) - epoch(orbp))
end

"""
    _propagate_vector!(
        to_Δt,
        orbp::OrbitPropagator{Tepoch, T},
        vt::AbstractVector,
        ntasks::Integer
    ) where {Tepoch <: Number, T <: Number} -> Vector{SVector{3, T}}, Vector{SVector{3, T}}

Propagate the orbit using `orbp` for every element of `vt`, which the function `to_Δt`
converts to an elapsed time [s] from the initial elements' epoch. The first and the last
elements are propagated with `orbp` itself, so that it is left at the last instant, whereas
the remaining ones are partitioned among `ntasks` parallel tasks, each one propagating a
copy of `orbp`. If `ntasks` is lower than 1, the propagation is sequential.

# Returns

- `Vector{SVector{3, T}}`: Array with the position vectors [m] in the inertial frame at each
    instant defined in `vt`.
- `Vector{SVector{3, T}}`: Array with the velocity vectors [m / s] in the inertial frame at
    each instant defined in `vt`.
"""
function _propagate_vector!(
    to_Δt, orbp::OrbitPropagator, vt::AbstractVector, ntasks::Integer
)
    # We need to perform the first propagation to obtain the return type of the propagator.
    r₀, v₀ = propagate!(orbp, to_Δt(first(vt)))

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

    # The first and the last instants are propagated separately. Hence, only `len_vt - 2`
    # instants are partitioned among the tasks. We must not create more tasks than that,
    # otherwise the surplus tasks would be assigned the same partition and would write
    # concurrently to the same output elements. We also must have at least one task,
    # otherwise no instant would be propagated at all.
    num_tasks = max(min(Int(ntasks), len_vt - 2), 1)

    # If we have only two instants in the time vector, we will not spawn any threads,
    # because the first and the last instants are propagated separately.
    if len_vt > 2
        # The propagation usually modifies the structure. Hence, we need one propagator per
        # task. We copy them here because `orbp` is mutated by the task that uses it, and
        # copying inside the loop would race with that mutation.
        corbps = [c == 1 ? orbp : copy(orbp) for c in 1:num_tasks]

        @maybe_threads num_tasks for c in 1:num_tasks
            # We already propagated for the first instant, and we must ensure we propagate
            # the last instant at the end of the function.
            i₀, i₁ = @views get_partition(c, inds[(1 + begin):(end - 1)], num_tasks)

            corbp = corbps[c]

            # The indices come from a partition of `eachindex(vt)`, and the output vectors
            # have the same length as `vt`. Hence, we can skip the bounds checking.
            @inbounds for i in i₀:i₁
                vr[i - Δi], vv[i - Δi] = propagate!(corbp, to_Δt(vt[i]))
            end
        end
    end

    # We must ensure that the last propagation instant is the one obtained at the end to
    # keep the internal data of the propagation consistent.
    vr[end], vv[end] = propagate!(orbp, to_Δt(last(vt)))

    return vr, vv
end

"""
    _to_julian_day(epoch::Number) -> Number
    _to_julian_day(epoch::DateTime) -> Float64

Convert `epoch`, which is either a Julian Day [UTC] or a `DateTime` [UTC], to a Julian Day
[UTC].
"""
_to_julian_day(epoch::Number) = epoch
_to_julian_day(epoch::DateTime) = datetime2julian(epoch)

"""
    _to_seconds(t::Number) -> Number
    _to_seconds(t::Union{Dates.Period, Dates.CompoundPeriod}) -> Float64

Convert `t`, which is either a number [s] or a period defined using **Dates.jl**, to
seconds.

`Dates.toms` is type unstable for compound periods (see
https://github.com/JuliaLang/julia/pull/54995), so the conversion is implemented here.
"""
_to_seconds(t::Number) = t
_to_seconds(p::Dates.Period) = Dates.toms(p) / 1000

function _to_seconds(p::Dates.CompoundPeriod)
    return (isempty(p.periods) ? 0.0 : Float64(sum(Dates.toms, p.periods))) / 1000
end

"""
    _to_state_vectors(
        to_epoch,
        vt::AbstractVector,
        vr_i::AbstractVector,
        vv_i::AbstractVector
    ) -> Vector{OrbitStateVector}

Pack the position vectors `vr_i` [m] and the velocity vectors `vv_i` [m / s] into orbit
state vectors whose epochs [Julian Day] are obtained by applying the function `to_epoch` to
the elements of `vt`. The vectors `vr_i` and `vv_i` are 1-based, whereas `vt` can have any
axes.
"""
function _to_state_vectors(
    to_epoch, vt::AbstractVector, vr_i::AbstractVector, vv_i::AbstractVector
)
    Δi = firstindex(vt) - 1

    return map(eachindex(vr_i)) do k
        return OrbitStateVector(to_epoch(vt[k + Δi]), vr_i[k], vv_i[k])
    end
end

end # module Propagators
