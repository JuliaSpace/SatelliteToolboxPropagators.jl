# API

**SatelliteToolboxPropagators.jl** has a propagator API to improve the usability. This API
has the following requirements, which must be met by all propagators that uses it.

Every propagator must have a structure derived from `OrbitPropagator` with the following
requirement:

```julia
struct OrbitPropagator<Propagator name>{Tepoch, T} <: OrbitPropagator{Tepoch, T}
    <Any field required by the propagator>
end
```

where `Tepoch` is the type used to represent the epoch of the input elements, whereas `T` is
the type used for the internal variables.

## Initialization

The initialization is performed by the function:

```julia
Propagators.init(T, args...; kwargs...)
```

where `T = Val(<Orbit propagator symbol>)`, and it must return and object of type
`OrbitPropagator<Propagator name>`. The arguments and keywords depends on the propagator and
must be documented in the docstring. The propagator must record the epoch during the
initialization, which must be kept constant during the entire object existence. It also
needs to record the instant of the last propagation.

## Epoch

Each propagator must return the initial element epoch in Julian Day by the function:

```julia
Propagators.epoch(orbp)
```

Notice that this value must never change during the existence of the object `orbp`.

## Last Instant

Each propagator must return the last propagation instant [s] measured from the epoch. This
action must be performed by the function:

```julia
Propagators.last_instant(orbp)
```

## Propagation

The following functions must be overloaded by each propagator.

```julia
Propagators.propagate!(orbp, t)
```

Propagate the orbit of propagator `orbp` by `t` [s] from the epoch. This function must
return the propagated position and velocity represented in the same reference frame used in
the initialization. **The output vector must use the SI**.

```julia
Propagators.step!(orbp, dt)
```

Propagate the orbit of propagator `orbp` by `dt` [s] from the instant of the last
propagation. This function must return the propagated position and velocity represented in
the same reference frame used in the initialization. **The output vector must use the SI**.

> **Note**
> The API provides a default implementation for `Propagators.step!`. Hence, strictly
> speaking, only the implementation of `Propagators.propagate!` is required for the API.
> However, there are some cases in which it is more accurate to implement
> `Propagators.step!` and use this algorithm to build the function `Propagators.propagate!`.
> In those cases, the user must overload both functions.

We also have the function `Propagators.propagate_to_epoch!`, but the default implementation
should work for all propagators.

## In-place initialization (Optional)

If the propagator supports in-place initialization, it must overload the following function:

```julia
Propagators.init!(orbp::OrbitPropagator<Propagator name>, args...; kwargs...)
```

## Mean Elements (Optional)

The function

```julia
Propagators.mean_elements(orbp)
```

should return the mean elements using the structure `KeplerianElements` related to the
latest propagation. Notice that this is an **optional** feature. If the propagator does not
implement it, it will return `nothing`.

## Name (Optional)

The propagator can overload the function:

```julia
Propagators.name(orbp)
```

to return its name. The system uses this information to display the object using the
function `show`. If the function is not provided, the structure name is used by default.

## Simultaneous Initialization and Propagation (Optional)

If the propagator supports, it can overload the functions:

```julia
Propagators.propagate(t, args...; kwargs...)
```

that simultaneously initialize and propagate the orbit to the instant `t` [s] after the
input elements. `args...` and `kwargs...` must be the same as in the initialization function
`Propagators.init`.

### Fitting Mean Elements (Optional)

The propagator can implement the following functions to fit a set of osculating state
vectors into mean elements for its theory:

```julia
Propagators.fit_mean_elements(T, vjd::AbstractVector{Tjd}, vr_i::AbstractVector{Tv}, vv_i::AbstractVector{Tv}; kwargs...) where {Tjd<:Number, Tv<:AbstractVector}
Propagators.fit_mean_elements!(orbp::OrbitPropagator<Propagator name>, vjd::AbstractVector{Tjd}, vr_i::AbstractVector{Tv}, vv_i::AbstractVector{Tv}; kwargs...) where {Tjd<:Number, Tv<:AbstractVector}
```

where `T = Val(<Orbit propagator symbol>)`, `vr_i` and `vv_i` are a set of position [m] and
velocity [m / s] vectors obtained at the instants in `vjd` [Julian Day]. Those functions
must return the mean elements used to initialize the propagator in the function
`Propagator.init`.

Each propagator type can define its own set of keyword arguments to configure the fitting
process.

The first signature will allocate a new propagator, whereas the second will use the
allocated one passed as the first argument. In the latter, the propagator needs to be
initialized with the fitted elements.
