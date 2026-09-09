# API

```@meta
CurrentModule = SatelliteToolboxPropagators
```

This document describes the design of the orbit propagators and the requirements a
propagator must meet to be used through the `Propagators` API.

## Low-Level Propagators

Each propagator has a set of functions that work directly with its data structure.

### Initialization

Each propagator must have an initialization function that receives the input elements and
returns a propagator structure. This entity contains all the initialized variables related
to the propagator.

The type of the input elements varies according to the propagator.

No restriction is imposed on the output structure, which can also contain variables to
reduce the computational burden.

This function shall be named as `<propagator identifier>_init`.

```julia
j2d = j2_init(orb; j2c = J2C_EGM2008)
```

The propagator can also implement the in-place initialization, where a structure will be
provided and must have its values re-initialized:

```julia
j2_init!(j2d, orb)
```

This version can be used to reduce the number of allocations.

### Propagation

Each propagator must have a propagation function that receives the propagator structure and
a time `t`. It must return the position and velocity after `t` [T] from the initialized
epoch. The unit [T] of `t` varies according to the propagator.

The propagation must be performed ideally in the same reference frame in which the input
elements were represented. If this is not the case, the documentation must clearly state it.

The output vectors must be in the same reference frame in which the input elements were
represented. Hence, if the propagator transforms the inputs to another frame, then it must
convert them back.

For the sake of simplification, it is advised, although not imposed, that the propagator
structure is updated with the latest orbital elements using any representation.

The propagator should return the vectors preferably in SI. If this is not the case, the
documentation must state clearly.

This function shall be named as `<propagator identifier>!`.

```julia
r_i, v_i = j2!(j2d, 10)
```

### Simultaneous Initialization and Propagation (Optional)

The propagator can define the following function to simultaneously initialize the
propagator and propagate the orbit: `<propagator identifier>`. Its first argument must be
the elapsed time from the epoch associated with the input elements, and the others must be
exactly the same as the initialization function.

All the considerations related to the propagation function also apply here.

This function must return the same result as the propagator function and also the
initialized propagator structure.

```julia
r_i, v_i, j2d = j2(10, orb; j2c = J2C_EGM2008)
```

## Propagators API

**SatelliteToolboxPropagators.jl** has a propagator API to improve the usability. This API
has the following requirements, which must be met by all propagators that use it.

Every propagator must have a structure derived from `OrbitPropagator` with the following
requirement:

```julia
struct OrbitPropagator<Propagator name>{Tepoch <: Number, T <: Number} <: OrbitPropagator{Tepoch, T}
    <Any field required by the propagator>
end
```

where `Tepoch` is the type used to represent the epoch of the input elements, whereas `T` is
the type used for the internal variables.

### Initialization

The initialization is performed by the function:

```julia
Propagators.init(T, args...; kwargs...)
```

where `T = Val(<Orbit propagator symbol>)`, and it must return an object of type
`OrbitPropagator<Propagator name>`. The arguments and keywords depend on the propagator and
must be documented in the docstring. The propagator must record the epoch during the
initialization, which must be kept constant during the entire object existence. It also
needs to record the instant of the last propagation.

### Epoch

Each propagator must return the initial element epoch in Julian Day by the function:

```julia
Propagators.epoch(orbp)
```

Notice that this value must never change during the existence of the object `orbp`.

### Last Instant

Each propagator must return the last propagation instant [s] measured from the epoch. This
action must be performed by the function:

```julia
Propagators.last_instant(orbp)
```

### Propagation

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

!!! note

    The API provides a default implementation for `Propagators.step!`. Hence, strictly
    speaking, only the implementation of `Propagators.propagate!` is required for the API.
    However, there are some cases in which it is more accurate to implement
    `Propagators.step!` and use this algorithm to build the function `Propagators.propagate!`.
    In those cases, the user must overload both functions.

The API also provides the vectorized versions of `Propagators.propagate!`, the functions
`Propagators.propagate_to_epoch!`, the support for `Dates` objects, and the output sinks.
Their default implementations, built on `Propagators.propagate!`, should work for all
propagators.

### In-place Initialization (Optional)

If the propagator supports in-place initialization, it must overload the following function:

```julia
Propagators.init!(orbp::OrbitPropagator<Propagator name>, args...; kwargs...)
```

### Mean Elements (Optional)

The function

```julia
Propagators.mean_elements(orbp)
```

should return the mean elements using the structure `KeplerianElements{MeanAnomaly}` related
to the latest propagation. Notice that this is an **optional** feature. If the propagator
does not implement it, it will return `nothing`.

### Name (Optional)

The propagator can overload the function:

```julia
Propagators.name(orbp)
```

to return its name. The system uses this information to display the object using the
function `show`. If the function is not provided, the structure name is used by default.

### Printing (Optional)

The function `show` prints the propagator using the type with its parameters, the name, and
the epoch. The rich representation, `show(io, MIME("text/plain"), orbp)`, prints the same
header followed by the body of the rich representation of the structure returned by the
function:

```julia
Propagators.propagator_data(orbp)
```

which is printed with `SatelliteToolboxBase.print_tree_body`. Hence, that structure must
overload this function, building its fields and sections with the printing helpers of
**SatelliteToolboxBase.jl**. If `Propagators.propagator_data` is not provided, it returns
`nothing` and only the epoch and the last propagation instant are printed.

The function:

```julia
Propagators.is_initialized(orbp)
```

should return `false` when the propagator was created without initial elements and cannot
be propagated yet, so that the compact representation prints the status instead of
accessing undefined fields. If the function is not provided, it returns `true`.

### Fitting Mean Elements (Optional)

The propagator can implement the following functions to fit a set of osculating state
vectors into mean elements for its theory:

```julia
Propagators.fit_mean_elements(T, vjd::AbstractVector{Tjd}, vr_i::AbstractVector{Tv}, vv_i::AbstractVector{Tv}; kwargs...) where {Tjd <: Number, Tv <: AbstractVector}
Propagators.fit_mean_elements!(orbp::OrbitPropagator<Propagator name>, vjd::AbstractVector{Tjd}, vr_i::AbstractVector{Tv}, vv_i::AbstractVector{Tv}; kwargs...) where {Tjd <: Number, Tv <: AbstractVector}
```

where `T = Val(<Orbit propagator symbol>)`, `vr_i` and `vv_i` are a set of position [m] and
velocity [m / s] vectors obtained at the instants in `vjd` [Julian Day]. Those functions
must return the mean elements used to initialize the propagator in the function
`Propagators.init`, the covariance matrix of the least-square algorithm, and a `NamedTuple`
with its statistics, whose fields are `converged::Bool`, `iterations::Int`,
`position_rmse` [m], `velocity_rmse` [m / s], and `total_rmse`. The propagators that use
Keplerian elements return them as `KeplerianElements{MeanAnomaly}`.

Each propagator type can define its own set of keyword arguments to configure the fitting
process.

The first signature will allocate a new propagator, whereas the second will use the
allocated one passed as the first argument. In the latter, the propagator needs to be
initialized with the fitted elements.
