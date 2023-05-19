API for the Orbit Propagators
=============================

This document describes the design of orbit propagators.

## Initialization

Each propagator must have an initialization function that receives the input elements and
returns a propagator structure. This entity contains all the initialized variables related
to the propagator.

The type of the input elements varies according to the propagator.

No restriction is imposed on the output structure, which can also contain variables to
reduce the computational burden.

This function shall be named as `<propagator identifier>_init`.

```julia
j2d = j2_init(orb; j2c = j2c_egm08)
```

The propagator can also implement the in-place initialization, where a structure will be
provided and must have its values re-initialized:

```julia
j2_init!(j2d, orb)
```

This version can be used to reduce the number of allocations.

## Propagation

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

## Simultaneous Initialization and Propagation (Optional)

The propagator can defined the following function to simultaneously initialize the
propagator and propagate the orbit: `<propagator identifier>`. Its first argument must be
the elapsed time from the epoch associated with the input elements, and the others must be
exaclty the same as the initialization function.

All the considerations related to the propagation function also apply here.

This function must return the same result as the propagator function and also the
initialized propagator structure.

```julia
r_i, v_i, j2d = j2(10, orb; j2c = j2c_egm08)
```

## API

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

### Initialization

The initialization is performed by the function:

```julia
Propagators.init(T, args...; kwargs...)
```

where `T = Val(<Orbit propagator symbol>)`, and it must return and object of type
`OrbitPropagator<Propagator name>`. The arguments and keywords depends on the propagator and
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

> **Note**
> The API provides a default implementation for `Propagators.step!`. Hence, strictly
> speaking, only the implementation of `Propagators.propagate!` is required for the API.
> However, there are some cases in which it is more accurate to implement
> `Propagators.step!` and use this algorithm to build the function `Propagators.propagate!`.
> In those cases, the user must overload both functions.

We also have the function `Propagators.propagate_to_epoch!`, but the default implementation
should work for all propagators.

### In-place initialization (Optional)

If the propagator supports in-place initialization, it must overload the following function:

```julia
Propagators.init!(orbp::OrbitPropagator<Propagator name>, args...; kwargs...)
```

### Mean Elements (Optional)

The function

```julia
Propagators.mean_elements(orbp)
```

should return the mean elements using the structure `KeplerianElements` related to the
latest propagation. Notice that this is an **optional** feature. If the propagator does not
implement it, it will return `nothing`.

### Name (Optional)

The propagator can overload the function:

```julia
Propagators.name(orbp)
```

to return its name. The system uses this information to display the object using the
function `show`. If the function is not provided, the structure name is used by default.

### Simultaneous Initialization and Propagation (Optional)

If the propagator supports, it can overload the functions:

```julia
Propagators.propagate(t, args...; kwargs...)
```

that simultaneously initialize and propagate the orbit to the instant `t` [s] after the
input elements. `args...` and `kwargs...` must be the same as in the initialization function
`Propagators.init`.
