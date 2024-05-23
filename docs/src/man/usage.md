# Usage

```@meta
CurrentModule = SatelliteToolboxPropagators
```

```@repl usage
using SatelliteToolboxPropagators
```

All the propagators can be accessed using the available API, which allows to initialize and
propagate the orbit.

All the API function are available inside the module `Propagators` that is exported by this
package.

## Initialization

We can initialize an orbit propagator using the function:

```julia
Propagators.init(::Val{:propagator}, args...; kwargs...)
```

It initializes a propagator of type `:propagator` using the arguments `args...` and keywords
`kwargs...` supported by it. Currently, the following algorithms are available:

|                       **Propagator Name** | **Symbol** |
|------------------------------------------:|:-----------|
|            J2 analytical orbit propagator | `:J2`      |
| J2 osculating analytical orbit propagator | `:J2osc`   |
|            J4 analytical orbit propagator | `:J4`      |
| J4 osculating analytical orbit propagator | `:J4osc`   |
|                SGP4/SDP4 orbit propagator | `:SGP4`    |
|      Two body analytical orbit propagator | `:TwoBody` |

See the documentation of each algorithm to verify the supported arguments and keywords.

For example, a J2 analytical orbit propagator can be initialized with a set of mean elements
as follows:

```@repl usage
orb = KeplerianElements(
    date_to_jd(2023, 1, 1, 0, 0, 0),
    7190.982e3,
    0.001111,
    98.405 |> deg2rad,
    100    |> deg2rad,
    90     |> deg2rad,
    19     |> deg2rad
)

orbp = Propagators.init(Val(:J2), orb)
```

We also have the function:

```julia
init!(orbp::OrbitPropagator, args...; kwargs...)
```

that initializes the propagator `orbp` in-place using the arguments `args...` and keywords
`kwargs...`. This function allows to reduce the memory allocation when performing multiple
initialization. However, notice that this function is **optional** in the propagator API.
Hence, some algorithms might not support it.

## Propagate the Orbit

After the initialization, we can propagate the orbit using some functions as follows.

```julia
propagate!(orbp::OrbitPropagator{Tepoch, T}, t::Number) where {Tepoch, T}
```

This function propagates the orbit using `orbp` by `t` [s] from the initial orbit epoch. It
returns two vectors: the position [m] and velocity [m/s]. Both are represented in the same
inertial reference frame used to describe the input data during the propagator
initialization.

For example, using the initialized propagator, we can propagate the orbit by 6000s as
follows:

```@repl usage
Propagators.propagate!(orbp, 6000)
```

---

```julia
propagate_to_epoch!(orbp::OrbitPropagator{Tepoch, T}, jd::Number) where {Tepoch, T}
```

This function propagates the orbit using `orbp` until the epoch `jd` [Julian Day]. It
returns two vectors: the position [m] and velocity [m/s]. Both are represented in the same
inertial reference frame used to describe the input data during the propagator
initialization.

For example, using the initialized propagator, we can propagate the orbit to
2023-01-02T00:00:00 as follows:

```@repl usage
Propagators.propagate_to_epoch!(orbp, date_to_jd(2023, 1, 2, 0, 0, 0))
```

---

```julia
step!(orbp::OrbitPropagator{Tepoch, T}, Δt::Number) where {Tepoch, T}
```

Every time we propagate the orbit, the propagation instant is recorded inside the structure.
Hence, we can use the function `step!` to propagate the orbit by `Δt` [s] from the current
orbit epoch. It returns two vectors: the position [m] and velocity [m/s]. Both are
represented in the same inertial reference frame used to describe the input data during the
propagator initialization.

For example, using the initialized propagator, we can advance the propagation by 60 s as
follows:

```@repl usage
orbp

Propagators.step!(orbp, 60)

orbp
```

## Simultaneous Initialization and Propagation

We can use the following functions to simultaneously initialize and propagate the orbit:

```julia
propagate(::Val{:propagator}, Δt::Number, args...; kwargs...)
propagate_to_epoch(::Val{:propagator}, jd::Number, args...; kwargs...)
```

The symbol `:propagator`, the arguments `args...`, and the keywords `kwargs...` are the same
as described in the section related to the propagator initialization.

The first function propagates the orbit by `Δt` [s] from the initial epoch, whereas the
second propagated the orbit until the instant `jd` [Julian Day].

Both functions returns the position vector [m], the velocity vector [m/s], and the
initialized propagator structure. Notice that both vectors are represented in the same
inertial reference frame used to describe the input data for the initialization.

However, notice that those functions are **optional** in the propagator API. Hence, some
propagators might not support them.

For example, we can propagate by 60 s a set of mean elements using a J2 analytical
propagator as follows:

```@repl usage
Propagators.propagate(Val(:J2), 60, orb)
```

## Helpers

We have the following functions that provide some useful information related to the
propagators.

```julia
epoch(orbp::OrbitPropagator{Tepoch, T}) where {Tepoch<:Number, T<:Number}
```

It returns the initial elements' epoch of the propagator `orbp` [JD].

```@repl usage
Propagators.epoch(orbp)
```

---

```julia
last_instant(orbp::OrbitPropagator{Tepoch, T}) where {Tepoch<:Number, T<:Number}
```

It returns the last propagation instant [s] measured from the input elements' epoch.

```@repl usage
Propagators.last_instant(orbp)
```

---

```julia
name(orbp::OrbitPropagator)
```

It returns the propagator name.

```@repl usage
Propagators.name(orbp)
```

---

```julia
mean_elements(orbp::OrbitPropagator)
```

It returns the mean elements using the structure `KeplerianElements` of the latest
propagation performed by `orbp`. Notice that this is an optional funciton in the
propagators' API. If a propagator does not support it, this function returns `nothing`.

```@repl usage
Propagators.mean_elements(orbp)
```
