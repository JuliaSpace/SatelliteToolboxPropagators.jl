# Usage

```@meta
CurrentModule = SatelliteToolboxPropagators
```

```@repl usage
using SatelliteToolboxPropagators
using Dates
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
propagate!(orbp::OrbitPropagator, Δt::Number[, sink = Tuple])
propagate!(orbp::OrbitPropagator, p::Union{Dates.Period, Dates.CompoundPeriod}[, sink = Tuple])
```

This function propagates the orbit using `orbp` by `t` [s] or by the period defined by `p`
from the initial orbit epoch .

The return type depends on the parameter `sink`. If `sink` is `Tuple` (default), the
function returns two vectors: the position [m] and velocity [m / s]. Both are represented in
the same inertial reference frame used to describe the input data during the propagator
initialization. On the other hand, if `sink` is `OrbitStateVector`, it returns an object of
type `OrbitStateVector` [SI] with the same information.

For example, using the initialized propagator, we can propagate the orbit by 6000s as
follows:

```@repl usage
Propagators.propagate!(orbp, 6000)

Propagators.propagate!(orbp, Dates.Hour(1))

Propagators.propagate!(orbp, Dates.Minute(60))
```

To obtain an `OrbitStateVector` as the output, we can use:

```@repl usage
Propagators.propagate!(orbp, 6000, OrbitStateVector)

Propagators.propagate!(orbp, Dates.Hour(1), OrbitStateVector)

Propagators.propagate!(orbp, Dates.Minute(60), OrbitStateVector)
```

The API also supports propagating the orbit for multiple instants if we pass a vector of
instants as follows:

```julia
propagate!(orbp::OrbitPropagator, vt::AbstractVector[, sink = Tuple]; kwargs...)
propagate!(orbp::OrbitPropagator, vp::AbstractVector{Union{Dates.Period, Dates.CompundPeriod}}[, sink = Tuple]; kwargs...)
```

Those functions propagate the orbit using `orbp` for every instant defined in `vt` [s] or
for every period defined in `vp` from the initial orbit epoch.

The return type depends on the parameter `sink`. If it is `Tuple` (default), the output is a
tuple with the arrays containing the position and velocity vectors. If it is
`OrbitStateVector`, the output is an array of `OrbitStateVector` [SI] objects.

In this case, the user can pass the keyword argument `ntasks` to specify the number of
parallel tasks to propagate the orbit. If it is set to a number equal or lower than 1, the
function will propagate the orbit sequentially. If `ntasks` is omitted, it defaults to the
number of available threads.

For example, using the initialized propagator, we can propagate the orbit every 5 min in an
hour as follows:

```@repl usage
Propagators.propagate!(orbp, 0:300:3600)

Propagators.propagate!(orbp, Dates.Minute(0):Dates.Minute(5):Dates.Minute(60))
```

To obtain an array of `OrbitStateVector` as the output, we can use:

```@repl usage
Propagators.propagate!(orbp, 0:300:3600, OrbitStateVector)

Propagators.propagate!(
    orbp,
    Dates.Minute(0):Dates.Minute(5):Dates.Minute(60),
    OrbitStateVector
)
```

---

```julia
propagate_to_epoch!(orbp::OrbitPropagator, jd::Number[, sink = Tuple])
propagate_to_epoch!(orbp::OrbitPropagator, dt::DateTime[, sink = Tuple])
```

This function propagates the orbit using `orbp` until the Julian Day `jd` [UTC] or until the
epoch defined by the `DateTime` object `dt` [UTC].

The return type depends on the parameter `sink`. If `sink` is `Tuple` (default), the
function returns two vectors: the position [m] and velocity [m / s]. Both are represented in
the same inertial reference frame used to describe the input data during the propagator
initialization. On the other hand, if `sink` is `OrbitStateVector`, it returns an object of
type `OrbitStateVector` [SI] with the same information.

For example, using the initialized propagator, we can propagate the orbit to
2023-01-02T00:00:00 as follows:

```@repl usage
Propagators.propagate_to_epoch!(orbp, date_to_jd(2023, 1, 2, 0, 0, 0))

Propagators.propagate_to_epoch!(orbp, DateTime("2023-01-02"))
```

To obtain an `OrbitStateVector` as the output, we can use:

```@repl usage
Propagators.propagate_to_epoch!(orbp, date_to_jd(2023, 1, 2, 0, 0, 0), OrbitStateVector)

Propagators.propagate_to_epoch!(orbp, DateTime("2023-01-02"), OrbitStateVector)
```

The API also supports propagating the orbit for multiple instants if we pass a vector of
epochs as follows:

```julia
propagate_to_epoch!(orbp::OrbitPropagator, vjd::AbstractVector[, sink = Tuple]; kwargs...)
propagate_to_epoch!(orbp::OrbitPropagator, vdt::AbstractVector{DateTime}[, sink = Tuple]; kwargs...)
```

Those functions propagate the orbit using `orbp` for every epoch defined in the vector of
Julian Days `vjd` [UTC] or in the vector of `DateTime` objects `vdt` [UTC].

The return type depends on the parameter `sink`. If it is `Tuple` (default), the output is
two arrays containing the position [m] and velocity vectors [m / s]. Both are represented in
the same inertial reference frame used to describe the input data during the propagator
initialization. If it is `OrbitStateVector`, the output is an array of `OrbitStateVector`
[SI] objects.

In this case, the user can pass the keyword argument `ntasks` to specify the number of
parallel tasks to propagate the orbit. If it is set to a number equal or lower than 1, the
function will propagate the orbit sequentially. If `ntasks` is omitted, it defaults to the
number of available threads.

For example, using the initialized propagator, we can propagate the orbit to midnight of
every day in January 2023 as follows:

```@repl usage
Propagators.propagate_to_epoch!(
    orbp,
    date_to_jd(2023, 1, 1, 0, 0, 0):date_to_jd(2023, 1, 31, 0, 0, 0)
)

Propagators.propagate_to_epoch!(
    orbp,
    DateTime("2023-01-01"):Dates.Day(1):DateTime("2023-01-31")
)
```

To obtain an array of `OrbitStateVector` as the output, we can use:

```@repl usage
Propagators.propagate_to_epoch!(
    orbp,
    date_to_jd(2023, 1, 1, 0, 0, 0):date_to_jd(2023, 1, 31, 0, 0, 0),
    OrbitStateVector
)

Propagators.propagate_to_epoch!(
    orbp,
    DateTime("2023-01-01"):Dates.Day(1):DateTime("2023-01-31"),
    OrbitStateVector
)
```

---

```julia
step!(orbp::OrbitPropagator, Δt::Number[, sink = Tuple])
```

Every time we propagate the orbit, the propagation instant is recorded inside the structure.
Hence, we can use the function `step!` to propagate the orbit by `Δt` [s] from the current
orbit epoch.

The return type depends on the parameter `sink`. If `sink` is `Tuple` (default), the
function returns two vectors: the position [m] and velocity [m / s]. Both are represented in
the same inertial reference frame used to describe the input data during the propagator
initialization. On the other hand, if `sink` is `OrbitStateVector`, it returns an object of
type `OrbitStateVector` [SI] with the same information.

For example, using the initialized propagator, we can advance the propagation by 60 s as
follows:

```@repl usage
orbp

Propagators.step!(orbp, 60)

orbp
```

To obtain an `OrbitStateVector` as the output, we can use:

```@repl usage
orbp

Propagators.step!(orbp, 60, OrbitStateVector)

orbp
```

## Simultaneous Initialization and Propagation

We can use the following functions to simultaneously initialize and propagate the orbit:

```julia
propagate([sink = Tuple, ]::Val{:propagator}, Δt::Number, args...; kwargs...)
propagate([sink = Tuple, ]::Val{:propagator}, p::Union{Dates.Period, Dates.CompundPeriod}, args...; kwargs...)
propagate_to_epoch([sink = Tuple, ]::Val{:propagator}, jd::Number, args...; kwargs...)
propagate_to_epoch([sink = Tuple, ]::Val{:propagator}, dt::DateTime, args...; kwargs...)
```

The symbol `:propagator`, the arguments `args...`, and the keywords `kwargs...` are the same
as described in the section related to the propagator initialization.

The two first functions propagates the orbit by `Δt` [s] or by the period defined by `p`
from the initial epoch, whereas the two last functions propagated the orbit until the
Julian Day `jd` [UTC] or until the epoch defined by the `DateTime` object `dt` [UTC]

The return type of both functions depends on the parameter `sink`. If `sink` is `Tuple`
(default), the function returns two vectors: the position [m] and velocity [m / s]. Both are
represented in the same inertial reference frame used to describe the input data during the
propagator initialization. On the other hand, if `sink` is `OrbitStateVector`, it returns an
object of type `OrbitStateVector` [SI] with the same information. Additionally to those
informations, the function also return the initialized orbit propagator.

For example, we can propagate by 60 s a set of mean elements using a J2 analytical
propagator as follows:

```@repl usage
r_i, v_i, orbp = Propagators.propagate(Val(:J2), 60, orb)

r_i, v_i, orbp = Propagators.propagate(Val(:J2), Dates.Second(60), orb)
```

To obtain an `OrbitStateVector` as the output, we can use:

```@repl usage
sv_i, orbp = Propagators.propagate(OrbitStateVector, Val(:J2), 60, orb)

sv_i, orbp = Propagators.propagate(OrbitStateVector, Val(:J2), Dates.Second(60), orb)
```

The following algorithm propagates a set of mean elements using a J4 analytical propagator
until the midnight of June 19, 2023:

```@repl usage
r_i, v_i, orbp = Propagators.propagate_to_epoch(
    Val(:J4),
    date_to_jd(2023, 6, 19, 0, 0, 0),
    orb
)

r_i, v_i, orbp = Propagators.propagate_to_epoch(Val(:J4), DateTime("2023-06-19"), orb)
```

To obtain an `OrbitStateVector` as the output, we can use:

```@repl usage
sv_i, orbp = Propagators.propagate_to_epoch(
    OrbitStateVector,
    Val(:J4),
    date_to_jd(2023, 6, 19, 0, 0, 0),
    orb
)

sv_i, orbp = Propagators.propagate_to_epoch(
    OrbitStateVector,
    Val(:J4),
    DateTime("2023-06-19"),
    orb
)
```

The API also supports propagating the orbit for multiple instants or epochs if we pass a
time vector as follows:

```julia
propagate([sink = Tuple, ]::Val{:propagator}, vt::AbstractVector, args...; kwargs...)
propagate([sink = Tuple, ]::Val{:propagator}, vp::AbstractVector{Union{Dates.Period, Dates.CompundPeriod}}, args...; kwargs...)
propagate_to_epoch([sink = Tuple, ]::Val{:propagator}, vjd::AbstractVector, args...; kwargs...)
propagate_to_epoch([sink = Tuple, ]::Val{:propagator}, vdt::DateTime, args...; kwargs...)
```

The symbol `:propagator`, the arguments `args...`, and the keywords `kwargs...` are the same
as described in the section related to the propagator initialization.

The two first functions propagate the orbit to each instant in the vector `vt` [s] or period
in the array `vp` from the initial epoch. The two last functions propagate the orbit to
each Julian Day in the vector `vjd` [UTC] or epoch defined by a `DateTime` object in the
array `vdt` [UTC].

The return type of both functions depends on the parameter `sink`. If `sink` is `Tuple`
(default), the output is two arrays containing the position [m] and velocity vectors [m /
s]. Both are represented in the same inertial reference frame used to describe the input
data during the propagator initialization. On the other hand, if `sink` is
`OrbitStateVector`, it returns an array of `OrbitStateVector` [SI] with the same
information. Additionally to those informations, the function also return the initialized
orbit propagator.

For example, we can propagate a set of mean elements every 5 min in an hour using the J2
analytical propagator as follows:

```@repl usage
vr_i, vv_i, orbp = Propagators.propagate(Val(:J2), 0:300:3600, orb)

vr_i, vv_i, orbp = Propagators.propagate(
    Val(:J2),
    Dates.Minute(0):Dates.Minute(5):Dates.Minute(60),
    orb
)
```

To obtain an array of `OrbitStateVector` as the output, we can use:

```@repl usage
vsv_i, orbp = Propagators.propagate(OrbitStateVector, Val(:J2), 0:300:3600, orb)

vsv_i, orbp = Propagators.propagate(
    OrbitStateVector,
    Val(:J2),
    Dates.Minute(0):Dates.Minute(5):Dates.Minute(60),
    orb
)
```

The following algorithm propagates a set of mean elements to midnight of every day in
January 2023 using the J4 analytical propagator:

```@repl usage
vr_i, vv_i, orbp = Propagators.propagate_to_epoch(
    Val(:J4),
    date_to_jd(2023, 1, 1, 0, 0, 0):date_to_jd(2023, 1, 31, 0, 0, 0),
    orb
)

vr_i, vv_i, orbp = Propagators.propagate_to_epoch(
    Val(:J4),
    DateTime("2023-01-01"):Dates.Day(1):DateTime("2023-01-31"),
    orb
)
```

To obtain an array of `OrbitStateVector` as the output, we can use:

```@repl usage
vsv_i, orbp = Propagators.propagate_to_epoch(
    OrbitStateVector,
    Val(:J4),
    date_to_jd(2023, 1, 1, 0, 0, 0):date_to_jd(2023, 1, 31, 0, 0, 0),
    orb
)

vsv_i, orbp = Propagators.propagate_to_epoch(
    OrbitStateVector,
    Val(:J4),
    DateTime("2023-01-01"):Dates.Day(1):DateTime("2023-01-31"),
    orb
)
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
