Usage
=====

```@meta
CurrentModule = SatelliteToolboxPropagators
DocTestSetup = quote
    using SatelliteToolboxPropagators
end
```

All the propagators can be accessed using the available API, which allows to initialize and
propagate the orbit.

All the API function are available inside the module [`Propagators`](@ref) that is exported
by this package.

## Initialization

We can initialize an orbit propagator using the function:

```julia
function Propagators.init(::Val{:propagator}, args...; kwargs...)
```

It initializes a propagator of type `:propagator` using the arguments `args...` and keywords
`kwargs...` supported by it. Currently, the following algorithms are available:

| **Propagator Name**                       | **Symbol** |
|------------------------------------------:|:-----------|
| J2 analytical orbit propagator            | `:J2`      |
| J2 osculating analytical orbit propagator | `:J2osc`   |
| J4 analytical orbit propagator            | `:J4`      |
| SGP4/SDP4 orbit propagator                | `:SGP4`    |

See the documentation of each algorithm to verify the supported arguments and keywords.

For example, a J2 analytical orbit propagator can be initialized with a set of mean elements
as follows:

```jldoctest J2
julia> orb = KeplerianElements(
           date_to_jd(2023, 1, 1, 0, 0, 0),
           7190.982e3,
           0.001111,
           98.405 |> deg2rad,
           100    |> deg2rad,
           90     |> deg2rad,
           19     |> deg2rad
       )
KeplerianElements{Float64, Float64}:
           Epoch :    2.45995e6 (2023-01-01T00:00:00)
 Semi-major axis : 7190.98     km
    Eccentricity :    0.001111
     Inclination :   98.405    °
            RAAN :  100.0      °
 Arg. of Perigee :   90.0      °
    True Anomaly :   19.0      °

julia> orbp = Propagators.init(Val(:J2), orb)
OrbitPropagatorJ2{Float64, Float64}:
   Propagator name : J2 Orbit Propagator
  Propagator epoch : 2023-01-01T00:00:00
  Last propagation : 2023-01-01T00:00:00
```

We also have the function:

```julia
function init!(orbp::OrbitPropagator, args...; kwargs...)
```

that initializes the propagator `orbp` in-place using the arguments `args...` and keywords
`kwargs...`. This function allows to reduce the memory allocation when performing multiple
initialization. However, notice that this function is **optional** in the propagator API.
Hence, some algorithms might not support it.

## Propagate the Orbit

After the initialization, we can propagate the orbit using some functions as follows.

```julia
function propagate!(orbp::OrbitPropagator{Tepoch, T}, t::Number) where {Tepoch, T}
```

This function propagates the orbit using `orbp` by `t` [s] from the initial orbit epoch. It
returns two vectors: the position [m] and velocity [m/s]. Both are represented in the same
inertial reference frame used to describe the input data during the propagator
initialization.

For example, using the initialized propagator, we can propagate the orbit by 6000s as
follows:

```jldoctest J2
julia> Propagators.propagate!(orbp, 6000)
([1.3152752189693751e6, -1.5933886600051078e6, 6.8797020668377e6], [993.0153602143952, -7152.788090016819, -1844.2917380794433])
```

---

```julia
function propagate_to_epoch!(orbp::OrbitPropagator{Tepoch, T}, jd::Number) where {Tepoch, T}
```

This function propagates the orbit using `orbp` until the epoch `jd` [Julian Day]. It
returns two vectors: the position [m] and velocity [m/s]. Both are represented in the same
inertial reference frame used to describe the input data during the propagator
initialization.

For example, using the initialized propagator, we can propagate the orbit to
2023-01-02T00:00:00 as follows:

```jldoctest J2
julia> Propagators.propagate_to_epoch!(orbp, date_to_jd(2023, 1, 2, 0, 0, 0))
([1.2007441457621562e6, -7.014291058654364e6, -1.0443532029853627e6], [-1262.9147609767083, 860.1551065186032, -7285.029142259778])
```

---

```julia
function step!(orbp::OrbitPropagator{Tepoch, T}, Δt::Number) where {Tepoch, T}
```

Every time we propagate the orbit, the propagation instant is recorded inside the structure.
Hence, we can use the function `step!` to propagate the orbit by `Δt` [s] from the current
orbit epoch. It returns two vectors: the position [m] and velocity [m/s]. Both are
represented in the same inertial reference frame used to describe the input data during the
propagator initialization.

For example, using the initialized propagator, we can advance the propagation by 60 s as
follows:

```jldoctest J2
julia> orbp
OrbitPropagatorJ2{Float64, Float64}:
   Propagator name : J2 Orbit Propagator
  Propagator epoch : 2023-01-01T00:00:00
  Last propagation : 2023-01-02T00:00:00

julia> Propagators.step!(orbp, 60)
([1.1228778966719273e6, -6.949273896376436e6, -1.4786561473080558e6], [-1337.5317183567313, 1308.4670676524293, -7204.023280682461])

julia> orbp
OrbitPropagatorJ2{Float64, Float64}:
   Propagator name : J2 Orbit Propagator
  Propagator epoch : 2023-01-01T00:00:00
  Last propagation : 2023-01-02T00:01:00
```

## Simultaneous Initialization and Propagation

We can use the following functions to simultaneously initialize and propagate the orbit:

```julia
function propagate(::Val{:propagator}, Δt::Number, args...; kwargs...)
function propagate_to_epoch(::Val{:propagator}, jd::Number, args...; kwargs...)
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

```jldoctest J2
julia> Propagators.propagate(Val(:J2), 60, orb)
([1.433577581999473e6, -2.5460150747451796e6, 6.562531225282293e6], [784.1942683097858, -6851.527297258612, -2825.9666374272706], J2 Orbit Propagator (Epoch = 2023-01-01T00:00:00, Δt = 60.0 s))
```

## Helpers

We have the following functions that provide some useful information related to the
propagators.

```julia
function epoch(orbp::OrbitPropagator{Tepoch, T}) where {Tepoch<:Number, T<:Number}
```

It returns the initial elements' epoch of the propagator `orbp` [JD].

```jldoctest J2
julia> Propagators.epoch(orbp)
2.4599455e6
```

---

```julia
function last_instant(orbp::OrbitPropagator{Tepoch, T}) where {Tepoch<:Number, T<:Number}
```

It returns the last propagation instant [s] measured from the input elements' epoch.

```jldoctest J2
julia> Propagators.last_instant(orbp)
86460.0
```

---

```julia
function name(orbp::OrbitPropagator)
```

It returns the propagator name.

```jldoctest J2
julia> Propagators.name(orbp)
"J2 Orbit Propagator"
```

---

```julia
function mean_elements(orbp::OrbitPropagator)
```

It returns the mean elements using the structure `KeplerianElements` of the latest
propagation performed by `orbp`. Notice that this is an optional funciton in the
propagators' API. If a propagator does not support it, this function returns `nothing`.

```jldoctest J2
julia> Propagators.mean_elements(orbp)
KeplerianElements{Float64, Float64}:
           Epoch :    2.45995e6 (2023-01-02T00:01:00)
 Semi-major axis : 7190.98     km
    Eccentricity :    0.001111
     Inclination :   98.405    °
            RAAN :  100.957    °
 Arg. of Perigee :   87.0755   °
    True Anomaly :  104.918    °
```
