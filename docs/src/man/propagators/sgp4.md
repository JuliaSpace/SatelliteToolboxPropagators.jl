SGP4/SDP4
=========

```@meta
CurrentModule = SatelliteToolboxPropagators
DocTestSetup = quote
    using SatelliteToolboxPropagators
end
```

This package implements the interface to the
[SGP4/SDP4](https://en.wikipedia.org/wiki/Simplified_perturbations_models) propagator
provided by
[**SatelliteToolboxSgp4.jl**](https://github.com/JuliaSpace/SatelliteToolboxSgp4.jl).

## Algorithm

The SGP4/SDP4 implementation was built using **[1, 2, 3]**.

## Initialization

We must initialize the SGP4/SDP4 propagator with a [two-line element set
(TLE)](https://en.wikipedia.org/wiki/Two-line_element_set) using the following function:

```julia
function Propagators.init(Val(:SGP4), tle::TLE; kwargs...)
```

which creates a SGP4/SDP4 propagator structure [`OrbitPropagatorSgp4`](@ref) with the `tle`.
The following keyword selects the gravitational constants for the propagation algorithm:

- `sgp4c::Sgp4Constants`: SGP4 orbit propagator constants (see [`Sgp4Constants`](@ref)).
    (**Default** = `sgp4c_wgs84`)
    
The package [**SatelliteToolboxSgp4.jl**] contains some pre-build constants for this
propagator:

| **SGP4/SDP4 Propagator Constants** | **Description**           | **Type**  |
|-----------------------------------:|:--------------------------|:----------|
|                      `sgp4c_wgs84` | Constants based on WGS-84 | `Float64` |
|                  `sgp4c_wgs84_f32` | Constants based on WGS-84 | `Float32` |
|                      `sgp4c_wgs72` | Constants based on WGS-72 | `Float64` |
|                  `sgp4c_wgs72_f32` | Constants based on WGS-72 | `Float32` |

!!! note
    The type used in the propagation will be the same as used to define the constants in the
    structure `sgp4c`.
    
!!! note
    The package
    [**SatelliteToolboxTle.jl**](https://github.com/JuliaSpace/SatelliteToolboxTle.jl)
    defines the type `TLE`, which is re-exported here. It contains some useful
    functionalities, such as TLE fetching from online services. For more information, refer
    to the [package
    documentation](https://juliaspace.github.io/SatelliteToolboxTle.jl/stable/).
    
```jldoctest
julia> tle = tle"""
       AMAZONIA 1
       1 47699U 21015A   23083.68657856 -.00000044  10000-8  43000-4 0  9990
       2 47699  98.4304 162.1097 0001247 136.2017 223.9283 14.40814394108652"""
TLE:
                     Name : AMAZONIA 1
         Satellite number : 47699
 International designator : 21015A
       Epoch (Year / Day) : 23 /  83.68657856 (2023-03-24T16:28:40.388)
       Element set number : 999
             Eccentricity :   0.00012470 deg
              Inclination :  98.43040000 deg
                     RAAN : 162.10970000 deg
      Argument of perigee : 136.20170000 deg
             Mean anomaly : 223.92830000 deg
          Mean motion (n) :  14.40814394 revs/day
        Revolution number : 10865
                       B* : 0.000043 1/[er]
                    ṅ / 2 : -0.000000 rev/day²
                    n̈ / 6 : 0.000000 rev/day³

julia> Propagators.init(Val(:SGP4), tle)
OrbitPropagatorSgp4{Float64, Float64}:
   Propagator name : SGP4 Orbit Propagator
  Propagator epoch : 2023-03-24T16:28:40.388
  Last propagation : 2023-03-24T16:28:40.388
```

We also support initializing an SGP4/SDP4 propagator passing the TLE information in
individual terms:

```
function Propagators.init(Val(:SGP4), epoch::Number, n₀::Number, e₀::Number, i₀::Number, Ω₀::Number, ω₀::Number, M₀::Number, bstar::Number; kwargs...)
```

where:

- `epoch::Number`: Epoch of the orbital elements [Julian Day].
- `n₀::Number`: SGP type "mean" mean motion at epoch [rad/s].
- `e₀::Number`: "Mean" eccentricity at epoch.
- `i₀::Number`: "Mean" inclination at epoch [rad].
- `Ω₀::Number`: "Mean" longitude of the ascending node at epoch [rad].
- `ω₀::Number`: "Mean" argument of perigee at epoch [rad].
- `M₀::Number`: "Mean" mean anomaly at epoch [rad].
- `bstar::Number`: Drag parameter (B*).

The keywords `kwargs...` are the same as in the first version of this function.

## References

- **[1]** **Hoots, F. R., Roehrich, R. L** (1980). *Models for Propagation of NORAD Elements
  Set*. **Spacetrack Report No. 3**.
- **[2]** **Vallado, D. A., Crawford, P., Hujsak, R., Kelso, T. S** (2006). *Revisiting
  Spacetrack Report #3: Rev1*. **AIAA**.
- **[3]** SGP4 Source code of [STRF](https://github.com/cbassa/strf), which the C code was
  converted by Paul. S. Crawford and Andrew R. Brooks.
