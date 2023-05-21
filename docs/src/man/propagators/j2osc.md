J2 Osculating Analytical Orbit Propagator
=========================================

```@meta
CurrentModule = SatelliteToolboxPropagators
DocTestSetup = quote
    using SatelliteToolboxPropagators
end
```

This algorithm uses the [J2 propagator](@ref J2_Propagator) to obtain the secular effects of
the Keplerian elements caused only by the J2 term of the geopotential field. Afterward, it
adds short-term perturbations. This model is useful when fitting an orbit to a set of mean
elements for the J2 orbit propagator. Hence, we can use it to verify, for example, how close
a satellite is from a Sun-Synchronous orbit.

## Algorithm

The algorithm implemented here is based on **[1]**.

## Initialization

We can initialize the J2 osculating analytical orbit propagator with the following function:

```julia
function Propagators.init(Val(:J2osc), orb₀::KeplerianElements, dn_o2::Number = 0, ddn_o6::Number = 0; kwargs...)
```

which creates a J2 osculating propagator structure [`OrbitPropagatorJ2Osculating`](@ref)
with the mean Keplerian elements `orb₀`. It supports the following positional arguments:

- `orb₀::KeplerianElements`: Initial mean Keplerian elements [SI units].
- `dn_o2::Number`: First time derivative of the mean motion divided by two [rad/s^2].
    (**Default** = 0)
- `ddn_o6::Number`: Second time derivative of the mean motion divided by six [rad/s^3].
    (**Default** = 0)
    
The following keyword selects the gravitational constants for the propagation algorithm:

- `j2c::J2PropagatorConstants`: J2 orbit propagator constants (see
  [`J2PropagatorConstants`](@ref)). (**Default** = `j2c_egm2008`)

This package contains some pre-built propagation constants for this propagator:

| **J2 Propagator Constant** | **Description**                  | **Type**  |
|---------------------------:|:---------------------------------|:----------|
|              `j2c_egm2008` | EGM-2008 gravitational constants | `Float64` |
|          `j2c_egm2008_f32` | EGM-2008 gravitational constants | `Float32` |
|              `j2c_egm1996` | EGM-1996 gravitational constants | `Float64` |
|          `j2c_egm1996_f32` | EGM-1996 gravitational constants | `Float32` |
|                `j2c_jgm02` | JGM-02 gravitational constants   | `Float64` |
|            `j2c_jgm02_f32` | JGM-02 gravitational constants   | `Float32` |
|                `j2c_jgm03` | JGM-03 gravitational constants   | `Float64` |
|            `j2c_jgm03_f32` | JGM-03 gravitational constants   | `Float32` |

!!! note
    The type used in the propagation will be the same as used to define the constants in the
    structure `j2c`.
    
```jldoctest
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

julia> orbp = Propagators.init(Val(:J2osc), orb)
OrbitPropagatorJ2Osculating{Float64, Float64}:
   Propagator name : J2 Osculating Orbit Propagator
  Propagator epoch : 2023-01-01T00:00:00
  Last propagation : 2023-01-01T00:00:00
```

## References

- **[1]** **Vallado, D. A** (2013). *Fundamentals of Astrodynamics and Applications*. 4th
  ed. **Microcosm Press**, Hawthorn, CA, USA.
