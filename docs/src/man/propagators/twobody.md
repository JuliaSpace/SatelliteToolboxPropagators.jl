Two Body Analytical Propagator
==============================

```@meta
CurrentModule = SatelliteToolboxPropagators
DocTestSetup = quote
    using SatelliteToolboxPropagators
end
```

The two-body analytical orbit propagator considers the Earth a perfect sphere with uniform
density. Hence, it propagates the orbit using the solution considering Newtonian gravity. It
has an extremely low precision but with a minimal computational burden. Thus, it is helpful
in some analysis that requires propagating the orbit many times for short periods.

## Algorithm

The algorithm implemented here is based on **[1]**.

Since we are considering a spherical Earth with uniform density, gravity points towards the
center of Earth. Thus, we propagate the orbit by updating the satellite mean anomaly since
all other Keplerian elements do not change. The equation to correct the mean anomaly is:

```math
M(t) = M_0 + \sqrt{\frac{\mu}{a_0^3}} \cdot \left(t - t_0\right)
```

where ``t_0`` is the initial mean elements' epoch, ``a_0`` is the mean semi-major axis, and
``\mu`` the Earth's standard gravitational parameter.

## Initialization

We can initialize the two-body analytical propagator with the following function:

```julia
function Propagators.init(Val(:TwoBody), orb₀::KeplerianElements; kwargs...)
```

which creates a two-body propagator structure [`OrbitPropagatorTwoBody`](@ref) with the mean
Keplerian elements `orb₀`. The following keyword selects the standard gravitational
parameter for the propagation algorithm:

- `μ::T`: Standard gravitational parameter of the central body [m³/s²].
    (**Default** = `tbc_m0`)
    
This package contains some pre-built gravitational parameters of the Earth for this
propagator:

| **Two-Body Propagator Constant** | **Description**                          | **Type**  |
|---------------------------------:|:-----------------------------------------|-----------|
|                         `tbc_m0` | Earth's standard gravitational parameter | `Float64` |
|                     `tbc_m0_f32` | Earth's standard gravitational parameter | `Float32` |

!!! note
    The type used in the propagation will be the same as used to define the gravitational
    constant `μ`.
    
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

julia> orbp = Propagators.init(Val(:TwoBody), orb)
OrbitPropagatorTwoBody{Float64, Float64}:
   Propagator name : Two Body Orbit Propagator
  Propagator epoch : 2023-01-01T00:00:00
  Last propagation : 2023-01-01T00:00:00
```

## References

- **[1]** **Vallado, D. A** (2013). *Fundamentals of Astrodynamics and Applications*. 4th
  ed. **Microcosm Press**, Hawthorn, CA, USA.
