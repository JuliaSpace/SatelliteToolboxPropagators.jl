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
function Propagators.init(Val(:J2osc), orb₀::KeplerianElements; kwargs...)
```

which creates a J2 osculating propagator structure [`OrbitPropagatorJ2Osculating`](@ref)
with the mean Keplerian elements `orb₀`.
    
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

## Fitting Mean Elements

We can use the function:

```julia
function Propagators.fit_mean_elements(::Val{:J2osc}, vjd::AbstractVector{Tjd}, vr_i::AbstractVector{Tv}, vv_i::AbstractVector{Tv}; kwargs...) -> KeplerianElements{Float64, Float64}, SMatrix{6, 6, Float64}
```

to fit a set of mean Keplerian elements for the J2 osculating orbit propagator using the
osculating elements represented by a set of position vectors `vr_i` [m] and a set of
velocity vectors `vv_i` [m / s] represented in an inertial reference frame at instants in
the array `vjd` [Julian Day].

It returns the fitted Keplerian elements and the final covariance matrix of the least-square
algorithm.

!!! note
    This algorithm version will allocate a new J2 propagator with the default constants
    `j2c_egm2008`. If another set of constants are required, use the function
    [`Propagators.fit_mean_elements!`](@ref) instead.
    
The following keywords are available to configure the fitting process:

- `atol::Number`: Tolerance for the residue absolute value. If the residue is lower than
    `atol` at any iteration, the computation loop stops. (**Default** = 2e-4)
- `rtol::Number`: Tolerance for the relative difference between the residues. If the
    relative difference between the residues in two consecutive iterations is lower than
    `rtol`, the computation loop stops. (**Default** = 2e-4)
- `initial_guess::Union{Nothing, KeplerianElements}`: Initial guess for the mean elements
    fitting process. If it is `nothing`, the algorithm will obtain an initial estimate from
    the osculating elements in `vr_i` and `vv_i`. (**Default** = nothing)
- `jacobian_perturbation::Number`: Initial state perturbation to compute the
    finite-difference when calculating the Jacobian matrix. (**Default** = 1e-3)
- `jacobian_perturbation_tol::Number`: Tolerance to accept the perturbation when calculating
    the Jacobian matrix. If the computed perturbation is lower than
    `jacobian_perturbation_tol`, we increase it until it absolute value is higher than
    `jacobian_perturbation_tol`. (**Default** = 1e-7)
- `max_iterations::Int`: Maximum number of iterations allowed for the least-square fitting.
    (**Default** = 50)
- `mean_elements_epoch::Number`: Epoch for the fitted mean elements.
    (**Default** = vjd[end])
- `verbose::Bool`: If `true`, the algorithm prints debugging information to `stdout`.
    (**Default** = true)
- `weight_vector::AbstractVector`: Vector with the measurements weights for the least-square
    algorithm. We assemble the weight matrix `W` as a diagonal matrix with the elements in
    `weight_vector` at its diagonal. (**Default** = `@SVector(ones(Bool, 6))`)
    
```julia
julia> vr_i = [
           [-6792.402703741442, 2192.6458461287293, 0.18851758695295118] .* 1000,
           [-6357.88873265975, 2391.9476768911686, 2181.838771262736] .* 1000
       ];

julia> vv_i = [
           [0.3445760107690598, 1.0395135806993514, 7.393686131436984] .* 1000,
           [2.5285015912807003, 0.27812476784300005, 7.030323100703928] .* 1000
       ];

julia> vjd = [
           2.46002818657856e6,
           2.460028190050782e6
       ];

julia> orb, P = Propagators.fit_mean_elements(Val(:J2), vjd, vr_i, vv_i)
ACTION:   Fitting the mean elements for the J2 propagator.
           Iteration        Position RMSE        Velocity RMSE           Total RMSE       RMSE Variation
                                     [km]             [km / s]                  [ ]
PROGRESS:          3          6.15798e-05           0.00975024              9.75043         -0.000203903 %

(KeplerianElements{Float64, Float64}: Epoch = 2.46003e6 (2023-03-24T16:33:40.388), [0.9999771881514358 3.332831176445557e-7 … -9.720588392821589e-5 -7.124383305247769e-5; 3.3327603024689727e-7 0.9999779181940862 … 0.0032646452155289567 1.7018081890699746e-5; … ; -9.720588394942824e-5 0.0032646452155281466 … 2.2082152640332316e-5 1.794499832376423e-8; -7.12438330606674e-5 1.7018081892440988e-5 … 1.7944998330774395e-8 2.1724736479160325e-5])

julia> orb
KeplerianElements{Float64, Float64}:
           Epoch :    2.46003e6 (2023-03-24T16:33:40.388)
 Semi-major axis : 7155.99       km
    Eccentricity :    0.00307382
     Inclination :   98.4357     °
            RAAN :  162.113      °
 Arg. of Perigee :   33.0518     °
    True Anomaly :  344.956      °
```

## References

- **[1]** **Vallado, D. A** (2013). *Fundamentals of Astrodynamics and Applications*. 4th
  ed. **Microcosm Press**, Hawthorn, CA, USA.
