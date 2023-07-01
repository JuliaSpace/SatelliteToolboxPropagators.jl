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
              Eccentricity :   0.00012470
               Inclination :  98.43040000 deg
                      RAAN : 162.10970000 deg
       Argument of perigee : 136.20170000 deg
              Mean anomaly : 223.92830000 deg
           Mean motion (n) :  14.40814394 revs / day
         Revolution number : 10865
                        B* :      4.3e-05 1 / er
                     ṅ / 2 :     -4.4e-07 rev / day²
                     n̈ / 6 :        1e-09 rev / day³

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

## Fitting TLEs

We can use the function:

```julia
function Propagators.fit_mean_elements(::Val{:SGP4}, vjd::AbstractVector{Tjd}, vr_teme::AbstractVector{Tv}, vv_teme::AbstractVector{Tv}; kwargs...) -> KeplerianElements{Float64, Float64}, SMatrix{6, 6, Float64}
```

to fit a Two-Line Element set (`TLE`) for the SGP4 orbit propagator using the osculating
elements represented by a set of position vectors `vr_teme` [m] and a set of velocity
vectors `vv_teme` [m / s] represented in the True-Equator, Mean-Equinox reference frame
(TEME) at instants in the array `vjd` [Julian Day].

It returns the fitted TLE and the final covariance matrix of the least-square algorithm.

This algorithm was based on **[4]**.

!!! notes
    This algorithm version will allocate a new SGP4 propagator with the default constants
    `sgp4c_wgs84`. If another set of constants are required, use the function
    [`Propagators.fit_mean_elements!`](@ref) instead.
    
The following keywords are available to configure the fitting process:

- `atol::Number`: Tolerance for the residue absolute value. If the residue is lower than
    `atol` at any iteration, the computation loop stops. (**Default** = 2e-4)
- `rtol::Number`: Tolerance for the relative difference between the residues. If the
    relative difference between the residues in two consecutive iterations is lower than
    `rtol`, the computation loop stops. (**Default** = 2e-4)
- `estimate_bstar::Bool`: If `true`, the algorithm will try to estimate the B* parameter.
    Otherwise, it will be set to 0 or to the value in initial guess (see section  **Initial
    Guess**). (**Default** = true)
- `initial_guess::Union{Nothing, AbstractVector, TLE}`: Initial guess for the TLE fitting
    process. If it is `nothing`, the algorithm will obtain an initial estimate from the
    osculating elements in `vr_teme` and `vv_teme`. For more information, see the section
    **Initial Guess**. (**Default** = nothing)
- `jacobian_perturbation::Number`: Initial state perturbation to compute the
    finite-difference when calculating the Jacobian matrix. (**Default** = 1e-3)
- `jacobian_perturbation_tol::Number`: Tolerance to accept the perturbation when calculating
    the Jacobian matrix. If the computed perturbation is lower than
    `jacobian_perturbation_tol`, we increase it until it absolute value is higher than
    `jacobian_perturbation_tol`. (**Default** = 1e-7)
- `max_iterations::Int`: Maximum number of iterations allowed for the least-square fitting.
    (**Default** = 50)
- `mean_elements_epoch::Number`: Epoch for the fitted TLE. (**Default** = vjd[end])
- `verbose::Bool`: If `true`, the algorithm prints debugging information to `stdout`.
    (**Default** = true)
- `weight_vector::AbstractVector`: Vector with the measurements weights for the least-square
    algorithm. We assemble the weight matrix `W` as a diagonal matrix with the elements in
    `weight_vector` at its diagonal. (**Default** = `@SVector(ones(Bool, 6))`)
- `classification::Char`: Satellite classification character for the output TLE.
    (**Default** = 'U')
- `element_set_number::Int`: Element set number for the output TLE. (**Default** = 0)
- `international_designator::String`: International designator string for the output TLE.
    (**Default** = "999999")
- `name::String`: Satellite name for the output TLE. (**Default** = "UNDEFINED")
- `revolution_number::Int`: Revolution number for the output TLE. (**Default** = 0)
- `satellite_number::Int`: Satellite number for the output TLE. (**Default** = 9999)
    
```julia
julia> vr_teme = [
           [-6792.402703741442, 2192.6458461287293, 0.18851758695295118] .* 1000,
           [-6357.88873265975, 2391.9476768911686, 2181.838771262736] .* 1000
       ];

julia> vv_teme = [
           [0.3445760107690598, 1.0395135806993514, 7.393686131436984] .* 1000,
           [2.5285015912807003, 0.27812476784300005, 7.030323100703928] .* 1000
       ];

julia> vjd = [
           2.46002818657856e6,
           2.460028190050782e6
       ];

julia> tle, P = Propagators.fit_mean_elements(Val(:SGP4), vjd, vr_teme, vv_teme; estimate_bstar = false)
ACTION:   Fitting the TLE.
           Iteration        Position RMSE        Velocity RMSE           Total RMSE       RMSE Variation
                                     [km]             [km / s]                  [ ]
PROGRESS:          3          5.80079e-09          4.38304e-07          4.38343e-07             -99.9514 %

(TLE: UNDEFINED (Epoch = 2023-03-24T16:33:40.388), [1.085542785763592 -0.02468955108588304 … -0.07505082051267041 -5691.217934788844; -0.024689551096344593 1.0180284581267773 … 0.023015858293558573 1745.6307022638805; … ; -0.07505082051315085 0.023015858293655808 … 0.06515845129137222 4946.150014964575; -5691.217934826224 1745.630702271188 … 4946.15001496459 3.7558624425042915e8])

julia> tle
TLE:
                      Name : UNDEFINED
          Satellite number : 9999
  International designator : 999999
        Epoch (Year / Day) : 23 /  83.69005079 (2023-03-24T16:33:40.388)
        Element set number : 0
              Eccentricity :   0.00012463
               Inclination :  98.43040000 deg
                      RAAN : 162.11312239 deg
       Argument of perigee : 136.15040637 deg
              Mean anomaly : 241.97934028 deg
           Mean motion (n) :  14.40814157 revs / day
         Revolution number : 0
                        B* :            0 1 / er
                     ṅ / 2 :            0 rev / day²
                     n̈ / 6 :            0 rev / day³
```

## References

- **[1]** **Hoots, F. R., Roehrich, R. L** (1980). *Models for Propagation of NORAD Elements
  Set*. **Spacetrack Report No. 3**.
- **[2]** **Vallado, D. A., Crawford, P., Hujsak, R., Kelso, T. S** (2006). *Revisiting
  Spacetrack Report #3: Rev1*. **AIAA**.
- **[3]** SGP4 Source code of [STRF](https://github.com/cbassa/strf), which the C code was
  converted by Paul. S. Crawford and Andrew R. Brooks.
- **[4]** **Vallado, D. A., Crawford, P** (2008). *SGP4 Orbit Determination*. **AIAA**.
