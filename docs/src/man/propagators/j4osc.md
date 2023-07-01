J4 Osculating Analytical Orbit Propagator
=========================================

```@meta
CurrentModule = SatelliteToolboxPropagators
DocTestSetup = quote
    using SatelliteToolboxPropagators
end
```

This algorithm uses the [J4 propagator](@ref J4_Propagator) to obtain the secular effects of
the Keplerian elements caused by the terms ``J_2``, ``J_2^2``, and ``J_4``of the
geopotential field. Afterward, it adds short-term perturbations using only the ``J_2`` term.
This model is useful when fitting an orbit to a set of mean elements for the J4 orbit
propagator.

## Algorithm

The algorithm implemented here is based on **[1]**.

## Initialization

We can initialize the J4 osculating analytical orbit propagator with the following function:

```julia
function Propagators.init(Val(:J4osc), orb₀::KeplerianElements, dn_o2::Number = 0, ddn_o6::Number = 0; kwargs...)
```

which creates a J4 osculating propagator structure [`OrbitPropagatorJ4Osculating`](@ref)
with the mean Keplerian elements `orb₀`. It supports the following positional arguments:

- `orb₀::KeplerianElements`: Initial mean Keplerian elements [SI units].
- `dn_o2::Number`: First time derivative of the mean motion divided by two [rad/s^2].
    (**Default** = 0)
- `ddn_o6::Number`: Second time derivative of the mean motion divided by six [rad/s^3].
    (**Default** = 0)
    
The following keyword selects the gravitational constants for the propagation algorithm:

- `j4c::J4PropagatorConstants`: J4 orbit propagator constants (see
  [`J4PropagatorConstants`](@ref)). (**Default** = `j4c_egm2008`)

This package contains some pre-built propagation constants for this propagator:

| **J4 Propagator Constant** | **Description**                  | **Type**  |
|---------------------------:|:---------------------------------|:----------|
|              `j4c_egm2008` | EGM-2008 gravitational constants | `Float64` |
|          `j4c_egm2008_f32` | EGM-2008 gravitational constants | `Float32` |
|              `j4c_egm1996` | EGM-1996 gravitational constants | `Float64` |
|          `j4c_egm1996_f32` | EGM-1996 gravitational constants | `Float32` |
|                `j4c_jgm02` | JGM-02 gravitational constants   | `Float64` |
|            `j4c_jgm02_f32` | JGM-02 gravitational constants   | `Float32` |
|                `j4c_jgm03` | JGM-03 gravitational constants   | `Float64` |
|            `j4c_jgm03_f32` | JGM-03 gravitational constants   | `Float32` |

!!! note
    The type used in the propagation will be the same as used to define the constants in the
    structure `j4c`.
    
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

julia> orbp = Propagators.init(Val(:J4osc), orb)
OrbitPropagatorJ4Osculating{Float64, Float64}:
   Propagator name : J4 Osculating Orbit Propagator
  Propagator epoch : 2023-01-01T00:00:00
  Last propagation : 2023-01-01T00:00:00
```

## Fitting Mean Elements

We can use the function:

```julia
function Propagators.fit_mean_elements(::Val{:J4osc}, vjd::AbstractVector{Tjd}, vr_i::AbstractVector{Tv}, vv_i::AbstractVector{Tv}; kwargs...) -> KeplerianElements{Float64, Float64}, SMatrix{6, 6, Float64}
```

to fit a set of mean Keplerian elements for the J4 osculating orbit propagator using the
osculating elements represented by a set of position vectors `vr_i` [m] and a set of
velocity vectors `vv_i` [m / s] represented in an inertial reference frame at instants in
the array `vjd` [Julian Day].

It returns the fitted Keplerian elements and the final covariance matrix of the least-square
algorithm.

!!! note
    This algorithm version will allocate a new J4 propagator with the default constants
    `j4c_egm2008`. If another set of constants are required, use the function
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

julia> orb, P = Propagators.fit_mean_elements(Val(:J4), vjd, vr_i, vv_i)
ACTION:   Fitting the mean elements for the J4 propagator.
           Iteration        Position RMSE        Velocity RMSE           Total RMSE       RMSE Variation
                                     [km]             [km / s]                  [ ]
PROGRESS:          3          6.13971e-05           0.00972202              9.72221         -0.000197872 %

(KeplerianElements{Float64, Float64}: Epoch = 2.46003e6 (2023-03-24T16:33:40.388), [0.9999771881617059 3.332632038373772e-7 … -9.720779655127568e-5 -7.122976863467179e-5; 3.332604662785857e-7 0.9999779183950048 … 0.0032646363435036604 1.7058462259742015e-5; … ; -9.720779655821982e-5 0.0032646363435034492 … 2.208202282175653e-5 1.8337551711581296e-8; -7.122976863572516e-5 1.7058462251039033e-5 … 1.8337551683284122e-8 2.172632163462607e-5])

julia> orb
KeplerianElements{Float64, Float64}:
           Epoch :    2.46003e6 (2023-03-24T16:33:40.388)
 Semi-major axis : 7155.94       km
    Eccentricity :    0.00306707
     Inclination :   98.4357     °
            RAAN :  162.113      °
 Arg. of Perigee :   33.108      °
    True Anomaly :  344.9        °
```

## References

- **[1]** **Vallado, D. A** (2013). *Fundamentals of Astrodynamics and Applications*. 4th
  ed. **Microcosm Press**, Hawthorn, CA, USA.
