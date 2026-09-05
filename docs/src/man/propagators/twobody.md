# Two Body Analytical Propagator

```@meta
CurrentModule = SatelliteToolboxPropagators
```

```@setup tb
using SatelliteToolboxPropagators
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
Propagators.init(Val(:TwoBody), orb₀::KeplerianElements; kwargs...) -> OrbitPropagatorTwoBody
```

which creates a two-body propagator structure [`OrbitPropagatorTwoBody`](@ref) with the mean
Keplerian elements `orb₀`. The following keyword selects the standard gravitational
parameter for the propagation algorithm:

- `m0::T`: Standard gravitational parameter of the central body [m³/s²].
    (**Default** = `TBC_M0`)

This package contains some pre-built gravitational parameters of the Earth for this
propagator:

| **Two-Body Propagator Constant** | **Description**                          | **Type**  |
|---------------------------------:|:-----------------------------------------|-----------|
|                         `TBC_M0` | Earth's standard gravitational parameter | `Float64` |
|                     `TBC_M0_F32` | Earth's standard gravitational parameter | `Float32` |

!!! note

    The type used in the propagation will be the same as used to define the gravitational
    constant `μ`.

```@repl tb
orb = KeplerianElements(
           date_to_jd(2023, 1, 1, 0, 0, 0),
           7190.982e3,
           0.001111,
           98.405 |> deg2rad,
           100    |> deg2rad,
           90     |> deg2rad,
           19     |> deg2rad
       )

orbp = Propagators.init(Val(:TwoBody), orb)
```

## Fitting Mean Elements

We can use the function:

```julia
Propagators.fit_mean_elements(::Val{:TwoBody}, vjd::AbstractVector{Tjd}, vr_i::AbstractVector{Tv}, vv_i::AbstractVector{Tv}; kwargs...) -> KeplerianElements{MeanAnomaly, Float64, Float64}, SMatrix{6, 6, Float64}
```

to fit a set of mean Keplerian elements for the two-body orbit propagator using the
osculating elements represented by a set of position vectors `vr_i` [m] and a set of
velocity vectors `vv_i` [m / s] represented in an inertial reference frame at instants in
the array `vjd` [Julian Day].

It returns the fitted Keplerian elements and the final covariance matrix of the least-square
algorithm.

!!! note

    This algorithm version will allocate a new two-body propagator with the default
    gravitational parameter `TBC_M0`. If another value is required, use the function
    [`Propagators.fit_mean_elements!`](@ref) instead.

The following keywords are available to configure the fitting process:

- `atol::Number`: Tolerance for the residual absolute value. If the residual is lower than
    `atol` at any iteration, the computation loop stops.
    (**Default** = 2e-4)
- `rtol::Number`: Tolerance for the relative difference between the residuals. If the
    relative difference between the residuals in two consecutive iterations is lower than
    `rtol`, the computation loop stops.
    (**Default** = 2e-4)
- `initial_guess::Union{Nothing, KeplerianElements}`: Initial guess for the mean elements
    fitting process. If it is `nothing`, the algorithm will obtain an initial estimate from
    the osculating elements in `vr_i` and `vv_i`.
    (**Default** = nothing)
- `jacobian_method::AbstractJacobianMethod`: Method used to compute the Jacobian matrix. It
    can be `FiniteDiffJacobian()` for finite differences or `ForwardDiffJacobian()` for
    `ForwardDiff.jl` automatic differentiation.
    (**Default** = `FiniteDiffJacobian()`)
- `jacobian_perturbation::Number`: Initial state perturbation to compute the
    finite-difference when calculating the Jacobian matrix.
    (**Default** = 1e-3)
- `jacobian_perturbation_tol::Number`: Tolerance to accept the perturbation when calculating
    the Jacobian matrix. If the computed perturbation is lower than
    `jacobian_perturbation_tol`, we increase it until its absolute value is higher than
    `jacobian_perturbation_tol`.
    (**Default** = 1e-7)
- `max_iterations::Int`: Maximum number of iterations allowed for the least-square fitting.
    (**Default** = 50)
- `mean_elements_epoch::Number`: Epoch for the fitted mean elements.
    (**Default** = vjd[end])
- `verbose::Bool`: If `true`, the algorithm prints debugging information to `stdout`.
    (**Default** = true)
- `weight_vector::AbstractVector`: Vector with the measurements weights for the least-square
    algorithm. We assemble the weight matrix `W` as a diagonal matrix with the elements in
    `weight_vector` at its diagonal.
    (**Default** = `@SVector(ones(Bool, 6))`)

```@repl tb
vr_i = [
    [-6792.402703741442, 2192.6458461287293, 0.18851758695295118] .* 1000,
    [-6357.88873265975, 2391.9476768911686, 2181.838771262736] .* 1000
];

vv_i = [
    [0.3445760107690598, 1.0395135806993514, 7.393686131436984] .* 1000,
    [2.5285015912807003, 0.27812476784300005, 7.030323100703928] .* 1000
];

vjd = [
    2.46002818657856e6,
    2.460028190050782e6
];

orb, P = Propagators.fit_mean_elements(Val(:TwoBody), vjd, vr_i, vv_i)

orb
```

## References

- **[1]** **Vallado, D. A** (2013). *Fundamentals of Astrodynamics and Applications*. 4th
  ed. **Microcosm Press**, Hawthorn, CA, USA.
