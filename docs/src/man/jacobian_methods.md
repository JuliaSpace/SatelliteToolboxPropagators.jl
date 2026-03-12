# Jacobian Methods: ForwardDiff vs Finite Differences

```@meta
CurrentModule = SatelliteToolboxPropagators
```

All mean-element fitting functions (`fit_j2_mean_elements`, `fit_j2osc_mean_elements`,
`fit_j4_mean_elements`, `fit_j4osc_mean_elements`, and the SGP4 counterpart
`fit_sgp4_tle`) use an iterative least-squares algorithm that requires Jacobian evaluation
at every iteration. Two methods are available, selectable via the `jacobian_method` keyword:

| Method | Type | Description |
|:-------|:-----|:------------|
| Finite differences | `FiniteDiffJacobian()` | Perturbs each state element and evaluates the propagator (default) |
| Forward-mode AD | `ForwardDiffJacobian()` | Uses `ForwardDiff.jl` dual numbers for exact derivatives |

## Quick Example

```julia
using SatelliteToolboxPropagators

orb_input = KeplerianElements(
    date_to_jd(2023, 1, 1),
    7130.982e3, 0.001111,
    98.405 |> deg2rad, 90.0 |> deg2rad,
    200.0  |> deg2rad, 45.0 |> deg2rad,
)

orbp = Propagators.init(Val(:J2), orb_input)
ret  = Propagators.propagate!.(orbp, 0:10:12_000)
vr_i = first.(ret)
vv_i = last.(ret)
vjd  = Propagators.epoch(orbp) .+ (0:10:12_000) ./ 86400

# Finite-difference Jacobian (default)
orb_fd, P_fd = fit_j2_mean_elements(vjd, vr_i, vv_i;
    mean_elements_epoch = vjd[begin],
    verbose = false,
)

# ForwardDiff Jacobian
orb_ad, P_ad = fit_j2_mean_elements(vjd, vr_i, vv_i;
    mean_elements_epoch = vjd[begin],
    jacobian_method     = ForwardDiffJacobian(),
    verbose = false,
)
```

## Performance Comparison

The table below shows representative median timings for the J2, J2 osculating, J4, and J4
osculating propagators on a LEO orbit (1201 data points, `atol = rtol = 2e-4`, default
settings).

| Propagator | FD Median | AD Median | Speedup |
|:-----------|----------:|----------:|--------:|
| J2         | ~2.5 ms   | ~1.8 ms   | ~1.4× AD |
| J2osc      | ~5.0 ms   | ~3.0 ms   | ~1.7× AD |
| J4         | ~3.0 ms   | ~2.0 ms   | ~1.5× AD |
| J4osc      | ~6.0 ms   | ~3.5 ms   | ~1.7× AD |

!!! note

    Exact numbers depend on hardware and Julia version. The key takeaway is that
    `ForwardDiffJacobian()` is typically faster for these propagators because it avoids
    `N+1` propagation evaluations per Jacobian (finite differences) and instead computes
    exact derivatives in a single forward pass with dual numbers.

### Running Your Own Benchmark

```julia
using SatelliteToolboxPropagators
using BenchmarkTools

orb_input = KeplerianElements(
    date_to_jd(2023, 1, 1),
    7130.982e3, 0.001111,
    98.405 |> deg2rad, 90.0 |> deg2rad,
    200.0  |> deg2rad, 45.0 |> deg2rad,
)

orbp = Propagators.init(Val(:J2), orb_input)
ret  = Propagators.propagate!.(orbp, 0:10:12_000)
vr_i = first.(ret)
vv_i = last.(ret)
vjd  = Propagators.epoch(orbp) .+ (0:10:12_000) ./ 86400

kw = (; mean_elements_epoch = vjd[begin], verbose = false)

println("FiniteDiffJacobian:")
@btime fit_j2_mean_elements($vjd, $vr_i, $vv_i;
    jacobian_method = FiniteDiffJacobian(), $kw...)

println("ForwardDiffJacobian:")
@btime fit_j2_mean_elements($vjd, $vr_i, $vv_i;
    jacobian_method = ForwardDiffJacobian(), $kw...)
```

## Accuracy

Both methods converge to the same fitted mean elements. The `ForwardDiffJacobian` computes
*exact* (to machine precision) derivatives, while `FiniteDiffJacobian` introduces a small
truncation error controlled by `jacobian_perturbation` and `jacobian_perturbation_tol`. In
practice, the difference in the final fitted elements is negligible for well-conditioned
problems. For ill-conditioned problems or very tight tolerances, the AD Jacobian can provide
better convergence behavior.

## Which Method to Choose?

- **`FiniteDiffJacobian()`** (default): 
  Good for most use cases. Allows fine-tuning via `jacobian_perturbation` and
  `jacobian_perturbation_tol`.
- **`ForwardDiffJacobian()`**: Recommended when performance matters, when fitting many
  TLEs/orbits in a batch, or when tighter convergence tolerances are needed. The exact
  Jacobian eliminates perturbation-related tuning entirely.
