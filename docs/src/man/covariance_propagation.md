# Covariance Propagation

```@meta
CurrentModule = SatelliteToolboxPropagators
```

This page demonstrates **linear covariance propagation** through the J2 analytical orbit
propagator using `ForwardDiff.jl` to compute the state transition Jacobian. The same
approach applies to J4, J2 osculating, and J4 osculating propagators.

## Background

Given an initial covariance \(P_0\) in Cartesian state space, the propagated covariance at
elapsed time \(\Delta t\) is:

```math
P(\Delta t) = J(\Delta t) \, P_0 \, J(\Delta t)^\top
```

where \(J = \partial [r; v] / \partial x_0\) is the \(6 \times 6\) Jacobian of the
propagator mapping from the initial state to the propagated state. `ForwardDiff.jl` computes
this Jacobian exactly (to machine precision) in a single forward pass.

## Step-by-Step Example

### 1. Set Up the Orbit

```julia
using SatelliteToolboxPropagators
using ForwardDiff
using LinearAlgebra
using StaticArrays
using Printf

orb_input = KeplerianElements(
    date_to_jd(2023, 1, 1),
    7130.982e3, 0.001111,
    98.405 |> deg2rad, 90.0 |> deg2rad,
    200.0  |> deg2rad, 45.0 |> deg2rad,
)
```

### 2. Define the Propagator Map

We define a function that maps the initial Cartesian state `x₀ = [r₀; v₀]` to the
propagated state `[r(Δt); v(Δt)]` through the J2 propagator:

```julia
function j2_map(x₀::AbstractVector, Δt, epoch; j2c = j2c_egm2008)
    r₀ = SVector{3}(x₀[1], x₀[2], x₀[3])
    v₀ = SVector{3}(x₀[4], x₀[5], x₀[6])

    T = eltype(x₀)
    j2d = J2Propagator{typeof(epoch), T}()
    j2d.j2c = J2PropagatorConstants{T}(T(j2c.R0), T(j2c.μm), T(j2c.J2))

    orb = rv_to_kepler(r₀, v₀, epoch)
    j2_init!(j2d, orb)
    r, v = j2!(j2d, Δt)
    return vcat(r, v)
end
```

### 3. Compute the Jacobian with ForwardDiff

```julia
epoch = date_to_jd(2023, 1, 1)

orbp = Propagators.init(Val(:J2), orb_input)
r₀, v₀ = Propagators.propagate!(orbp, 0.0)
x₀ = vcat(SVector{3}(r₀), SVector{3}(v₀))

Δt = 3600.0  # 1 hour

J = ForwardDiff.jacobian(x -> j2_map(x, Δt, epoch), x₀)
```

### 4. Fabricate an Initial Covariance

In practice, \(P_0\) comes from an orbit determination solution (least-squares fit, Kalman
filter, etc.). Here we use representative diagonal uncertainties for demonstration:

```julia
σ_pos = 50.0    # [m] position 1σ
σ_vel = 0.05    # [m/s] velocity 1σ

P₀ = Diagonal(SVector{6}(
    σ_pos^2, σ_pos^2, σ_pos^2,
    σ_vel^2, σ_vel^2, σ_vel^2,
))
```

### 5. Propagate the Covariance

```julia
P_Δt = J * P₀ * J'

σ_r = sqrt.(diag(P_Δt)[1:3])
σ_v = sqrt.(diag(P_Δt)[4:6])

@printf("Position 1σ at Δt = %.0f s:\n", Δt)
@printf("  σ_x = %.3f m, σ_y = %.3f m, σ_z = %.3f m\n", σ_r...)
@printf("Velocity 1σ at Δt = %.0f s:\n", Δt)
@printf("  σ_vx = %.6f m/s, σ_vy = %.6f m/s, σ_vz = %.6f m/s\n", σ_v...)
```

### 6. Covariance Time History

We can propagate the covariance at multiple times to observe how uncertainty grows:

```julia
println("─" ^ 78)
@printf("  %8s  %12s  %12s  %12s  %12s  %12s  %12s\n",
    "Δt [s]", "σ_x [m]", "σ_y [m]", "σ_z [m]",
    "σ_vx [m/s]", "σ_vy [m/s]", "σ_vz [m/s]")
println("─" ^ 78)

for Δt in [0.0, 60.0, 300.0, 600.0, 1800.0, 3600.0, 7200.0, 14400.0, 43200.0, 86400.0]
    J_t = ForwardDiff.jacobian(x -> j2_map(x, Δt, epoch), x₀)

    P_t = J_t * P₀ * J_t'

    σ_r = sqrt.(diag(P_t)[1:3])
    σ_v = sqrt.(diag(P_t)[4:6])

    @printf("  %8.0f  %12.3f  %12.3f  %12.3f  %12.6f  %12.6f  %12.6f\n",
        Δt, σ_r..., σ_v...)
end
```

## Notes

- **Linear approximation**: The formula \(P(\Delta t) = J P_0 J^\top\) is valid when the
  uncertainties are small enough that the propagator can be linearized around the nominal
  trajectory. For large uncertainties or long propagation times, consider sigma-point or
  Monte Carlo methods.

- **Mean-element covariance**: The `fit_j2_mean_elements` function with
  `ForwardDiffJacobian()` returns the covariance of the *fit residuals* in mean-element
  space, which is a different quantity from the propagated Cartesian covariance shown here.

- **Other propagators**: Replace `J2Propagator`, `J2PropagatorConstants`, `j2_init!`, and
  `j2!` with their J4 or osculating counterparts to propagate covariance through those
  models. The structure is identical.
