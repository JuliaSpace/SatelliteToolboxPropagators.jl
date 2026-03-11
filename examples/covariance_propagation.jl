# ==========================================================================================
#              Covariance Propagation via J2 Propagator + ForwardDiff Jacobian
# ==========================================================================================
#
# Demonstrates linear covariance propagation through the J2 analytical orbit propagator
# using ForwardDiff to compute the state transition Jacobian.
#
# Given an initial covariance P₀ in Cartesian state space, the propagated covariance at
# elapsed time Δt is:
#
#   P(Δt) = J(Δt) P₀ J(Δt)ᵀ
#
# where J = ∂[r; v]/∂x₀ is the 6×6 Jacobian of the J2 propagator mapping from the
# initial Cartesian state to the propagated Cartesian state.
#
# Usage:
#   julia --project examples/covariance_propagation.jl
#
# ==========================================================================================

using SatelliteToolboxPropagators
using ForwardDiff
using LinearAlgebra
using StaticArrays
using Printf

# ------------------------------------------------------------------------------------------
# Orbit definition
# ------------------------------------------------------------------------------------------

orb_input = KeplerianElements(
    date_to_jd(2023, 1, 1),
    7130.982e3, 0.001111,
    98.405 |> deg2rad, 90.0 |> deg2rad,
    200.0  |> deg2rad, 45.0 |> deg2rad,
)

epoch = orb_input.t

# ------------------------------------------------------------------------------------------
# Define the J2 propagator map: x₀ = [r₀; v₀] → [r(Δt); v(Δt)]
# ------------------------------------------------------------------------------------------

function j2_map(x₀::AbstractVector, Δt, epoch_jd; j2c = j2c_egm2008)
    r₀ = SVector{3}(x₀[1], x₀[2], x₀[3])
    v₀ = SVector{3}(x₀[4], x₀[5], x₀[6])

    T = eltype(x₀)
    j2d = J2Propagator{typeof(epoch_jd), T}()
    j2d.j2c = J2PropagatorConstants{T}(T(j2c.R0), T(j2c.μm), T(j2c.J2))

    orb = rv_to_kepler(r₀, v₀, epoch_jd)
    j2_init!(j2d, orb)
    r, v = j2!(j2d, Δt)
    return vcat(r, v)
end

# ------------------------------------------------------------------------------------------
# Compute initial Cartesian state
# ------------------------------------------------------------------------------------------

orbp = Propagators.init(Val(:J2), orb_input)
r₀, v₀ = Propagators.propagate!(orbp, 0.0)
x₀ = Vector(vcat(SVector{3}(r₀), SVector{3}(v₀)))

# ------------------------------------------------------------------------------------------
# Fabricate an initial covariance P₀
# ------------------------------------------------------------------------------------------
#
# In practice, P₀ would come from an orbit determination solution. Here we use
# representative diagonal uncertainties for demonstration.

σ_pos = 50.0    # [m] position 1σ
σ_vel = 0.05    # [m/s] velocity 1σ

P₀ = Diagonal([
    σ_pos^2, σ_pos^2, σ_pos^2,
    σ_vel^2, σ_vel^2, σ_vel^2,
])

# ------------------------------------------------------------------------------------------
# Propagate covariance at several times
# ------------------------------------------------------------------------------------------

println("=" ^ 86)
println("    Covariance Propagation via J2 Propagator + ForwardDiff Jacobian")
println("=" ^ 86)
println()
@printf("  Orbit:  a = %.3f km, e = %.6f, i = %.3f°\n",
    orb_input.a / 1e3, orb_input.e, rad2deg(orb_input.i))
@printf("  Epoch:  JD %.8f\n", epoch)
@printf("  State:  x₀ = [r₀; v₀]  (Cartesian, inertial)\n")
@printf("  P₀:     diag(σ_pos² = %.0f m², σ_vel² = %.4f m²/s²)\n", σ_pos^2, σ_vel^2)
println()

println("── Propagated 1σ uncertainties ────────────────────────────────────────────────────")
println()
@printf("  %10s  %12s  %12s  %12s  %14s  %14s  %14s\n",
    "Δt [s]", "σ_x [m]", "σ_y [m]", "σ_z [m]",
    "σ_vx [m/s]", "σ_vy [m/s]", "σ_vz [m/s]")
@printf("  %10s  %12s  %12s  %12s  %14s  %14s  %14s\n",
    "──────────", "────────────", "────────────", "────────────",
    "──────────────", "──────────────", "──────────────")

for Δt in [0.0, 60.0, 300.0, 600.0, 1800.0, 3600.0, 7200.0, 14400.0, 43200.0, 86400.0]
    J = ForwardDiff.jacobian(x -> j2_map(x, Δt, epoch), x₀)

    P_t = J * P₀ * J'

    σ_r = sqrt.(diag(P_t)[1:3])
    σ_v = sqrt.(diag(P_t)[4:6])

    @printf("  %10.0f  %12.3f  %12.3f  %12.3f  %14.6e  %14.6e  %14.6e\n",
        Δt, σ_r[1], σ_r[2], σ_r[3], σ_v[1], σ_v[2], σ_v[3])
end

# ------------------------------------------------------------------------------------------
# Detailed view at Δt = 3600 s
# ------------------------------------------------------------------------------------------

Δt_detail = 3600.0

println()
println("── Detailed view at Δt = $(Int(Δt_detail)) s ──────────────────────────────────────")
println()

y = j2_map(x₀, Δt_detail, epoch)
J = ForwardDiff.jacobian(x -> j2_map(x, Δt_detail, epoch), x₀)

@printf("  Position:  [%+14.3f, %+14.3f, %+14.3f] m\n", y[1], y[2], y[3])
@printf("  Velocity:  [%+14.6f, %+14.6f, %+14.6f] m/s\n", y[4], y[5], y[6])
println()

println("  Jacobian J = ∂[r;v]/∂x₀  (6×6):")
for i in 1:6
    @printf("    [")
    for j in 1:6
        @printf(" %+13.5e", J[i, j])
        j < 6 && @printf(",")
    end
    @printf(" ]\n")
end

P_t = J * P₀ * J'

println()
println("  Propagated covariance P(Δt) (6×6):")
for i in 1:6
    @printf("    [")
    for j in 1:6
        @printf(" %+13.5e", P_t[i, j])
        j < 6 && @printf(",")
    end
    @printf(" ]\n")
end

σ_pos_rss = sqrt(tr(P_t[1:3, 1:3]))
σ_vel_rss = sqrt(tr(P_t[4:6, 4:6]))
println()
@printf("  RSS position uncertainty: %.3f m\n", σ_pos_rss)
@printf("  RSS velocity uncertainty: %.6e m/s\n", σ_vel_rss)

println()
println("=" ^ 86)
