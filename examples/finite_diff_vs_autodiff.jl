# ==========================================================================================
#              Finite-Difference vs ForwardDiff Jacobian Comparison
# ==========================================================================================
#
# Benchmarks FiniteDiffJacobian() vs ForwardDiffJacobian() for mean-element fitting across
# the J2, J2 osculating, J4, and J4 osculating propagators.
#
# Usage:
#   julia --project examples/finite_diff_vs_autodiff.jl
#
# ==========================================================================================

using SatelliteToolboxPropagators
using BenchmarkTools
using Printf
using Dates

const REPORT_PATH = joinpath(@__DIR__, "finite_diff_vs_autodiff_report.md")

const PROPAGATORS = [
    (name = "J2",     sym = Val(:J2)),
    (name = "J2osc",  sym = Val(:J2osc)),
    (name = "J4",     sym = Val(:J4)),
    (name = "J4osc",  sym = Val(:J4osc)),
]

# ------------------------------------------------------------------------------------------
# Orbit definition
# ------------------------------------------------------------------------------------------

const ORB_INPUT = KeplerianElements(
    date_to_jd(2023, 1, 1),
    7130.982e3,
    0.001111,
    98.405 |> deg2rad,
    90.0   |> deg2rad,
    200.0  |> deg2rad,
    45.0   |> deg2rad,
)

function generate_osc_data(sym, orb_input, Δt_range)
    orbp = Propagators.init(sym, orb_input)
    ret  = Propagators.propagate!.(orbp, Δt_range)
    vr_i = first.(ret)
    vv_i = last.(ret)
    vjd  = Propagators.epoch(orbp) .+ collect(Δt_range) ./ 86400
    return vjd, vr_i, vv_i
end

function run_propagator(prop, Δt_range)
    vjd, vr_i, vv_i = generate_osc_data(prop.sym, ORB_INPUT, Δt_range)
    kw = (; mean_elements_epoch = vjd[begin], verbose = false)

    println("\n  Benchmarking FiniteDiffJacobian...")
    b_fd = @benchmark Propagators.fit_mean_elements(
        $(prop.sym), $vjd, $vr_i, $vv_i;
        jacobian_method = FiniteDiffJacobian(), $kw...
    )
    display(b_fd)

    println("\n  Benchmarking ForwardDiffJacobian...")
    b_ad = @benchmark Propagators.fit_mean_elements(
        $(prop.sym), $vjd, $vr_i, $vv_i;
        jacobian_method = ForwardDiffJacobian(), $kw...
    )
    display(b_ad)

    t_fd     = median(b_fd).time / 1e6
    t_ad     = median(b_ad).time / 1e6
    alloc_fd = median(b_fd).allocs
    alloc_ad = median(b_ad).allocs
    mem_fd   = median(b_fd).memory / 1024
    mem_ad   = median(b_ad).memory / 1024

    @printf("\n  Median: FD = %.1f ms, AD = %.1f ms\n\n", t_fd, t_ad)

    return (;
        name = prop.name,
        t_fd, t_ad,
        alloc_fd, alloc_ad,
        mem_fd, mem_ad,
    )
end

# ------------------------------------------------------------------------------------------
# Warmup
# ------------------------------------------------------------------------------------------

println("Warmup pass: compiling all code paths...")
Δt_range = 0:10:12_000

for prop in PROPAGATORS
    vjd, vr_i, vv_i = generate_osc_data(prop.sym, ORB_INPUT, Δt_range)
    kw = (; mean_elements_epoch = vjd[begin], verbose = false)
    Propagators.fit_mean_elements(prop.sym, vjd, vr_i, vv_i; jacobian_method = FiniteDiffJacobian(),  kw...)
    Propagators.fit_mean_elements(prop.sym, vjd, vr_i, vv_i; jacobian_method = ForwardDiffJacobian(), kw...)
end

println("Warmup complete.\n")

# ------------------------------------------------------------------------------------------
# Run all propagators
# ------------------------------------------------------------------------------------------

results = []

for prop in PROPAGATORS
    println("=" ^ 80)
    println("Propagator: $(prop.name)")
    println("=" ^ 80)
    push!(results, run_propagator(prop, Δt_range))
end

# ------------------------------------------------------------------------------------------
# Generate Markdown report
# ------------------------------------------------------------------------------------------

open(REPORT_PATH, "w") do io
    println(io, "# Finite-Difference vs ForwardDiff Jacobian — Benchmark Report")
    println(io)
    println(io, "> Auto-generated on $(Dates.format(now(), "yyyy-mm-dd HH:MM:SS"))")
    println(io, "> Julia $(VERSION) — $(Sys.MACHINE)")
    println(io)

    println(io, "## Orbit")
    println(io)
    println(io, "| Parameter | Value |")
    println(io, "|:----------|------:|")
    @printf(io, "| a  | %.3f km |\n", ORB_INPUT.a / 1e3)
    @printf(io, "| e  | %.6f |\n", ORB_INPUT.e)
    @printf(io, "| i  | %.3f° |\n", rad2deg(ORB_INPUT.i))
    @printf(io, "| Ω  | %.1f° |\n", rad2deg(ORB_INPUT.Ω))
    @printf(io, "| ω  | %.1f° |\n", rad2deg(ORB_INPUT.ω))
    @printf(io, "| f  | %.1f° |\n", rad2deg(ORB_INPUT.f))
    println(io)
    @printf(io, "Data points: %d (Δt = 0:10:%d s)\n", length(Δt_range), last(Δt_range))
    println(io)

    println(io, "## Performance")
    println(io)
    println(io, "| Propagator | FD Median (ms) | AD Median (ms) | Speedup | FD Allocs | AD Allocs | FD Mem (KiB) | AD Mem (KiB) |")
    println(io, "|:-----------|---------------:|---------------:|--------:|----------:|----------:|-------------:|-------------:|")

    for r in results
        ratio = r.t_fd / r.t_ad
        arrow = ratio > 1 ? "AD" : "FD"
        @printf(io, "| %s | %.1f | %.1f | %.2f× %s | %d | %d | %.0f | %.0f |\n",
            r.name, r.t_fd, r.t_ad, max(ratio, 1/ratio), arrow,
            r.alloc_fd, r.alloc_ad, r.mem_fd, r.mem_ad)
    end
    println(io)
end

println("=" ^ 80)
println("Report written to: $REPORT_PATH")
println("=" ^ 80)
