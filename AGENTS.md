# Repository Instructions

## Package Structure

- This package provides analytical and SGP4/SDP4 orbit propagators for the SatelliteToolbox.jl ecosystem and supports Julia 1.10 or newer within Julia 1.x.
- `src/SatelliteToolboxPropagators.jl` is the module entrypoint. It loads the generic `Propagators` API, shared types, propagator-specific API methods, numerical kernels, and finally precompile workloads; preserve this dependency order when adding code.
- `src/api/Propagators.jl` defines the common `OrbitPropagator` interface. Files under `src/api/` connect concrete propagators to that interface, while files under `src/propagators/` contain the J2, J2-osculating, J4, J4-osculating, and two-body kernels. SGP4 delegates its kernel to `SatelliteToolboxSgp4` and therefore has no local file under `src/propagators/`.
- `test/runtests.jl` includes feature tests from `test/j2.jl`, `test/j2osc.jl`, `test/j4.jl`, `test/j4osc.jl`, `test/sgp4.jl`, `test/twobody.jl`, and `test/api.jl`, then runs the quality and allocation checks in `test/performance.jl` on stable Julia releases.
- `docs/` is a separate Documenter environment. User-facing guides live under `docs/src/man/`, and API references live in `docs/src/lib/library.md`.
- Root and environment-specific manifest files are ignored by `.gitignore`; do not commit generated `Manifest.toml` files.

## Commands

- Instantiate the package environment: `julia --project=. -e 'using Pkg; Pkg.instantiate()'`.
- Load the package as a quick smoke test: `julia --project=. -e 'using SatelliteToolboxPropagators'`.
- Run the full test suite: `julia --project=. -e 'using Pkg; Pkg.test()'`.
- Run one feature test directly: `julia --project=. -e 'using Test, Dates, StaticArrays, SatelliteToolboxPropagators; include("test/twobody.jl")'`. Replace `twobody.jl` with any non-performance test file included by `test/runtests.jl`.
- Bootstrap the docs environment after dependency changes or in a fresh checkout: `julia --project=docs -e 'using Pkg; Pkg.develop(path="."); Pkg.instantiate()'`.
- Build documentation locally: `julia --project=docs docs/make.jl local`.
- Format Julia sources with JuliaFormatter installed in the default environment: `julia -e 'using JuliaFormatter; format(".")'`. Do not add `--project=.` because JuliaFormatter is not a package dependency.
- Check that formatting produced no diff: run the formatter, then run `git diff --exit-code`.
- Use generous timeouts for first runs. Instantiation, precompilation, and the full suite can be quiet for several minutes.

## Testing

- Run the nearest focused test while iterating, then run the full suite before finishing behavior changes.
- The full stable-Julia test run dynamically adds Aqua, JET, and AllocCheck to the temporary test environment. JET and allocation checks are skipped on Julia 1.12 and newer where the repository records upstream incompatibilities; all performance checks are skipped on nightly.
- Match the existing nested `@testset "..." verbose = true begin` organization and add regression tests alongside the affected propagator or in `test/api.jl` for shared API behavior.
- Preserve coverage for both `Float64` and `Float32`. Exercise `ForwardDiffJacobian` when changing fitting or Jacobian paths.
- Keep reference scenarios, their sources, units, and numerical tolerances visible in tests. Do not replace externally validated results with values generated only by the implementation under test.
- There is no test-name selector. Directly include a feature test for focused work; `test/performance.jl` depends on setup performed by `test/runtests.jl` and is not a standalone focused test.

## Code Style

- Follow BlueStyle and `.JuliaFormatter.toml`; its options are the formatting source of truth.
- Preserve the established four-space indentation, explicit `return` statements, trailing commas in multiline calls, aligned assignments where helpful, and section separators used throughout `src/` and `test/`.
- Use descriptive domain notation already established in the codebase, including symbols such as `Δt`, `μ`, `Ω`, and `ω`, when it improves consistency with orbital mechanics formulas.
- Keep implementations generic over numeric and epoch types. Avoid introducing `Float64` conversions that would break `Float32`, `ForwardDiff.Dual`, or other `Number` subtypes unless an API explicitly requires them.
- Add or update docstrings for public API changes. State reference frames, timescales, units, mutation, keywords, and return shapes explicitly.
- Keep exported API methods in the relevant `src/api/` file and numerical details in the matching `src/propagators/` file rather than bypassing the `Propagators` interface.

## Behavioral Constraints

- Preserve SI conventions: position is in meters, velocity in meters per second, propagation offsets are in seconds, angular elements are in radians, and numeric epochs are Julian Days. Epoch-oriented APIs interpret epochs as UTC.
- Distinguish propagation from stepping: `Propagators.propagate!` evaluates an offset from the initial epoch, while `Propagators.step!` advances relative to the current propagation instant.
- Propagators are mutable. Vector propagation copies state for parallel tasks and leaves the supplied propagator at the last requested instant; preserve this deterministic final state when changing threading code.
- Preserve the sink API: tuple sinks return position and velocity separately, while `OrbitStateVector` sinks include the corresponding epoch.
- Keep `Val` dispatch tags and concrete wrapper types consistent across initialization, fitting, propagation, copying, display, tests, docs, and precompile workloads.
- Treat allocation limits in `test/performance.jl` as part of expected behavior for hot propagation and Jacobian paths.

## Continuous Integration

- Standard CI tests the oldest supported Julia minor (1.10) and latest stable Julia on Linux x64, macOS arm64, and Windows x64; nightly CI exercises the same platform set separately.
- CI runs `julia-buildpkg` before `julia-runtest`, processes coverage on standard CI, and builds Documenter output on latest stable Julia under Ubuntu.
- No formatter check is configured in CI. Run JuliaFormatter locally for Julia source changes.

## Not Configured

- There is no `deps/build.jl`; do not add a separate `Pkg.build()` step to local instructions.
- There are no package extensions, pre-commit hooks, or dedicated linter configuration.
- JuliaFormatter, Aqua, JET, and AllocCheck are not regular project dependencies. Do not add them to `[deps]` merely to run developer checks.
