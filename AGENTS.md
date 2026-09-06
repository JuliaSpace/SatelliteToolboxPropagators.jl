# Repository Instructions

## Package Structure

- This package provides analytical and SGP4/SDP4 orbit propagators for the SatelliteToolbox.jl ecosystem and supports Julia 1.10 or newer within Julia 1.x (`[compat] julia = "1.10, 1.11, 1.12"`). It depends on **SatelliteToolboxBase.jl** v2, whose `KeplerianElements{Tanomaly, Tepoch, T}` selects the stored anomaly.
- `src/SatelliteToolboxPropagators.jl` is the module entrypoint. It loads, in this order, the generic `Propagators` API (`src/api/Propagators.jl`), the shared types (`src/types.jl`), the propagator-specific API methods (`src/api/*.jl`), the numerical kernels (`src/propagators/*.jl`), and finally the precompile workloads (`src/precompile.jl`). Preserve this dependency order when adding code.
- `src/api/Propagators.jl` defines the common `OrbitPropagator` interface and the time and sink conversions shared by every propagator. Files under `src/api/` connect the concrete propagators to that interface, while files under `src/propagators/` contain the J2, J2-osculating, J4, J4-osculating, and two-body kernels, the J2 short-period corrections shared by the osculating propagators (`osculating.jl`), and the least-square algorithm shared by every mean elements fitting function (`fit.jl`). SGP4 delegates its kernel to `SatelliteToolboxSgp4` and therefore has no local file under `src/propagators/`.
- `test/runtests.jl` includes the feature tests `test/j2.jl`, `test/j2osc.jl`, `test/j4.jl`, `test/j4osc.jl`, `test/sgp4.jl`, `test/twobody.jl`, and `test/api.jl`, then, on stable Julia releases only, the quality checks in `test/quality.jl` (Aqua, JET) and the allocation checks in `test/performance.jl` (AllocCheck). Both tool groups are skipped on prereleases.
- Test-only dependencies (AllocCheck, Aqua, JET, Test) are declared in `[extras]` and `[targets]` of `Project.toml`; there is no `test/Project.toml`.
- `docs/` is a separate Documenter environment. User-facing guides live under `docs/src/man/`, the propagator design is described in `docs/src/man/API.md`, and the API reference is generated in `docs/src/lib/library.md`.
- Root and environment-specific manifest files are ignored by `.gitignore`; do not commit generated `Manifest.toml` files.

## Commands

- Instantiate the package environment: `julia --project=. -e 'using Pkg; Pkg.instantiate()'`.
- Load the package as a quick smoke test: `julia --project=. -e 'using SatelliteToolboxPropagators'`.
- Run the full test suite: `julia --project=. -e 'using Pkg; Pkg.test()'`.
- Run one feature test directly: `julia --project=. -e 'using Test, Dates, StaticArrays, SatelliteToolboxPropagators; include("test/twobody.jl")'`. Replace `twobody.jl` with any feature test file included by `test/runtests.jl`. `test/quality.jl` and `test/performance.jl` need the test dependencies and are only runnable through `Pkg.test()`.
- Bootstrap the docs environment after dependency changes or in a fresh checkout: `julia --project=docs -e 'using Pkg; Pkg.develop(path="."); Pkg.instantiate()'`.
- Build documentation locally: `julia --project=docs docs/make.jl local`.
- Format Julia sources with JuliaFormatter installed in the default environment: `julia -e 'using JuliaFormatter; format(".")'`. Do not add `--project=.` because JuliaFormatter is not a package dependency. Review the diff before committing: the formatter drops the alignment of the first line of an aligned block of one-line methods when a docstring precedes it, and that alignment is kept on purpose.
- Use generous timeouts for first runs. Instantiation, precompilation, and the full suite can be quiet for several minutes.

## Testing

- Run the nearest focused test while iterating, then run the full suite before finishing behavior changes. The allocation site limits of the fitting functions are only checked on Julia 1.10 and 1.11, so run the suite on one of them when changing the fitting or the Jacobian paths.
- Match the existing nested `@testset "..." verbose = true begin` organization and add regression tests alongside the affected propagator or in `test/api.jl` for shared API behavior.
- Preserve coverage for both `Float64` and `Float32`. Exercise `ForwardDiffJacobian` when changing fitting or Jacobian paths.
- Keep reference scenarios, their sources, units, and numerical tolerances visible in tests. Do not replace externally validated results with values generated only by the implementation under test.
- There is no test-name selector. Directly include a feature test for focused work; `test/quality.jl` and `test/performance.jl` depend on setup performed by `test/runtests.jl` and are not standalone focused tests.

## Code Style

- Follow BlueStyle and `.JuliaFormatter.toml`; its options are the formatting source of truth, except for the alignment case noted in Commands.
- Preserve the established four-space indentation, explicit `return` statements, trailing commas in multiline calls, aligned assignments where helpful, and section separators used throughout `src/` and `test/`. Private helpers use `function ... end` rather than the one-line form.
- Every function has a docstring, including the private ones, except the overloads of `Base` functions (`copy`, `convert`, `show`, `iterate`, ...), which keep the `Base` documentation and get a comment when their behavior needs explanation.
- Docstring signatures longer than 92 columns are broken like code definitions: one argument per line, no trailing comma after the last one, and `) where {...} -> Type` on the closing line. Several signatures of one docstring are separated by a blank line. Keyword defaults are written `(**Default**: value)` on their own line.
- Functions that return `KeplerianElements` state the anomaly type in the docstring return, e.g. `-> KeplerianElements{MeanAnomaly, Tepoch, T}`.
- Type names never start with an underscore, including the `const` unions used as type aliases (`PropagatorData`, `PropagationInstant`). Private functions and constants do start with an underscore.
- Exported constants use `SCREAMING_SNAKE_CASE`, e.g. `J2C_EGM2008`, `J4C_JGM03_F32`, and `TBC_M0`.
- Use descriptive domain notation already established in the codebase, including symbols such as `Δt`, `μ`, `Ω`, and `ω`, when it improves consistency with orbital mechanics formulas.
- Keep implementations generic over numeric and epoch types. Avoid introducing `Float64` conversions that would break `Float32`, `ForwardDiff.Dual`, or other `Number` subtypes unless an API explicitly requires them.
- Add or update docstrings for public API changes. State reference frames, timescales, units, mutation, keywords, and return shapes explicitly. Regenerate the `julia-repl` examples by executing them when the output changes.
- Keep exported API methods in the relevant `src/api/` file and numerical details in the matching `src/propagators/` file rather than bypassing the `Propagators` interface.

## Behavioral Constraints

- Preserve SI conventions: position is in meters, velocity in meters per second, propagation offsets are in seconds, angular elements are in radians, and numeric epochs are Julian Days. Epoch-oriented APIs interpret epochs as UTC.
- Distinguish propagation from stepping: `Propagators.propagate!` evaluates an offset from the initial epoch, while `Propagators.step!` advances relative to the current propagation instant.
- Propagators are mutable. Vector propagation copies state for parallel tasks and leaves the supplied propagator at the last requested instant; preserve this deterministic final state when changing the shared kernel `_propagate_vector!`.
- Preserve the sink API: tuple sinks return position and velocity separately, while `OrbitStateVector` sinks include the corresponding epoch.
- Every mean elements output (`Propagators.mean_elements`, the fitting functions, and the epoch updates) returns `KeplerianElements{MeanAnomaly}` for every propagator, including SGP4, whereas the SGP4 fitting functions return a `TLE`. The osculating propagators store their osculating elements as `KeplerianElements{TrueAnomaly}`.
- The propagator element type is chosen by `_propagator_eltype`: the constants type wins for floating-point elements, and the types are promoted otherwise, which keeps the propagation differentiable with `ForwardDiff.Dual` elements.
- Adding a propagator to the shared fitting algorithm only requires the dispatch methods at the top of `src/propagators/fit.jl` and a `_similar_propagator` method.
- Keep `Val` dispatch tags and concrete wrapper types consistent across initialization, fitting, propagation, copying, display, tests, docs, and precompile workloads.
- Treat allocation limits in `test/performance.jl` as part of expected behavior: the kernels and the Jacobians must not allocate, and the fitting functions must not add allocation sites outside their verbose printing path, whose decorated strings are rendered once at load time for that reason.

## Continuous Integration

- Standard CI tests Julia 1.10 (pinned minor) and the latest stable Julia 1.x on Linux x64, macOS arm64, and Windows x64; nightly CI exercises the same platform set separately.
- CI runs `julia-buildpkg` before `julia-runtest`, processes coverage on standard CI, and builds Documenter output on the latest stable Julia under Ubuntu.
- No formatter check is configured in CI. Run JuliaFormatter locally for Julia source changes.

## Not Configured

- There is no `deps/build.jl`; do not add a separate `Pkg.build()` step to local instructions.
- There are no package extensions, pre-commit hooks, or dedicated linter configuration.
- JuliaFormatter is not a project dependency and must not be added to `[deps]`. Aqua, JET, and AllocCheck are test dependencies only; do not add them to `[deps]`.
